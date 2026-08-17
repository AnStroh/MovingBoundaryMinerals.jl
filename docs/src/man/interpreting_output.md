```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Interpreting output](@id interpreting-output)

Every example prints and/or returns a similar set of diagnostics. This page explains what each one *is* and how to reason about it — it is not a substitute for judging your own results.

!!! note "There is no universal pass/fail threshold"
    What counts as an acceptable value for `MB_Error`, `Residual`, etc. depends on your specific problem: whether you work in SI or non-dimensional units (see [`scaling`](@ref)/[`rescale`](@ref)), the geometry and scale of your domain, whether the system is open or closed, and how large the quantities involved are to begin with. The numbers below (e.g. "≪ 1%", "machine precision") are typical for well-resolved SI-unit examples in this package — they are starting points for building intuition, not general rules. Always interpret them in the context of your own setup, and cross-check against the composition profile plot rather than a single number in isolation.

## "The total mass difference is X% (relevant for closed systems)"

This message comes from [`calc_mass_err`](@ref), which computes

```math
\mathrm{ErrM} = \frac{\mathrm{Mass_{end}} - \mathrm{Mass_0}}{\mathrm{Mass_0}} \times 100\%
```

i.e. the *relative* change in total mass between the start and the current state, via [`calc_mass_vol`](@ref)/[`calc_mass_vol_simple_diff`](@ref). It is only meaningful for closed systems (no mass entering/leaving the modelling domain, e.g. Neumann outer BCs).

!!! warning "Large percentages can be a scaling artefact, not a bug"
    Because this is a *relative* error, it blows up whenever the reference mass `Mass0` (the initial mass) is close to zero in absolute terms — independent of whether anything is actually wrong. This happens, for example, in examples like A2 where the initial composition is homogeneously zero except at a single boundary node: `Mass0` is then tiny, so even a small absolute change in mass produces a huge percentage. If you see a large value, check the absolute mass values (`Mass`, `Mass0`) and the composition profile plot before assuming there is a conservation problem.

For problems with a substantial, non-trivial initial mass (most diffusion-couple examples in SI units), this percentage staying small (e.g. ≪ 1%, though what counts as "small enough" is problem-dependent) over the run is a reasonable sanity check that the time step and grid resolution are adequate — not proof that the result is physically correct.

## `Residual`

Examples with a moving interface push `(V_ip*dC - (JR - JL))^2` into a `Residual` array at every time step — a self-consistency check on the interface-velocity equation used to enforce the Stefan condition (see [The Stefan Problem](@ref stefan-problem)). In SI-unit examples this typically stays near machine precision (``\sim 10^{-20}`` or smaller); in a non-dimensional or rescaled setup the natural magnitude will differ. What matters more than the absolute value is the trend: a `Residual` that grows over the course of a run usually indicates the time step is too large or the grid near the interface is too coarse, regardless of the units used.

## Other commonly returned diagnostics

- `KD_sim`, `T_sim`: the simulated distribution coefficient and temperature at the interface at each stored time step — useful for plotting the simulated K_D(T) path against the thermodynamic K_D(T) curve (`KDlin`, `Tlin`).
- `C_left_check`, `C_right_check`, `T_check`: interface compositions and temperature recorded at every iteration, for debugging/plotting the interface history at a finer cadence than the other diagnostics.
- `Mass`, `Mass2`: total mass computed two different ways (via [`calc_mass_vol`](@ref) and via [`trapezoidal_integration`](@ref) directly) — they should agree; a persistent mismatch between them (rather than just a large relative-to-`Mass0` percentage) is worth investigating.

## Common setup errors

- `"Please change the size of the system. Increase Ri[2]."` / `"Please change the resolution of the system. res[2] >= res[1]."`: these are deliberate input checks (interface radius must be inside the domain; the right-phase resolution must not be lower than the left-phase resolution), not solver failures — adjust `Ri`/`res` accordingly.
- Oscillating or negative compositions: usually the CFL number is too high for the resolution, or the grid is too coarse at the interface. Lower `CFL` and/or increase resolution/refinement (`MRefin`, `RefineLevel`, `RefineCond`, `nPoints` depending on `RefineMethod`) before suspecting a physics/setup error.