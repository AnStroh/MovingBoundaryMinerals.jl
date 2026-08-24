```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Getting started](@id getting-started)

This page walks through running your first simulation with `MovingBoundaryMinerals.jl` and reading its output. We use [`examples/Simple_Diff.jl`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/Simple_Diff.jl) — the simplest example in the package — before pointing you to the more advanced diffusion-couple and moving-boundary examples.

!!! tip "Would rather not write any code at all?"
    See the [GUI](@ref gui) instead — a local, browser-based interface covering three of the package's core models with pre-filled forms, no Julia code required.

## Running the example

`Simple_Diff.jl` has no external file dependencies, so it can be run from anywhere once the package is installed (see [Installation](@ref)):

```julia-repl
julia> using MovingBoundaryMinerals
julia> include(joinpath(pkgdir(MovingBoundaryMinerals), "examples", "Simple_Diff.jl"))
```

!!! warning "Examples that read phase-diagram data"
    Examples [D1](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D1.jl) and [D2](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D2.jl) load phase-diagram coefficients from CSV files using paths relative to the repository root (e.g. `"examples/Examples_phase_diagram/Coefficients_Reaction_lines.csv"`). To run those, clone the repository and run Julia with the repository root as your working directory, e.g. `julia --project=. examples/D1.jl` from a terminal opened at the repository root. All other examples, including this one, are self-contained.

This solves diffusion of a single component into a homogeneous, semi-infinite planar crystal (no moving boundary, no second phase) — the simplest problem the package can solve, and a good place to see the pieces that every example shares.

## Anatomy of an example

Every example in `examples/` follows the same overall structure, whether it models a single crystal or a full moving-boundary diffusion couple:

1. **Physics** — material and process parameters: diffusivities (or Arrhenius parameters to compute them), domain size `L`, geometry exponent `n` (1 = planar, 2 = cylindrical, 3 = spherical), initial/boundary compositions, and a temperature-time history.
2. **Numerics** — grid resolution, a CFL number controlling the (adaptive) time step, and outer boundary condition types (`BCout`: `0` = Neumann, `1` = Dirichlet).
3. **Initial condition** — the starting composition profile and grid.
4. **Time loop** — repeated until `t_tot` is reached:
   - compute a stable time step ([`calculate_dt`](@ref) or a similar CFL-based expression),
   - advance time ([`update_time!`](@ref)),
   - update any time-dependent parameters, e.g. temperature and diffusivity, by interpolating along the user-supplied history ([`update_t_dependent_param_simple!`](@ref), or [`linear_interpolation_1D`](@ref) directly),
   - assemble and solve the finite-element system for this time step ([`fill_matrix!`](@ref), [`set_outer_bc!`](@ref)),
   - optionally plot or record diagnostics (mass balance, composition, etc.).
5. **Post-processing** — the final state is returned and, if enabled, plotted.

Examples that include a moving interface (letters B, C, D — see [List of examples](@ref)) add an extra step per iteration: the interface condition is applied (flux-balance, total mass-balance, or the thermodynamically constrained Stefan condition, depending on the family — see [Boundary Conditions](@ref boundary-conditions) for how these differ, including which ones prescribe the interface velocity `V_ip` and which one solves for it), the interface is advected, and the grid is regridded around its new position before the diffusion step. See [Numerical Approach](@ref numerical-approach) for details.

### Temperature (and `K_D`) as a history, not a single number

Note the `T_ar`/`t_ar` pair (or `Tpath`/`tpath`, `T_user`/`t_user` in other examples): temperature is not a single fixed value but a piecewise-linear history in time, evaluated with [`linear_interpolation_1D`](@ref) at every step. This history does not need to be monotonic — see example [D2](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D2.jl) for a temperature path that heats up and cools down between user-defined nodes. The same pattern applies to the distribution coefficient `K_D` where relevant. See [General Remarks](@ref) for the conventions this follows across examples.

## Reading the output

With `plot_sim = true` (or `plot_end = true`), `Simple_Diff.jl` plots the composition profile `C` against position `x`, together with the initial profile `C0`/`x0` for comparison. The function also returns `x, C, x0, C0, D, t, t_tot`, so you can inspect or replot the result yourself, e.g.:

```julia
using Plots
plot(x, C, label="final")
plot!(x0, C0, label="initial")
```

Every example also tracks a mass-balance error (`MB_Error` or similar) internally via [`calc_mass_err`](@ref) as a sanity check, along with a few other recurring diagnostics (`Residual`, `KD_sim`, `T_sim`, ...) — see [Interpreting Output](@ref interpreting-output) for what these mean and, importantly, when a large-looking mass-error percentage is actually expected rather than a sign of a problem.

## Where to go next

The examples increase in complexity roughly in this order (see [List of examples](@ref) for the full picture and naming convention):

1. `Simple_Diff` — single-crystal diffusion (this page).
2. `Diff_couple_no_interaction` — two single crystals side by side, no interface reaction.
3. `Diff_couple_Flux` / `Diff_couple_MB` — a true diffusion couple with a moving interface, using the local flux-balance or total mass-balance interface condition respectively (see [Boundary Conditions](@ref boundary-conditions)).
4. `Diff_couple_Flux_growth` / `Diff_couple_MB_growth` — the same, with simultaneous interface growth.
5. `A1`–`C2` — benchmarked versions of the above against analytical/semi-analytical solutions.
6. `D1`/`D2` — thermodynamically constrained crystal growth/resorption, where interface compositions come from a digitized phase diagram rather than a fixed `K_D` (see [Digitization](@ref digitization)).
