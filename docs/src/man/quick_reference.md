```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Quick reference](@id quick-reference)

A compact lookup for the input parameters that recur across almost every example, for when you already know the package and just need a reminder. For the full explanation of any of these, see [General Remarks](@ref) and [The Stefan Problem](@ref stefan-problem).

## Geometry and domain

| Parameter | Meaning |
|:--|:--|
| `n` | Geometry exponent: `1` = planar, `2` = cylindrical, `3` = spherical. |
| `Ri` | Radii `[interface_position, total_domain_length]` in \[m\] (`Ri[1]` is where phase A ends and phase B begins). |
| `CompInt` | Bulk composition of interest of the solid solution/assemblage; sets the starting composition (D-family/Stefan examples). |

## Boundary conditions

See [Boundary Conditions](@ref boundary-conditions) for the equations behind each of these.

| Parameter | Meaning |
|:--|:--|
| `BCout` | Outer BC `[left right]`: `0` = Neumann, `1` = Dirichlet. |
| Inner (interface) BC | Set by which function you call, not a flag: [`set_inner_bc_flux!`](@ref) (flux-balance condition, "B" examples), [`set_inner_bc_mb!`](@ref) (total mass-balance condition, "C" examples), [`set_inner_bc_stefan!`](@ref) (thermodynamically constrained Stefan condition, "D" examples). |

## Grid resolution and refinement

See [Mesh Refinement](@ref mesh-refinement) for the equations behind each `RefineMethod`.

| Parameter | Meaning |
|:--|:--|
| `res` / `resmin` | Number of grid nodes `[left right]` / minimum allowed after regridding. |
| `RefineMethod` | `1` = m-refinement (smooth geometric grading, via [`create_grid!`](@ref)); `2` = h-refinement stopped after a fixed number of levels (via [`h_refinement1`](@ref)); `3` = h-refinement stopped by a resolution condition (via [`h_refinement2`](@ref)). |
| `MRefin` | Refinement factor — only used when `RefineMethod == 1`. |
| `RefineLevel` | Number of refinement levels — only used when `RefineMethod == 2`. |
| `RefineCond` | Stops refining once the last `dx` on the left ``\le`` `RefineCond * Ri[1]` — only used when `RefineMethod == 3`. |
| `nPoints` | Initial number of grid points for h-refinement (`RefineMethod` `2` or `3`); a single number applies to both sides, or `[left right]` for different resolutions per side. |
| `CFL` | Courant–Friedrichs–Lewy number controlling the adaptive time step (via [`calculate_dt`](@ref)/[`find_dt`](@ref)). |

## Diffusion coefficient

| Parameter | Meaning |
|:--|:--|
| `Di` | Diffusion coefficient(s) `[left right]` in [m``^2``/s]. Set both entries to `-1.0` to compute D from the Arrhenius relationship instead of a constant. |
| `D0`, `Ea` | Arrhenius pre-exponential factor [m``^2``/s] and activation energy [J/mol], used when `Di = [-1.0 -1.0]`. |

## Temperature and phase equilibrium (D-family / Stefan examples)

| Parameter | Meaning |
|:--|:--|
| `Tstart`, `Tstop` | Start/end temperature of the growth process, in [K]. |
| `TMIN`, `TMAX` | Bounds of the temperature range covered by the phase-diagram data; must satisfy `Tstart, Tstop` ``\in`` `[TMIN, TMAX]`. Not the same as `Tstart`/`Tstop` — see [The Stefan Problem](@ref stefan-problem). |
| `eq_values` / `coeff_up`, `coeff_do` | Quadratic phase-transformation-line coefficients (a, b, c per line); see [Digitization](@ref digitization) and [`coeff_trans_line`](@ref)/[`composition`](@ref). |
| `K_D` | Distribution/partition coefficient; can be a constant or a time-varying vector (first entry = initial value, last = value at `t_tot`) — see [General Remarks](@ref). |

## Common diagnostics

See [Interpreting Output](@ref interpreting-output) for `MB_Error`, `Residual`, `KD_sim`/`T_sim`, and common setup-error messages.
