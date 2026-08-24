```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Configuration Options](@id configuration-options)

Beyond the physical parameters covered in [Quick Reference](@ref quick-reference), several inputs act as switches that change which code path the package actually takes. This page walks through each one, what it changes, and when to pick which — with pointers to the pages that cover the underlying equations in full.

## Interface (inner) boundary condition

Set by which function you call, not a flag: [`set_inner_bc_flux!`](@ref) (prescribed-rate growth/resorption via the flux-balance condition), [`set_inner_bc_mb!`](@ref) (the total mass-balance condition), or [`set_inner_bc_stefan!`](@ref) (the thermodynamically constrained Stefan condition, via a digitized phase diagram). See [Boundary Conditions](@ref boundary-conditions) for the equations behind each and which example family uses which.

## Outer boundary conditions (`BCout`)

`BCout = [left, right]`, each entry `0` (Neumann, zero flux) or `1` (Dirichlet, fixed composition) — see [Boundary Conditions](@ref boundary-conditions) for the equations.

!!! warning "Cylindrical/spherical geometry requires a Neumann left boundary"
    For `n = 2` (cylindrical) or `n = 3` (spherical), `BCout[1]` must be `0`: the left phase sits at the centre of the geometry (`x = 0`) and grows/resorbs from the inside as the right phase wraps around it, so a Dirichlet condition there is not physically meaningful. The templates in `main_codes/` check for this and raise an error otherwise.

## Geometry (`n`)

`n = 1` planar (Cartesian slab), `n = 2` cylindrical, `n = 3` spherical — sets the geometry exponent in the governing diffusion equation (see [Numerical Approach](@ref numerical-approach)) and interacts with the outer-BC restriction above.

## Diffusion coefficient: constant or Arrhenius (`Di`)

`Di = [D_left, D_right]` as fixed values, or `Di = [-1.0, -1.0]` to compute both from the Arrhenius relationship instead, using `D0` (pre-exponential factor) and `Ea` (activation energy) — see [Numerical Approach](@ref numerical-approach) for the equation. The same switch applies to single-crystal diffusion. `Chemical_Stefan_problem_XT.jl` is the one exception: diffusivities there are always constant (see [General Remarks](@ref)).

## Time-varying vs. constant `K_D`/temperature

`K_D` and temperature can each be given as a vector describing their evolution over the run, rather than a single fixed value — the same mechanism covers isothermal/non-isothermal and constant/time-varying `K_D` problems. The first entry is the initial value, the last is the value at `t_tot`; equal first/last entries hold the parameter constant throughout (see [General Remarks](@ref)).

## Grid refinement strategy (`RefineMethod`)

`RefineMethod = 1` (smooth geometric grading), `2` (h-refinement, fixed number of levels), or `3` (h-refinement, stopped by a resolution condition). See [Mesh Refinement](@ref mesh-refinement) for the full equations and guidance on which to pick.

## Units and non-dimensionalization (`scaling`/`rescale`)

You can work directly in SI units throughout, or non-dimensionalize before solving and rescale back afterwards via [`scaling`](@ref)/[`rescale`](@ref) — most templates in `main_codes/` and `examples/` do the latter. `scaling` fixes a length scale `Lsc = 1e-3` m and a diffusion scale `Dsc` (the average of `Di`, or of `D0` when using the Arrhenius relationship), derives the dependent time and velocity scales from them, and non-dimensionalizes `Ri`, `Di`, `D0`, `V_ip`, `t_tot`, and `t_ar` accordingly; `rescale` inverts this once the time loop finishes. This is why many function docstrings note "units may differ from SI units if non-dimensionalisation has been performed" — those functions operate on whatever units they're handed, dimensional or not, so the choice is yours. Non-dimensionalizing mainly helps numerically when the physical values involved (grid spacing, diffusivities) span many orders of magnitude.

## Verbosity (`verbose`)

A `Bool` accepted by the grid-construction and regridding functions ([`create_grid!`](@ref), [`define_new_grid`](@ref), [`regrid!`](@ref), [`newton_solver`](@ref)). Off by default; turn it on to print grid statistics and Newton-solver convergence information — useful when debugging a grading factor (`MRefin`) that fails to converge, or an unexpectedly coarse/fine interface cell.
