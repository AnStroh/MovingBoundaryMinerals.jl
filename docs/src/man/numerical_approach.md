```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Numerical Approach](@id numerical-approach)

The moving interface is tracked explicitly rather than fixed in space. At every time step:

1. **Interface condition** — applied via [`set_inner_bc_flux!`](@ref), [`set_inner_bc_mb!`](@ref), or [`set_inner_bc_stefan!`](@ref), depending on the problem family.
2. **Advection and regridding** — the interface is advected with the resulting velocity and the grid on both sides is rebuilt around its new position ([`advect_interface_regrid!`](@ref), [`regrid!`](@ref)).
3. **Diffusion solve** — diffusion in both phases is solved implicitly with a Galerkin finite-element scheme ([`construct_matrix_fem`](@ref), [`solve_soe`](@ref)).

For the full derivation and benchmarks, see [Stroh2025](@cite).
