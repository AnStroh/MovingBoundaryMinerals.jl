```@meta
CurrentModule = MovingBoundaryMinerals
```
# [The (chemical) Stefan problem](@id stefan-problem)

The Stefan problem describes the movement of a reaction front in a thermodynamically constrained problem (e.g., the propagation of an ice front or the crystallization/resorption of minerals) [Rubinstein1971](@cite). We included the option to calculate the growth/dissolution of crystals based on a thermodynamically constrained data set (phase diagram). Here, we want to include a short description and additional information on the code. For the full derivation and benchmarks, see [Stroh2025](@cite).

## Input parameters

Most of the input parameters are self-explanatory. However, in this section we explain some more details on the individual input parameters. The following input parameters are required: the initial position of the interface, which also determines the size of the left phase in the composition profile, and the composition of the assemblage of interest `CompInt`. For the calculations, information on the start and end temperature (`Tstart`, `Tstop`) of the growth process as well as maximum and minimum temperatures (`TMIN`, `TMAX`) of the chosen section of the phase diagram are required. `TMIN` and `TMAX` are calculation parameters and define a temperature range in which the model can be used. They are not to be confused with `Tstart` and `Tstop` in the model, whereby `Tstart, Tstop` ``\in`` `[TMIN, TMAX]`. All temperatures must be given in K. Coefficients for the phase transformation lines can be stored in the parameter `eq_values`. `CompInt` stores the composition of interest of the solid solution and therefore determines the start composition.

## Determining the interface composition

By digitizing two adjacent reaction lines of the phase diagram, we create a binary phase diagram to which we can apply our code (see [Digitization](@ref digitization)). The reaction lines are described by two second-degree polynomials:

```math
X(T) = a T^2 + b T + c \qquad (1)
```

We use these polynomials ([`composition`](@ref)) to determine the compositions of the two phases as a function of temperature. These compositions are linked to the composition of the solid solution/assemblage via the lever rule. The initial composition profile is a step function with two homogeneous parts referring to the two materials. The respective compositions are based on the data from the phase diagram for the composition of `CompInt` at the starting temperature `Tstart`. Using this information, we calculate the length of the whole modelling domain.

## Numerical approach

The moving interface is tracked explicitly rather than fixed in space. At every time step:

1. Equilibrium compositions at the interface are evaluated from the phase diagram (or a fixed partition coefficient `K_D`) at the current temperature.
2. The interface velocity follows from the mass-balance (Stefan) condition, which equates the jump in diffusive flux across the interface to the rate of mass consumed/released by the moving boundary — see [`set_inner_bc_stefan!`](@ref), [`set_inner_bc_mb!`](@ref), and [`set_inner_bc_flux!`](@ref) for the different ways this condition can be discretized.
3. The interface is advected with this velocity and the grid on both sides is regridded around the new interface position ([`advect_interface_regrid!`](@ref), [`regrid!`](@ref)) so that resolution stays concentrated near the steep compositional gradient at the interface.
4. Diffusion in both phases is then solved implicitly on the new grid with a Galerkin finite-element discretization ([`construct_matrix_fem`](@ref), [`solve_soe`](@ref)).

This avoids prescribing the growth/resorption rate a priori: it is a direct consequence of the diffusive fluxes and the thermodynamic interface condition.
