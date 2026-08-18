```@meta
CurrentModule = MovingBoundaryMinerals
```
# [The (chemical) Stefan problem](@id stefan-problem)

The Stefan problem describes the movement of a reaction front in a thermodynamically constrained system — the propagation of an ice front, or the crystallization/resorption of minerals, are classic examples [Rubinstein1971](@cite). `MovingBoundaryMinerals.jl` solves this directly for mineral growth and resorption, using a digitized phase diagram to constrain the interface composition rather than a fixed partition coefficient. This page covers the physics and the input parameters specific to that mode; for the full derivation and benchmarks, see [Stroh2025](@cite).

## Input parameters

This section covers the input parameters that need more explanation than their names alone provide. Required inputs are: the initial interface position (which also sets the size of the left phase in the composition profile); the composition of interest `CompInt`, which fixes the starting composition; the start and end temperature of the growth process (`Tstart`, `Tstop`); and the temperature bounds of the digitized phase-diagram section (`TMIN`, `TMAX`), which must satisfy `Tstart, Tstop` ``\in`` `[TMIN, TMAX]` — don't confuse the two pairs. All temperatures are given in K. The phase-transformation-line coefficients go in `eq_values`.

## Determining the interface composition

Digitizing two adjacent reaction lines of the phase diagram reduces it to the binary, two-phase system our code solves for (see [Digitization](@ref digitization)). Each line is fit as a quadratic function of temperature:

```math
X(T) = a + b T + c T^2 \qquad (1)
```

We use these polynomials ([`composition`](@ref)) to determine the compositions of the two phases as a function of temperature, linked to the composition of the solid solution/assemblage via the lever rule. The initial composition profile is a step function — two homogeneous regions, one per phase — with each region's composition read from the phase diagram at `CompInt` and the starting temperature `Tstart`; this also fixes the length of the modelling domain.

## Numerical approach

Equilibrium compositions at the interface are evaluated from the phase diagram (or a fixed partition coefficient `K_D`) at the current temperature and applied as a Dirichlet condition ([`set_inner_bc_stefan!`](@ref)). This is the only part of the scheme specific to the Stefan problem; for the shared regridding and diffusion-solve steps, see [Numerical Approach](@ref numerical-approach).
