```@meta
CurrentModule = MovingBoundaryMinerals
```
# General Remarks

In this section we would like to make a few generally useful comments on the structure of the codes.

## Structure
We've kept the individual files as uniform as possible and stripped unnecessary parameters from the examples, so each one is quick to read and understand. The code works throughout in SI units — converting your own inputs is left to you, and the package does not check or convert units on your behalf. Most codes handle diffusion couples and consistently refer to the left (phase A) and right (phase B) side: wherever a variable stores two numbers, the first is always the left value and the second the right. See below for how the reusable `main_codes` templates relate to the worked examples.

## The `main_codes` templates

The [`main_codes`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/tree/main/main_codes) folder holds four bare-bones templates, one per topic, meant as a starting point for your own model rather than a worked example. Each corresponds to one of the [example](@ref "List of examples") families:

| Template | Problem family | Corresponding examples |
|:--|:--|:--|
| `Simple_Diff_1D_general.jl` | Single-crystal diffusion, no interface | `Simple_Diff`, A1, A2 |
| `Diff_couple_1D_general_Flux.jl` | Diffusion couple, flux-balance interface condition | `Diff_couple_Flux[_growth]`, B1–B7 |
| `Diff_couple_1D_general_MB.jl` | Diffusion couple, total mass-balance interface condition | `Diff_couple_MB[_growth]`, C1–C2 |
| `Chemical_Stefan_problem_XT.jl` | Thermodynamically constrained (Stefan) growth/resorption | D1, D2 |

## Compositions and concentrations
By default, concentrations are specified in mol; however, in the case where volume changes are negliable, density can be assumed to be unity (details in appenix A in [Stroh et al. (2025)](https://gmd.copernicus.org/articles/18/10203/2025/)), compositions/concentrations are equal. Therefore, model can handle either wt% or mol%. Whichever unit you choose, related parameters (activation energy, densities, ...) must be corrected accordingly — consistency across the whole model matters, since the units of concentration or composition only cancel out of the diffusion equation (Eq. 1) if applied uniformly. See Appendices A and B of [Stroh2025](@cite) for the full derivation.

## Temperature and distribution coefficient
Temperature and the distribution coefficient (`K_D`) can both be specified as vectors describing their evolution over time, rather than as single fixed values — this covers isothermal and non-isothermal problems, and constant or time-varying `K_D`, through the same mechanism. In each vector, the first entry is the initial value and the last is the value at the end of the simulation; if both are equal, the parameter is held constant throughout.

## Inner and outer boundaries
The interface (inner) boundary condition is set by which function you call rather than a flag — flux balance, total mass balance, or the thermodynamically constrained Stefan condition (see [Quick Reference](@ref quick-reference) for all three). Outer boundary conditions at the edges of the modelling domain are set independently via `BCout`, to either Dirichlet or Neumann conditions. See [Boundary Conditions](@ref boundary-conditions) for the equations behind each.

## Calculation of the diffusion coefficient
There are two ways to set the diffusion coefficient: (a) a constant value per side, given directly via `Di`; or (b) setting both entries of `Di` to `-1.0`, which computes D from the Arrhenius relationship instead — the same applies to single-crystal diffusion. The one exception is `Chemical_Stefan_problem_XT.jl`, where diffusivities are always constant: example D1 uses experimentally determined Fe–Mg diffusivities in olivine and melt [DohmenChakraborty2007a](@cite), [DohmenChakraborty2007b](@cite), [ZhangCherniak2010](@cite), combined via the effective-diffusivity evaluations of [Crank1975](@cite). You can substitute your own calculation method for the diffusion coefficient at any time.
