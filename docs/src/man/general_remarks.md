```@meta
CurrentModule = MovingBoundaryMinerals
```
# General Remarks

In this section we would like to make a few generally useful comments on the structure of the codes.

## Structure
We have tried to structure the individual files as uniformly as possible and to delete unnecessary parameters in the examples so that the user gets a good overview as quickly as possible and understands the individual examples well. Furthermore, the code works with SI units. Therefore, users are asked to adjust their values or to change the units consequently at their own risk. Most of the codes are created to handle diffusion couples. In these codes we always refer to the left (phase A) and right side (phase B). If two numbers are stored in a variable, the first refers to the left material and the second to the right material. We have also stored generally valid codes in the main_codes folder. The examples can be found in the examples folder. Parameters in the main_codes folder have no physical meaning and act as placeholders such that the code works out of the box. In contrast, the examples show codes with real or non-dimensional values.

## The `main_codes` templates

The [`main_codes`](https://github.com/AnStroh/MovingBoundaryMinerals.jl/tree/main/main_codes) folder holds four bare-bones templates, one per problem family, meant as a starting point for your own model rather than a worked example. Each corresponds to one of the [example](@ref "List of examples") families:

| Template | Problem family | Corresponding examples |
|:--|:--|:--|
| `Simple_Diff_1D_general.jl` | Single-crystal diffusion, no interface | `Simple_Diff`, A1, A2 |
| `Diff_couple_1D_general_Flux.jl` | Diffusion couple, flux-balance interface condition | `Diff_couple_Flux[_growth]`, B1–B7 |
| `Diff_couple_1D_general_MB.jl` | Diffusion couple, total mass-balance interface condition | `Diff_couple_MB[_growth]`, C1–C2 |
| `Chemical_Stefan_problem_XT.jl` | Thermodynamically constrained (Stefan) growth/resorption | D1, D2 |

## Compositions and concentrations
In our Julia package MovingBoundaryMinerals.jl, we specify the concentration in mol. However, since we include any swelling or shrinking processes due to density differences, the composition can also be given in wt% or mol%. The respective parameters such as activation energy or densities must be corrected according to the utilized units. We want to emphasize that consistency is extremely important here. The units of concentrations (e.g., mol) or compositions (e.g., wt%) are shortened at the end in the diffusion equation (Eq. 1). For more details, see Appendices A and B of [Stroh2025](@cite).

## Temperature and distribution coefficient
In general, it is always possible to specify the temporal evolution of the distribution coefficient (KD) and the temperature as vectors. This makes it very easy to handle isothermal as well as non-isothermal problems and deal with constant or changing K_D values. Within the temperature and the KD vector, the first value defines the initial value. The last value defines the value at the end of the simulation. If the first and the last entry are the same, the respective parameter is constant.

## Inner and outer boundaries
The inner boundary at the interface can either be described by the flux balance approach or with total mass balance. Outer boundary conditions at the edges of the modeling domain can be set to Dirichlet or Neumann conditions using BCout.

## Calculation of the diffusion coefficient
There are two methods for calculating the diffusion coefficient in our package: (a) a constant diffusion coefficient is used, given by the respective values for the left and right side under the variable `Di`; and (b) if both values in `Di` are replaced with -1.0, D is calculated using the Arrhenius relationship. The same also applies to diffusion processes in single crystals. However, there is one exception. Within Chemical_Stefan_problem.jl, the diffusivities are implemented as constant values. Within example D1, we specify the diffusion coefficient for Fe-Mg in olivine and in the melt based on experimental values [DohmenChakraborty2007a](@cite), [DohmenChakraborty2007b](@cite), [ZhangCherniak2010](@cite) and effective evaluations [Crank1975](@cite). The user can customize the calculation method of the diffusion coefficients at any time.
