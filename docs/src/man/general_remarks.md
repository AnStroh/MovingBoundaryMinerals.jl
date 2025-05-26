# General remarks

In this section we would like to make a few generally useful comments on the structure of the codes. 

## Structure
We have tried to structure the individual files as uniformly as possible and to delete unnecessary parameters in the examples so that the user gets a good overview as quickly as possible and understands the individual examples well. Furthermore, the code works with SI units. Subsequently, users are asked to adjust their values or to change the units consequently at their own risk.Most of the codes are created to handle diffusion couples. In these codes we always refer to the left (phase A) and right side (phase B). If two numbers are stored in a variable, the first refers to the left material and the second to the right material. We have also stored generally valid codes in the main_codes folder. The examples can be found in the examples folder. Parameters in the main_codes folder have no physical meaning and act as placeholders such that the code works out of the box. In contrast, the examples show codes with real or non-dimensional values.

## Compositions and concentrations
In our Julia package MOBILE.jl, we specify the concentration in mol. However, since we include any swelling or shrinking processes due to density differences, the composition can also be given in wt% or mol%. The respective parameters such as activation energy or densities must be corrected according to the utilized units. We want to emphasize, that consistence is extremely important here. The units of concentrations (e.g., mol) or compositions (e.g., wt%) are shortened at the end in the diffusion equation (Eq. 1). For more details: Appendices A and B in the related paper to this package (ADD REFERENCE!!!)

## Temperature and distribution coefficient
In general, it is always possible to specify the temporal evolution of the distribution coefficient (KD) and the temperature as vectors. This makes it very easy to handle isothermal as well as non-isothermal problems and deal with constant or changing K_D values. Within the temperature and the KD vector, the first value defines the initial value. The last value defines value at the end of the simulation. If the first and the last entry are the same, the respective parameter is constant. 

## Inner and outer boundaries
The inner boundary at the interface can either be described by the flux balance approach or with total mass balance. Outer boundary conditions at the edges of the modeling domain can be set to Dirichlet or Neumann conditions using BCout. 

## Calculation of the diffusion coefficient
There are two methods for calculating the diffusion coefficient in our package a) a constant diffusion coefficient is used by recording the respective values for the left and right side under the variable Di and b) if both values in Di are replaced with -1.0, D is calculated using the Arrhenius relationship. The same also applies to diffsion processes in single crystals. However, there is one exception. Within Chemical_Stefan_problem.jl, the diffusivities are implemented as constant values.  Within example D1, we specify the diffusion coefficient for Fe-Mg in olivine and in the melt based on experimental values (Dohmen and Chakraborty, 2007a, b; Zhang and Cherniak, 2010) and effective evaluations (Crank, 1975). The user can customize the calculation method of the diffusion coefficients at any time.

## References
Crank, J.: The mathematics of diffusion, 2d ed., Clarendon Press, Oxford, [Eng], 414 pp., 1975.
Dohmen, R. and Chakraborty, S.: Fe–Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine, Phys. Chem. Miner., 34, 597–598, https://doi.org/10.1007/s00269-007-0185-3, 2007a.
Dohmen, R. and Chakraborty, S.: Fe–Mg diffusion in olivine II: point defect chemistry, change of diffusion mechanisms and a model for calculation of diffusion coefficients in natural olivine, Phys. Chem. Miner., 34, 409–430, https://doi.org/10.1007/s00269-007-0158-6, 2007b.
Zhang, Y. and Cherniak, D. J.: Diffusion in Minerals and Melts: Introduction, Rev. Mineral. Geochem., 72, 1–4, https://doi.org/10.2138/rmg.2010.72.1, 2010.

