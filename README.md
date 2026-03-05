[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anstroh.github.io/MovingBoundaryMinerals.jl/dev/)
[![CI](https://github.com/AnStroh/MovingBoundaryMinerals.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/AnStroh/MovingBoundaryMinerals.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15535732.svg)](https://doi.org/10.5281/zenodo.15535732)

# MovingBoundaryMinerals.jl

[MovingBoundaryMinerals.jl](https://github.com/AnStroh/MovingBoundaryMinerals.jl) is a software to model diffusion-controlled growth/resorption processes in mineral couples. 
We tested our software with various analytical and semi-analytical solutions (examples A1-A2, B1-B5, C1). In addition, we present some mineral growth/resorption examples (B6-B7, C2, D1) showing compositional profiles, which can be observed in natural samples. For further details, we recommend reading [Stroh et al. (2025)](https://doi.org/10.5194/egusphere-2025-2511).

> [!NOTE]
> This package is still under active development. 
>  - The benchmarks and examples are working and provide the user with an insight into the capabilities of the package.

## Installation

`MovingBoundaryMinerals.jl` is a registered package and can be added as follows:


```julia
using Pkg; Pkg.add("MovingBoundaryMinerals")
```
or
```julia
julia> ]
pkg> add MovingBoundaryMinerals
```

However, the package is developing and not every feature leads to a new release, one can also do `add MovingBoundaryMinerals#main`, which will clone the main branch of the repository. After installation, you can test the package by running the following commands:

```julia
using MovingBoundaryMinerals
julia> ]
  pkg> test MovingBoundaryMinerals
```
The test can take a while.

## How to cite
Stroh, A., Aellig, P. S., and Moulas, E.: Numerical modelling of diffusion-limited mineral growth for geospeedometry applications, Geosci. Model Dev., 18, 10203–10220, https://doi.org/10.5194/gmd-18-10203-2025, 2025. 

## Funding
The development of this package is supported by the DFG project 524829125 (VECTOR) and by the European Research Council through the MAGMA project, ERC Consolidator Grant \#771143.
