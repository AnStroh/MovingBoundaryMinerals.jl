```@meta
CurrentModule = MovingBoundaryMinerals
```
# MovingBoundaryMinerals.jl

[MovingBoundaryMinerals.jl](https://github.com/AnStroh/MovingBoundaryMinerals.jl) is a software package to model diffusion-controlled growth/resorption processes in mineral couples. It solves the moving-boundary (Stefan) problem directly, tracking the interface between two phases as part of the solution rather than prescribing its position or growth rate a priori. See [General Remarks](@ref) for an overview and [The Stefan Problem](@ref stefan-problem) for the underlying physics and numerics.

For further details, we recommend reading [Stroh2025](@cite).

# Installation

`MovingBoundaryMinerals.jl` requires Julia 1.10 or later (tested on the latest LTS release, 1.11, and pre-release Julia on Linux, macOS, and Windows — see the [CI workflow](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/.github/workflows/CI.yml) for the exact matrix). It is a registered package and can be added as follows:

```julia
using Pkg; Pkg.add("MovingBoundaryMinerals")
```
or
```julia-repl
julia> ]

(@v1.10) pkg> add MovingBoundaryMinerals
```

!!! info "Install from a specific branch"
    As the API is changing and not every new feature leads to a new release, one can also clone the main branch of the repository directly:
    ```julia
    add MovingBoundaryMinerals#main
    ```

After installation, you can test the package by running the following commands:
```julia-repl
using MovingBoundaryMinerals

julia> ]

(@v1.10) pkg> test MovingBoundaryMinerals
```
The test suite can take a while to run.

# Getting started

New to the package? [Getting Started](@ref getting-started) walks through running your first example end to end. Already know the package and just need a reminder of what a parameter does? See the [Quick Reference](@ref quick-reference). Don't want to write or edit Julia code at all? See the [GUI](@ref gui).

# AI use

ChatGPT and Copilot were used in translating the original MATLAB© codes into Julia. Furthermore, Copilot was used in expanding the documentation of the functions. Moreover, we asked Claude to help with the documentation, finding/fixing bugs, increasing readability, and creating a GUI. All the results were checked and are approved by the authors.
