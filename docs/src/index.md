```@meta
CurrentModule = Diff_Coupled
```# Diff_Coupled.jl

This is the documentation page of the package Diffusion coupled growth.
# Installation

[Diffusion-coupled-growth](https://github.com/AnStroh/Diffusion-coupled-growth) is not yet registered and can be added as follows:

```julia
using Pkg; Pkg.add("https://github.com/AnStroh/Diffusion-coupled-growth.git")
```
or
```julia-repl
julia> ]

(@v1.10) pkg> add https://github.com/AnStroh/Diffusion-coupled-growth.git
```

!!! info "Install from a specific branch"
    However, as the API is changing and not every new feature leads to a new release, one can also clone the main branch of the repository:
    ```julia
    add Jushttps://github.com/AnStroh/Diffusion-coupled-growth.git#main
    ```

After installation, you can test the package by running the following commands:
```julia-repl
using Diff_Coupled

julia> ]

(@v1.10) pkg> test Diff_Coupled
```
