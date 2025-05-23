```@meta
CurrentModule = MOBILE
```# MOBILE.jl

This is the documentation page of the package Diffusion coupled growth.
# Installation

[MOBILE.jl](https://github.com/AnStroh/MOBILE.jl) is not yet registered and can be added as follows:

```julia
using Pkg; Pkg.add("https://github.com/AnStroh/MOBILE.jl.git")
```
or
```julia-repl
julia> ]

(@v1.10) pkg> add https://github.com/AnStroh/MOBILE.jl.git
```

!!! info "Install from a specific branch"
    However, as the API is changing and not every new feature leads to a new release, one can also clone the main branch of the repository:
    ```julia
    add https://github.com/AnStroh/MOBILE.jl.git#main
    ```

After installation, you can test the package by running the following commands:
```julia-repl
using MOBILE

julia> ]

(@v1.10) pkg> test MOBILE
```
