```@meta
CurrentModule = MovingBoundaryMinerals
```
# [Digitizing the phase diagram and extracting coefficients](@id digitization)

## Phase diagrams

The phase diagram, which shows the composition of the system as a function of temperature, includes all important information about chemical reactions between mineral assemblages. Therefore, it is necessary to describe the equilibrium conditions at the interface. Equations for the phase transformation lines contain all thermodynamically important data. These equations can be described by a quadratic equation (Eq. 1). Within our code, two transformation lines are always required to determine the composition of the components at the edges of a stable assemblage field. The user needs to define coefficients (a, b, c) for both transformation lines as input prior to the start of the model.

## Work flow

We have used the `PlotDigitizer.jl` package in Julia [StrohFrasunkiewicz2025](@cite) ([https://github.com/AnStroh/PlotDigitizer.jl](https://github.com/AnStroh/PlotDigitizer.jl), version v0.1.0) to digitize the phase diagram. A shortened version of the package can be found in the `additional_codes` folder under `digitizePlot.jl`. The package can be used to digitize reaction lines and export the X-T coordinates. The coordinates can then be used to calculate the coefficients for the reaction lines. We provide a short guideline here.

!!! warning "Extra dependency required"
    `digitizePlot.jl` uses `GLMakie` for the interactive image viewer. `GLMakie` is **not** a dependency of `MovingBoundaryMinerals.jl` itself (it is heavy and only needed for this optional, interactive digitizing step), so you need to install it yourself first: `using Pkg; Pkg.add("GLMakie")`. `FileIO` and `DelimitedFiles`, the other two packages `digitizePlot.jl` uses, are already part of `MovingBoundaryMinerals.jl`'s dependencies.

1. Define the path to the image of the phase diagram (`.png` or `.jpg` file) as a variable. It is important that the image only shows the phase diagram (like `examples/Examples_phase_diagram/Ol_Phase_diagram_without_framework.png`) and has no border. The phase diagram should show the composition on the x-axis and the temperature on the y-axis.
2. Define the limits of the axes in `X_BC` and `Y_BC`.
3. Run `digitizePlot(X_BC, Y_BC, file_name)` to set or delete points on the two reaction lines, switch between lines, and export the digitized X-T coordinates to CSV.
4. Use `CalculateReactionLine.jl` (also in `additional_codes`) to fit the digitized points via least squares to

```math
X(T) = a T^2 + b T + c \qquad (1)
```

where `X` is the composition at the phase transition, `T` is the temperature in K, and `a`, `b`, `c` are the coefficients of the quadratic equation. The fitted coefficients are written to a CSV file (see `examples/Examples_phase_diagram/Coefficients_Reaction_lines.csv` for the file used by example [D1](https://github.com/AnStroh/MovingBoundaryMinerals.jl/blob/main/examples/D1.jl)) and can be loaded with [`coeff_trans_line`](@ref) and evaluated at any temperature with [`composition`](@ref).

!!! tip "Fewer points can be more accurate"
    CSV files with fewer, carefully placed points tend to give a more accurate quadratic fit near the ends of the composition range (X close to 0 or 1) than files with many densely digitized points. Always check the fit against the original phase diagram near its boundaries.

This step can be skipped entirely if the coefficients `a`, `b`, `c` are already known: in that case they can be entered directly, e.g. via `eq_values` in the example scripts.
