using Plots

const FONT_SIZE = 12.0
const PLOT_DPI  = 300

"""Package + Zenodo software-archive reference stamped on every plot and saved file, so a
shared/downloaded result stays traceable back to its source even without the surrounding
folder. The version is read from the loaded package itself (not hardcoded), so this never
goes stale."""
const CITATION_LINE = "MovingBoundaryMinerals.jl v$(pkgversion(MovingBoundaryMinerals)) - Stroh et al. (2025), doi:10.5281/zenodo.15535732"

"""Appends `CITATION_LINE` as a thin, borderless footer strip below the axes, so it sits
truly outside the plotted data (not overlapping any curve) and outside the axes box itself.
Plots.jl's `title` attribute is the only plot-level "outside the axes" text mechanism, but it
can only render above the plot - never below - so a bottom-of-figure caption has to be built
as its own tiny blank subplot stacked underneath the real one via a layout, rather than as a
single-plot attribute."""
function stamp_citation(p)
    footer = plot(xlim = (0, 1), ylim = (0, 1), framestyle = :none, grid = false, legend = false,
                   background_color_inside = :match)
    annotate!(footer, 1.0, 0.5, Plots.text(CITATION_LINE, 6, :grey50, :right, :vcenter))
    return plot(p, footer, layout = grid(2, 1, heights = [0.94, 0.06]), dpi = PLOT_DPI)
end

"""Render a Plots.jl figure to PNG bytes, in-memory (no temp file, no display() call)."""
function render_png(p)
    io = IOBuffer()
    show(io, MIME"image/png"(), p)
    return take!(io)
end

function plot_single_crystal(res)
    (; x, C, x0, C0) = res
    p = plot(x .* 1000, C, lw = 2, label = "Current composition")
    p = plot!(x0 .* 1000, C0, label = "Initial composition", color = :black, linestyle = :dash,
              xlabel = "x [mm]", ylabel = "C [-]", lw = 1.5, grid = :on, dpi = PLOT_DPI,
              legendfontsize = FONT_SIZE - 2, guidefontsize = FONT_SIZE, tickfontsize = FONT_SIZE - 1,
              legend_foreground_color = :transparent)
    return stamp_citation(p)
end

function plot_diffusion_couple(res)
    (; x_left, x_right, x0, C_left, C_right, C0, Ri) = res
    maxC = maximum([maximum(C_left), maximum(C_right)])
    p = plot(x_left .* 1000, C_left, lw = 2, label = "Left side")
    p = plot!(x_right .* 1000, C_right, lw = 2, label = "Right side")
    p = plot!(x0 .* 1000, C0, label = "Initial composition", color = :black, linestyle = :dash,
              xlabel = "x [mm]", ylabel = "C [-]", lw = 1.5, grid = :on, dpi = PLOT_DPI,
              legendfontsize = FONT_SIZE - 2, guidefontsize = FONT_SIZE, tickfontsize = FONT_SIZE - 1,
              legend_foreground_color = :transparent)
    p = plot!([Ri[1]; Ri[1]] .* 1000, [0; 1] .* maxC, color = :grey68, linestyle = :dashdot, lw = 2, label = "Interface")
    return stamp_citation(p)
end

function plot_thermo_growth(res)
    (; x_left, x_right, x0, C_left, C_right, C0) = res
    maxC = maximum([maximum(C_left), maximum(C_right)])
    p = plot(x_left .* 1e3, C_left, lw = 2, label = "Left side")
    p = plot!(x_right .* 1e3, C_right, lw = 2, label = "Right side")
    p = plot!(x0 .* 1e3, C0, label = "Initial composition", color = :black, linestyle = :dash,
              xlabel = "x [mm]", ylabel = "X_Fe [-]", lw = 1.5, grid = :on, dpi = PLOT_DPI,
              legendfontsize = FONT_SIZE - 2, guidefontsize = FONT_SIZE, tickfontsize = FONT_SIZE - 1,
              legend_foreground_color = :transparent)
    p = plot!([x_left[end]; x_left[end]] .* 1e3, [0; 1 * (maxC + 0.01)], color = :grey68, linestyle = :dashdot,
              lw = 2, label = "Interface", ylim = [C0[1] - 0.05; 1 * (maxC + 0.01)])
    return stamp_citation(p)
end

# ------------------------------------------------------------------
# Raw (distance, composition) profile data, for exporting alongside the plot so users can
# work with the numbers directly. The initial and final grids are kept as separate
# (x, C) pairs rather than one merged table: modes with a moving interface (thermo growth)
# regrid during the run, so the final grid can have a different number of points than the
# initial one - a single shared-length table would silently misalign or truncate data.
# ------------------------------------------------------------------

"""(x_initial, C_initial, x_final, C_final) for the single-crystal mode (fixed grid)."""
profile_data_single_crystal(res) = (x_initial = vec(res.x0), C_initial = vec(res.C0),
                                     x_final = vec(res.x), C_final = vec(res.C))

"""(x_initial, C_initial, x_final, C_final) for the two-phase modes (diffusion couple /
thermo growth). The final grid is read from `x_left`/`x_right` directly (not assumed to
match the initial grid `x0` in length), since thermo-growth regrids as the interface moves."""
function profile_data_two_phase(res)
    return (x_initial = vec(res.x0), C_initial = vec(res.C0),
            x_final = vcat(vec(res.x_left), vec(res.x_right)),
            C_final = vcat(vec(res.C_left), vec(res.C_right)))
end
