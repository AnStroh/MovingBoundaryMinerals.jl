using Dates, DelimitedFiles, TOML, Plots, XLSX, JLD2

"""Where completed runs are saved. Overridable via `MBM_GUI_RESULTS_DIR` so `GUI/test/` can
point it at a throwaway directory instead of polluting a real `GUI/results/` folder; unset,
it defaults to `GUI/results/` exactly as before."""
const RESULTS_DIR = get(ENV, "MBM_GUI_RESULTS_DIR", joinpath(@__DIR__, "results"))

"""Turn a free-text run name into something safe to use as part of a folder name."""
function sanitize_name(name::AbstractString)
    cleaned = replace(strip(name), r"[^A-Za-z0-9_\-]+" => "_")
    return strip(cleaned, '_')
end

"""Write a (x, C) profile as a tab-delimited file with a header row. Starts with a `#`-prefixed
comment line naming the source package/paper, so the file stays traceable if shared on its own -
most tools (Excel, pandas, numpy) either skip `#` lines automatically or make it trivial to."""
function write_profile(path::String, x::Vector{Float64}, C::Vector{Float64})
    open(path, "w") do io
        println(io, "# ", CITATION_LINE)
        writedlm(io, ["distance_m" "composition"], '\t')
        writedlm(io, hcat(x, C), '\t')
    end
end

"""Write the thermodynamic-growth mode's temperature-time path as a tab-delimited file, in the
same "time [days], temperature [°C]" units the form's textarea takes - not the seconds/Kelvin
`t_user`/`T_user` are rescaled to internally - so it can be read back into the form's path field
directly (copy-paste) to reproduce or tweak this run's T-t history."""
function write_path(path::String, t_user::AbstractVector, T_user::AbstractVector)
    open(path, "w") do io
        println(io, "# ", CITATION_LINE)
        writedlm(io, ["time_days" "temperature_C"], '\t')
        writedlm(io, hcat(t_user ./ (24 * 60 * 60), T_user .- 273.0), '\t')
    end
end

"""Write the initial/final profiles plus every input parameter and result into Julia's own
native `.jld2` format (via JLD2.jl) - `source`, `mode`, `run_name`, `x_initial`, `C_initial`,
`x_final`, `C_final`, `parameters`, `results` (the last two as nested `Dict`s), read back with
`load("profiles.jld2")` - for users staying in Julia who want the data back as-is, without also
opening `inputs.toml`."""
function write_jld2(path::String, profile, kwargs, extra_results::Dict, mode::String, run_name::AbstractString)
    jldsave(path;
        source = CITATION_LINE,
        mode = mode,
        run_name = run_name,
        x_initial = profile.x_initial,
        C_initial = profile.C_initial,
        x_final = profile.x_final,
        C_final = profile.C_final,
        parameters = Dict(string(k) => v for (k, v) in pairs(kwargs)),
        results = extra_results,
    )
end

"""Write the initial/final profiles plus every input parameter and result into a single
Excel `.xlsx` workbook - a "Profile" sheet (initial and final columns side by side; they may
have different lengths since some modes regrid) and a "Parameters" sheet (mode, run name,
every input, then every result), so Excel users have everything needed to inspect or
reproduce the run without also opening `inputs.toml`."""
function write_xlsx(path::String, profile, kwargs, extra_results::Dict, mode::String, run_name::AbstractString)
    XLSX.openxlsx(path; mode = "w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Profile")
        sheet["A1"] = CITATION_LINE
        sheet["A3"] = "distance_initial_m"
        sheet["B3"] = "composition_initial"
        sheet["C3"] = "distance_final_m"
        sheet["D3"] = "composition_final"
        sheet["A4", dim = 1] = profile.x_initial
        sheet["B4", dim = 1] = profile.C_initial
        sheet["C4", dim = 1] = profile.x_final
        sheet["D4", dim = 1] = profile.C_final

        params = XLSX.addsheet!(xf, "Parameters")
        params["A1"] = CITATION_LINE
        params["A3"] = "mode"
        params["B3"] = mode
        params["A4"] = "run_name"
        params["B4"] = run_name
        row = 5
        for (k, v) in pairs(kwargs)
            params["A$(row)"] = string(k)
            # A Vector assigned directly to a single cell silently spills across the following
            # rows/cells instead of erroring, which would corrupt this one-row-per-parameter
            # layout (e.g. thermo-growth's t_user/T_user path nodes) - write it as text instead.
            params["B$(row)"] = v isa AbstractVector ? join(v, ", ") : v
            row += 1
        end
        for (k, v) in extra_results
            params["A$(row)"] = string(k)
            params["B$(row)"] = v
            row += 1
        end
    end
end

"""
    save_run_outputs(mode, p, profile, kwargs, extra_results; run_name="") -> folder_name

Saves everything from one GUI run into a new timestamped folder under `GUI/results/`, so
results are never silently overwritten and every run is independently reproducible:

- `plot.png`, `plot.pdf` - the result figure at 300 dpi
- `profile_initial.tab`, `profile_final.tab` - the raw (distance, composition) data, tab-delimited
- `inputs.toml` - the mode, every input parameter used, and key scalar results (final time,
  final mass-balance error), so the run can be reproduced exactly
- `profiles.xlsx`, `profiles.jld2` - the same initial/final profile data and inputs.toml
  content, bundled into one file each for users who work in Excel or Julia (`profiles.jld2`,
  via JLD2.jl) rather than plain text

`profile` is a `(; x_initial, C_initial, x_final, C_final)` named tuple (see `plotting.jl`).
`extra_results` is a Dict of additional scalar results to record (e.g. final time, MB error).
If `run_name` is non-empty, it's appended to the timestamp to make the folder easier to find.

For the thermodynamic growth mode, `kwargs` additionally carries `t_user`/`T_user` (the
temperature-time path, in seconds/Kelvin after the "Total time" field's rescaling) - when both
are present, an extra `path.tab` is written alongside the profiles, giving that path back in the
same days/°C units and format the form's path field takes, for reuse or reproduction.
"""
function save_run_outputs(mode::String, p, profile, kwargs, extra_results::Dict; run_name::AbstractString = "")
    timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    safe_name = sanitize_name(run_name)
    folder_name = isempty(safe_name) ? timestamp : "$(timestamp)_$(safe_name)"
    folder = joinpath(RESULTS_DIR, folder_name)
    mkpath(folder)

    savefig(p, joinpath(folder, "plot.png"))
    savefig(p, joinpath(folder, "plot.pdf"))

    write_profile(joinpath(folder, "profile_initial.tab"), profile.x_initial, profile.C_initial)
    write_profile(joinpath(folder, "profile_final.tab"), profile.x_final, profile.C_final)
    if haskey(kwargs, :t_user) && haskey(kwargs, :T_user)
        write_path(joinpath(folder, "path.tab"), kwargs.t_user, kwargs.T_user)
    end
    write_xlsx(joinpath(folder, "profiles.xlsx"), profile, kwargs, extra_results, mode, run_name)
    write_jld2(joinpath(folder, "profiles.jld2"), profile, kwargs, extra_results, mode, run_name)

    record = Dict(
        "source" => CITATION_LINE,
        "mode" => mode,
        "timestamp" => timestamp,
        "run_name" => run_name,
        "parameters" => Dict(string(k) => v for (k, v) in pairs(kwargs)),
        "results" => extra_results,
    )
    open(joinpath(folder, "inputs.toml"), "w") do io
        TOML.print(io, record)
    end

    return folder_name
end
