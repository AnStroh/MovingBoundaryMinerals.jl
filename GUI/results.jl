using Dates, DelimitedFiles, TOML, Plots

const RESULTS_DIR = joinpath(@__DIR__, "results")

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

"""
    save_run_outputs(mode, p, profile, kwargs, extra_results; run_name="") -> folder_name

Saves everything from one GUI run into a new timestamped folder under `GUI/results/`, so
results are never silently overwritten and every run is independently reproducible:

- `plot.png`, `plot.pdf` - the result figure at 300 dpi
- `profile_initial.tab`, `profile_final.tab` - the raw (distance, composition) data, tab-delimited
- `inputs.toml` - the mode, every input parameter used, and key scalar results (final time,
  final mass-balance error), so the run can be reproduced exactly

`profile` is a `(; x_initial, C_initial, x_final, C_final)` named tuple (see `plotting.jl`).
`extra_results` is a Dict of additional scalar results to record (e.g. final time, MB error).
If `run_name` is non-empty, it's appended to the timestamp to make the folder easier to find.
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

    record = Dict(
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
