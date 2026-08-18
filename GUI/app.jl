ENV["GKSwstype"] = "100"   # headless GR backend - must be set before `using Plots`

using Oxygen
using HTTP
using Plots
using MovingBoundaryMinerals
using TOML

const GUI_DIR = @__DIR__

include(joinpath(GUI_DIR, "backends", "single_crystal.jl"))
include(joinpath(GUI_DIR, "backends", "diffusion_couple.jl"))
include(joinpath(GUI_DIR, "backends", "thermo_growth.jl"))
include(joinpath(GUI_DIR, "jobs.jl"))
include(joinpath(GUI_DIR, "plotting.jl"))
include(joinpath(GUI_DIR, "results.jl"))
include(joinpath(GUI_DIR, "templates", "layout.jl"))
include(joinpath(GUI_DIR, "templates", "single_crystal.jl"))
include(joinpath(GUI_DIR, "templates", "diffusion_couple.jl"))
include(joinpath(GUI_DIR, "templates", "thermo_growth.jl"))
include(joinpath(GUI_DIR, "templates", "history.jl"))

if Threads.nthreads() == 1
    @warn "Julia was started with only 1 thread. Long simulations (especially thermodynamic growth) will freeze the whole GUI while running. Launch with `julia -t auto` instead (the provided run_gui scripts already do this)."
end

# ------------------------------------------------------------------
# Form field parsing helpers
#
# `label` is the same text shown above the field on the form, so a bad value produces a
# message the user can actually act on (e.g. "'Total time [years]' must be a number (got
# '').") instead of one generic "check your inputs" for every possible mistake.
# ------------------------------------------------------------------
function num(form, key, label)
    raw = get(form, key, "")
    v = tryparse(Float64, raw)
    v === nothing && error("'$(label)' must be a number (got '$(raw)').")
    return v
end
function intg(form, key, label)
    raw = get(form, key, "")
    v = tryparse(Int, raw)
    v === nothing && error("'$(label)' must be a whole number (got '$(raw)').")
    return v
end
run_name_of(form) = get(form, "run_name", "")

"""Extracts a user-facing message from a caught input error. `num`/`intg` above always throw
a plain `ErrorException` with a specific message; anything else (a missing form key, say)
falls back to a generic message rather than leaking an internal exception string."""
error_message(e) = e isa ErrorException ? e.msg : "Invalid input - please check the values entered."

"""Builds the {status:"done", ...} JSON payload for a finished job, with download URLs
pointing at its saved results folder."""
function done_payload(job::Job)
    base = "/results/$(job.result_dir)"
    return Dict(
        "status" => "done",
        "summary" => job.summary,
        "folder" => job.result_dir,
        "png_url" => "$(base)/plot.png",
        "pdf_url" => "$(base)/plot.pdf",
        "data_url" => "$(base)/profile_final.tab",
        "folder_url" => "$(base)/",
    )
end

# ------------------------------------------------------------------
# Pages
# ------------------------------------------------------------------
@get "/" function(req::HTTP.Request)
    return HTTP.Response(302, ["Location" => "/single-crystal"])
end

@get "/single-crystal" function(req::HTTP.Request)
    return html(single_crystal_page())
end

@get "/diffusion-couple" function(req::HTTP.Request)
    return html(diffusion_couple_page())
end

@get "/thermo-growth" function(req::HTTP.Request)
    return html(thermo_growth_page())
end

@get "/history" function(req::HTTP.Request)
    runs = []
    if isdir(RESULTS_DIR)
        for folder in sort(readdir(RESULTS_DIR); rev = true)
            inputs_path = joinpath(RESULTS_DIR, folder, "inputs.toml")
            isfile(inputs_path) || continue
            try
                record = TOML.parsefile(inputs_path)
                push!(runs, (folder = folder, mode = record["mode"], timestamp = record["timestamp"]))
            catch
                continue   # skip anything unreadable/partially written rather than fail the whole page
            end
        end
    end
    return html(history_page(runs))
end

# ------------------------------------------------------------------
# Run endpoints - each starts a background job and returns its id
# ------------------------------------------------------------------
@post "/single-crystal/run" function(req::HTTP.Request)
    form = formdata(req)
    local kwargs
    try
        kwargs = (
            D0 = num(form, "D0", "Pre-exponential factor D0 [m²/s]"),
            Ea1 = num(form, "Ea1", "Activation energy [J/mol]"),
            L = num(form, "L", "Domain length [m]"),
            Cstart = num(form, "Cstart", "Initial composition [-]"),
            Cinf = num(form, "Cinf", "Composition at infinity [-]"),
            rho = num(form, "rho", "Density [kg/m³]"),
            n = intg(form, "n", "Geometry"),
            Tstart_C = num(form, "Tstart_C", "Starting temperature [°C]"),
            Tstop_C = num(form, "Tstop_C", "End temperature [°C]"),
            t_tot_years = num(form, "t_tot_years", "Total time [years]"),
            nx = intg(form, "nx", "Number of nodes"),
            CFL = num(form, "CFL", "CFL number"),
        )
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:single_crystal) do should_stop
        res = run_single_crystal(; kwargs..., should_stop)
        p = plot_single_crystal(res)
        profile = profile_data_single_crystal(res)
        extra = Dict("final_time_s" => res.t_tot, "final_time_years" => res.t_tot / (365.25 * 24 * 60 * 60))
        folder = save_run_outputs("single_crystal", p, profile, kwargs, extra; run_name)
        summary = "Done. Simulated $(round(res.t_tot / (365.25*24*60*60), digits = 3)) years."
        (folder, summary)
    end
    id === nothing && return json(Dict("error" => "A simulation is already running."), status = 409)
    return json(Dict("id" => id))
end

@post "/diffusion-couple/run" function(req::HTTP.Request)
    form = formdata(req)
    local kwargs
    try
        kwargs = (
            D0_left = num(form, "D0_left", "Pre-exponential factor, left [m²/s]"),
            D0_right = num(form, "D0_right", "Pre-exponential factor, right [m²/s]"),
            Ea_left = num(form, "Ea_left", "Activation energy, left [J/mol]"),
            Ea_right = num(form, "Ea_right", "Activation energy, right [J/mol]"),
            Ri_interface = num(form, "Ri_interface", "Initial interface position [m]"),
            Ri_total = num(form, "Ri_total", "Total domain length [m]"),
            Cl_i = num(form, "Cl_i", "Initial composition, left [-]"),
            Cr_i = num(form, "Cr_i", "Initial composition, right [-]"),
            KD_start = num(form, "KD_start", "Distribution coefficient K_D at t=0"),
            KD_end = num(form, "KD_end", "Distribution coefficient K_D at end"),
            n = intg(form, "n", "Geometry"),
            Tstart_C = num(form, "Tstart_C", "Starting temperature [°C]"),
            Tstop_C = num(form, "Tstop_C", "End temperature [°C]"),
            t_tot_years = num(form, "t_tot_years", "Total time [years]"),
        )
        if kwargs.Ri_interface >= kwargs.Ri_total
            error("'Initial interface position [m]' must be smaller than 'Total domain length [m]'.")
        end
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:diffusion_couple) do should_stop
        res = run_diffusion_couple(; kwargs..., should_stop)
        p = plot_diffusion_couple(res)
        profile = profile_data_two_phase(res)
        extra = Dict("final_time_s" => res.t_tot, "mass_balance_error_pct" => res.MB_Error[end])
        folder = save_run_outputs("diffusion_couple", p, profile, kwargs, extra; run_name)
        summary = "Done. Final mass-balance error: $(round(res.MB_Error[end], digits = 4))%."
        (folder, summary)
    end
    id === nothing && return json(Dict("error" => "A simulation is already running."), status = 409)
    return json(Dict("id" => id))
end

@post "/thermo-growth/run" function(req::HTTP.Request)
    form = formdata(req)
    local kwargs
    try
        kwargs = (
            Ri = num(form, "Ri", "Initial interface radius [m]"),
            n = intg(form, "n", "Geometry"),
            Tstart_C = num(form, "Tstart_C", "Starting temperature [°C]"),
            Tstop_C = num(form, "Tstop_C", "End temperature [°C]"),
            t_tot_days = num(form, "t_tot_days", "Total time [days]"),
            P = num(form, "P", "Pressure [Pa]"),
            D0 = num(form, "D0", "Pre-exponential factor D0 [m²/s]"),
            Ea = num(form, "Ea", "Activation energy [J/mol]"),
            CompInt = num(form, "CompInt", "Bulk composition of interest [-]"),
        )
        if !(0.0 < kwargs.CompInt < 1.0)
            error("'Bulk composition of interest [-]' must be between 0 and 1.")
        end
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:thermo_growth) do should_stop
        res = run_thermo_growth(; kwargs..., should_stop)
        p = plot_thermo_growth(res)
        profile = profile_data_two_phase(res)
        extra = Dict("final_time_s" => res.t_tot, "mass_balance_error_pct" => res.MB_Error[end])
        folder = save_run_outputs("thermo_growth", p, profile, kwargs, extra; run_name)
        summary = "Done. Final mass-balance error: $(round(res.MB_Error[end], digits = 4))%."
        (folder, summary)
    end
    id === nothing && return json(Dict("error" => "A simulation is already running."), status = 409)
    return json(Dict("id" => id))
end

# ------------------------------------------------------------------
# Job polling and cancellation
# ------------------------------------------------------------------
@get "/jobs/{id}/status" function(req::HTTP.Request, id::String)
    job = get_job(id)
    job === nothing && return json(Dict("status" => "failed", "error" => "Unknown job id."), status = 404)
    if job.status == :running
        return json(Dict("status" => "running"))
    elseif job.status == :done
        return json(done_payload(job))
    elseif job.status == :cancelled
        return json(Dict("status" => "cancelled"))
    else
        return json(Dict("status" => "failed", "error" => job.error))
    end
end

@post "/jobs/{id}/cancel" function(req::HTTP.Request, id::String)
    ok = cancel_job!(id)
    return json(Dict("ok" => ok))
end

mkpath(RESULTS_DIR)   # the mount below requires the folder to exist even before any run has happened

staticfiles(joinpath(GUI_DIR, "static"), "static")
# dynamicfiles (not staticfiles!): results/ gains new files at runtime as simulations
# finish, and staticfiles only scans its folder once, at server startup.
dynamicfiles(RESULTS_DIR, "results")

# ------------------------------------------------------------------
# Launch a chromeless "app window" browser view pointing at the server
# ------------------------------------------------------------------
function open_app_window(url::String)
    try
        @static if Sys.islinux()
            for browser in ("google-chrome", "google-chrome-stable", "chromium-browser", "chromium")
                if !isnothing(Sys.which(browser))
                    run(`$browser --app=$url --window-size=1200,800`; wait = false)
                    return
                end
            end
            run(`xdg-open $url`; wait = false)
        elseif Sys.isapple()
            for browser in ("Google Chrome", "Microsoft Edge")
                try
                    run(`open -na $browser --args --app=$url`; wait = false)
                    return
                catch
                end
            end
            run(`open $url`; wait = false)
        elseif Sys.iswindows()
            for browser in ("chrome", "msedge")
                try
                    run(`cmd /c start $browser --app=$url`; wait = false)
                    return
                catch
                end
            end
            run(`cmd /c start $url`; wait = false)
        end
    catch e
        @warn "Could not open a browser window automatically." exception=e
    end
    println("If a browser window did not open automatically, go to: $url")
end

function main()
    host = "127.0.0.1"
    port = 8811
    url = "http://$host:$port"
    println("Starting MovingBoundaryMinerals.jl GUI at $url ...")
    server = serve(host = host, port = port, async = true, show_banner = false)
    for _ in 1:50
        try
            HTTP.get(url; readtimeout = 1)
            break
        catch
            sleep(0.2)
        end
    end
    open_app_window(url)
    wait(server)
end

main()
