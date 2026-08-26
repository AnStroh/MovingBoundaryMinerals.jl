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
"""`min`/`strict` mirror (and, for `strict`, tighten) the client-side `min=` already declared
on the `<input>` in the templates - the browser's own HTML5 validation stops an ordinary user
from ever submitting an out-of-range value, but this endpoint is also a plain HTTP API that
nothing stops a raw request (or a browser with JS disabled) from hitting directly, so the
`num`/`intg` -> `num()`/`solve_soe()` path they eventually feed into needs to be able to trust
the numbers itself. `strict = true` for quantities where the boundary value (`0`, typically) is
itself degenerate - a `dx == 0` domain, a `CFL == 0` that never advances `dt` and hangs the
solver forever - rather than merely uninteresting."""
function num(form, key, label; min = nothing, strict = false)
    raw = get(form, key, "")
    v = tryparse(Float64, raw)
    v === nothing && error("'$(label)' must be a number (got '$(raw)').")
    if min !== nothing
        ok = strict ? (v > min) : (v >= min)
        ok || error("'$(label)' must be $(strict ? "greater than" : "at least") $(min) (got $(v)).")
    end
    return v
end
function intg(form, key, label; min = nothing)
    raw = get(form, key, "")
    v = tryparse(Int, raw)
    v === nothing && error("'$(label)' must be a whole number (got '$(raw)').")
    min !== nothing && v < min && error("'$(label)' must be at least $(min) (got $(v)).")
    return v
end
run_name_of(form) = get(form, "run_name", "")

"""Parses a temperature-time path textarea (one 'time_days, temperature_C' pair per line) into
`(t_user, T_user)` in the [s]/[K] units the physics functions expect - a 2-line path reproduces
D1.jl's linear cooling, more lines reproduce D2.jl's non-monotonic paths. Blank lines are
skipped so the field can be padded for readability. Validates the same invariants D2.jl itself
checks (matching lengths is automatic from the pairwise parsing here): at least 2 points, time
strictly increasing, and the first time value is 0."""
function path(form, key, label)
    raw = get(form, key, "")
    t_days = Float64[]
    T_C = Float64[]
    for (i, line) in enumerate(split(raw, '\n'))
        stripped = strip(line)
        isempty(stripped) && continue
        parts = split(stripped, ',')
        length(parts) != 2 &&
            error("'$(label)': line $(i) ('$(stripped)') must be 'time, temperature' with a single comma.")
        t = tryparse(Float64, strip(parts[1]))
        Tc = tryparse(Float64, strip(parts[2]))
        (t === nothing || Tc === nothing) &&
            error("'$(label)': line $(i) ('$(stripped)') must be two numbers separated by a comma.")
        push!(t_days, t)
        push!(T_C, Tc)
    end
    length(t_days) < 2 && error("'$(label)' needs at least 2 points (a start and an end).")
    issorted(t_days) || error("'$(label)': time values must strictly increase from one line to the next.")
    t_days[1] != 0.0 && error("'$(label)': the first time value must be 0.")
    return t_days .* (24 * 60 * 60), T_C .+ 273.0
end

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
        "xlsx_url" => "$(base)/profiles.xlsx",
        "jld2_url" => "$(base)/profiles.jld2",
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
                push!(runs, (
                    folder = folder,
                    mode = record["mode"],
                    timestamp = record["timestamp"],
                    run_name = get(record, "run_name", ""),
                ))
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
            D0 = num(form, "D0", "Pre-exponential factor D0 [m²/s]"; min = 0, strict = true),
            Ea1 = num(form, "Ea1", "Activation energy [J/mol]"; min = 0),
            L = num(form, "L", "Domain length [m]"; min = 0, strict = true),
            Cstart = num(form, "Cstart", "Initial composition [-]"; min = 0),
            Cinf = num(form, "Cinf", "Composition at infinity [-]"; min = 0),
            n = intg(form, "n", "Geometry"),
            Tstart_C = num(form, "Tstart_C", "Starting temperature [°C]"),
            Tstop_C = num(form, "Tstop_C", "End temperature [°C]"),
            t_tot_years = num(form, "t_tot_years", "Total time [years]"; min = 0, strict = true),
            nx = intg(form, "nx", "Number of nodes"; min = 2),
            CFL = num(form, "CFL", "CFL number"; min = 0, strict = true),
        )
        kwargs.n in (1, 2, 3) || error("'Geometry' must be 1, 2, or 3.")
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:single_crystal) do should_stop, report_progress
        res = run_single_crystal(; kwargs..., should_stop, report_progress)
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
            D0_left = num(form, "D0_left", "Pre-exponential factor, left [m²/s]"; min = 0, strict = true),
            D0_right = num(form, "D0_right", "Pre-exponential factor, right [m²/s]"; min = 0, strict = true),
            Ea_left = num(form, "Ea_left", "Activation energy, left [J/mol]"; min = 0),
            Ea_right = num(form, "Ea_right", "Activation energy, right [J/mol]"; min = 0),
            Ri_interface = num(form, "Ri_interface", "Initial interface position [m]"; min = 0, strict = true),
            Ri_total = num(form, "Ri_total", "Total domain length [m]"; min = 0, strict = true),
            Cl_i = num(form, "Cl_i", "Initial composition, left [-]"; min = 0),
            Cr_i = num(form, "Cr_i", "Initial composition, right [-]"; min = 0),
            KD_start = num(form, "KD_start", "Distribution coefficient K_D at t=0"; min = 0, strict = true),
            KD_end = num(form, "KD_end", "Distribution coefficient K_D at end"; min = 0, strict = true),
            V_ip = num(form, "V_ip", "Interface velocity [m/s]"),
            n = intg(form, "n", "Geometry"),
            Tstart_C = num(form, "Tstart_C", "Starting temperature [°C]"),
            Tstop_C = num(form, "Tstop_C", "End temperature [°C]"),
            t_tot_years = num(form, "t_tot_years", "Total time [years]"; min = 0, strict = true),
            nx_left = intg(form, "nx_left", "Number of nodes, left"; min = 2),
            nx_right = intg(form, "nx_right", "Number of nodes, right"; min = 2),
            CFL = num(form, "CFL", "CFL number"; min = 0, strict = true),
        )
        kwargs.n in (1, 2, 3) || error("'Geometry' must be 1, 2, or 3.")
        if kwargs.Ri_interface >= kwargs.Ri_total
            error("'Initial interface position [m]' must be smaller than 'Total domain length [m]'.")
        end
        if kwargs.nx_left > kwargs.nx_right
            error("'Number of nodes, left' must be at most 'Number of nodes, right'.")
        end
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:diffusion_couple) do should_stop, report_progress
        res = run_diffusion_couple(; kwargs..., should_stop, report_progress)
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
        t_user, T_user = path(form, "path", "Temperature-time path")
        t_tot_days = num(form, "t_tot_days", "Total time [days]"; min = 0, strict = true)
        # The path's own time column sets the *shape* of the T-t history (relative spacing between
        # nodes); this field is the "easy" total-duration knob the other two modes have via
        # t_tot_years, without having to hand-edit every row just to stretch/compress the run.
        # Scaling by t_user[end] (rather than replacing it outright) preserves that shape exactly.
        t_user = t_user .* (t_tot_days * 24 * 60 * 60 * inv(t_user[end]))
        kwargs = (
            Ri = num(form, "Ri", "Initial interface radius [m]"; min = 0, strict = true),
            n = intg(form, "n", "Geometry"),
            t_user = t_user,
            T_user = T_user,
            P = num(form, "P", "Pressure [Pa]"; min = 0, strict = true),
            D0 = num(form, "D0", "Pre-exponential factor D0, crystal side [m²/s]"; min = 0, strict = true),
            Ea = num(form, "Ea", "Activation energy, crystal side [J/mol]"; min = 0),
            deltaV = num(form, "deltaV", "Activation volume [m³/mol]"; min = 0),
            D0_right = num(form, "D0_right", "Pre-exponential factor D0, melt/fluid side [m²/s]"; min = 0, strict = true),
            Ea_right = num(form, "Ea_right", "Activation energy, melt/fluid side [J/mol]"; min = 0),
            alpha = num(form, "alpha", "Angle to [001] axis [°]"),
            beta = num(form, "beta", "Angle to [010] axis [°]"),
            gamma = num(form, "gamma", "Angle to [100] axis [°]"),
            CompInt = num(form, "CompInt", "Bulk composition of interest [-]"),
            nx_left = intg(form, "nx_left", "Number of nodes, left"; min = 2),
            nx_right = intg(form, "nx_right", "Number of nodes, right"; min = 2),
            CFL = num(form, "CFL", "CFL number"; min = 0, strict = true),
        )
        kwargs.n in (1, 2, 3) || error("'Geometry' must be 1, 2, or 3.")
        if !(0.0 < kwargs.CompInt < 1.0)
            error("'Bulk composition of interest [-]' must be between 0 and 1.")
        end
        if kwargs.nx_left > kwargs.nx_right
            error("'Number of nodes, left' must not be greater than 'Number of nodes, right'.")
        end
    catch e
        return json(Dict("error" => error_message(e)), status = 400)
    end
    run_name = run_name_of(form)
    id = start_job!(:thermo_growth) do should_stop, report_progress
        res = run_thermo_growth(; kwargs..., should_stop, report_progress)
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
# Lets a freshly-loaded page (a new tab, a reload, or reopening the app window after closing
# it - see gui.md's "Accidentally closed the window?") discover a simulation that's still running
# from *before* this page load, so it can reattach its progress bar/Cancel button instead of
# showing a blank idle form with no way to see progress or cancel. Without this, the job keeps
# running server-side exactly as documented, but there's no way back to its UI short of guessing
# its id - effectively a second, subtler version of the dead-end problem the other results-page
# fixes addressed.
@get "/jobs/current" function(req::HTTP.Request)
    job = CURRENT_JOB[]
    (job === nothing || job.status != :running) && return json(Dict("running" => false))
    return json(Dict("running" => true, "id" => job.id, "mode" => string(job.mode),
                      "progress" => job.progress, "started_at" => job.started_at))
end

@get "/jobs/{id}/status" function(req::HTTP.Request, id::String)
    job = get_job(id)
    job === nothing && return json(Dict("status" => "failed", "error" => "Unknown job id."), status = 404)
    if job.status == :running
        return json(Dict("status" => "running", "progress" => job.progress))
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

# ------------------------------------------------------------------
# History cleanup
# ------------------------------------------------------------------
@post "/results/{folder}/delete" function(req::HTTP.Request, folder::String)
    dir = joinpath(RESULTS_DIR, basename(folder))   # basename() blocks path traversal via "../"
    if !isdir(dir)
        return json(Dict("ok" => false, "error" => "Run not found."), status = 404)
    end
    rm(dir; recursive = true)
    return json(Dict("ok" => true))
end

@post "/results/delete-all" function(req::HTTP.Request)
    if isdir(RESULTS_DIR)
        for folder in readdir(RESULTS_DIR)
            dir = joinpath(RESULTS_DIR, folder)
            isdir(dir) && rm(dir; recursive = true)
        end
    end
    return json(Dict("ok" => true))
end

mkpath(RESULTS_DIR)

staticfiles(joinpath(GUI_DIR, "static"), "static")

# `dynamicfiles` re-reads each mounted file's *content* fresh on every request, but the set of
# routes it registers is still a one-time snapshot of the folder taken when it's called (see
# Oxygen's `mountfolder`) - it does not add routes for files created afterwards. Since
# RESULTS_DIR is empty at this point (mkpath above just ensures it exists) and every run's
# files are written well after the server has started, `dynamicfiles(RESULTS_DIR, "results")`
# would 404 on every single result a real session produces, and the "Open full results folder"
# link 404s unconditionally too (there's no index.html in a run's folder for it to serve). These
# two routes serve the same content but resolve the actual path fresh on every request instead.
@get "/results/{folder}/" function(req::HTTP.Request, folder::String)
    dir = joinpath(RESULTS_DIR, basename(folder))   # basename() blocks path traversal via "../"
    if !isdir(dir)
        # A bare 404 body would be a worse dead end than the ones above - no nav at all - and
        # this is reachable by ordinary use (e.g. a stale "Folder" link to a run deleted since
        # the history page was loaded), not just a mistyped URL.
        body = """<p>This run's results folder no longer exists (it may have been deleted).</p><p><a href="/history">&larr; Back to Past runs</a></p>"""
        return html(page("Results not found", "", body); status = 404)
    end
    # download (not a plain link): without it, clicking e.g. plot.pdf navigates the tab to
    # Chrome's built-in PDF viewer - a dead end in the chromeless app window this runs in, with
    # no address bar/back button to escape it (see gui.md's "Accidentally closed the window?").
    items = join(["""<li><a href="/results/$(folder)/$(entry)" download>$(entry)</a></li>"""
                  for entry in sort(readdir(dir))], "\n")
    # Uses the shared page shell (nav bar, styling) rather than a bare fragment: the GUI opens in
    # a chromeless app window with no address bar or back button (see gui.md), so a page reached
    # by clicking through from "Open full results folder" needs its own way back to the rest of
    # the app - a bare <ul> of links would otherwise be a dead end once you're on it.
    body = """<p><a href="/history">&larr; Back to Past runs</a></p><ul>$(items)</ul>"""
    return html(page("Results: $(folder)", "", body))
end
@get "/results/{folder}/{filename}" function(req::HTTP.Request, folder::String, filename::String)
    path = joinpath(RESULTS_DIR, basename(folder), basename(filename))
    isfile(path) || return HTTP.Response(404, "Not found")
    return file(path)
end

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

# Guarded so `GUI/test/` can `include("../app.jl")` to register every route above and test them
# via `Oxygen.internalrequest` without starting a real server - `main()` only runs when this file
# is launched directly (`julia ... GUI/app.jl`, exactly what `run_gui.*` and the docs do), not
# when some other file merely includes it.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
