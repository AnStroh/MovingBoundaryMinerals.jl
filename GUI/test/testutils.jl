using Oxygen, HTTP, JSON3, Test

"""POST a form-encoded request to `path` with the given `Dict` of params (values are
stringified and URI-escaped), simulated in-process via `Oxygen.internalrequest` - no real
HTTP server involved. Returns `(response, parsed_json_body)`."""
function post_form(path::String, params::AbstractDict)
    body = join(["$(k)=$(HTTP.escapeuri(string(v)))" for (k, v) in params], "&")
    req = Oxygen.Request("POST", path, ["Content-Type" => "application/x-www-form-urlencoded"], Vector{UInt8}(body))
    resp = Oxygen.internalrequest(req)
    return resp, JSON3.read(String(resp.body))
end

"""GET `path` in-process and return `(response, parsed_json_body)`."""
function get_json(path::String)
    resp = Oxygen.internalrequest(Oxygen.Request("GET", path, [], UInt8[]))
    return resp, JSON3.read(String(resp.body))
end

"""GET `path` in-process and return the raw response (for HTML pages, not JSON endpoints)."""
get_page(path::String) = Oxygen.internalrequest(Oxygen.Request("GET", path, [], UInt8[]))

"""
    poll_until_done(id; timeout_s=120) -> status_payload

Polls `/jobs/{id}/status` every 0.1s until it's no longer `"running"`. The default timeout is
generous on purpose: the *first* run of each backend mode in a fresh Julia process pays a
one-off JIT compilation cost that alone can take 20+ seconds - exactly the situation in CI.
Fails the test (via `@test false`) and returns the last-seen payload if the deadline passes
without the job finishing, rather than hanging forever.
"""
function poll_until_done(id::AbstractString; timeout_s::Real = 120)
    deadline = time() + timeout_s
    local data
    while time() < deadline
        _, data = get_json("/jobs/$(id)/status")
        data.status == "running" || return data
        sleep(0.1)
    end
    @test false   # timed out - surfaces as a normal test failure with a useful location
    return data
end

"""A known-good parameter set for `/single-crystal/run`, reused across test files so each one
doesn't have to repeat every field name."""
const VALID_SINGLE_CRYSTAL_PARAMS = Dict(
    "D0" => "2.75e-6", "Ea1" => "292880.0", "L" => "0.001", "Cstart" => "4.0", "Cinf" => "0.0",
    "n" => "1", "Tstart_C" => "1000.0", "Tstop_C" => "650.0",
    "t_tot_years" => "100.0", "nx" => "100", "CFL" => "0.99",
)

"""A known-good parameter set for `/diffusion-couple/run`."""
const VALID_DIFFUSION_COUPLE_PARAMS = Dict(
    "D0_left" => "2.75e-6", "D0_right" => "3.9e-7", "Ea_left" => "292879.6767", "Ea_right" => "360660.4018",
    "Ri_interface" => "0.0002", "Ri_total" => "0.0005", "Cl_i" => "0.6", "Cr_i" => "0.3",
    "KD_start" => "1.0", "KD_end" => "0.7", "V_ip" => "0.0", "n" => "3", "Tstart_C" => "1000.0", "Tstop_C" => "700.0",
    "t_tot_years" => "1000.0", "nx_left" => "100", "nx_right" => "150", "CFL" => "0.5",
)

"""A known-good parameter set for `/thermo-growth/run`, with the path's total time cut down
drastically from the form's own default (30 days) - the underlying physics needs a great many
small time steps regardless of the requested duration, so 30 days would take minutes even fully
JIT-warmed. 0.001 days is still a real (if tiny) simulation, just a fast one, suitable for CI."""
const VALID_THERMO_GROWTH_PARAMS = Dict(
    "Ri" => "0.0001", "n" => "1", "path" => "0, 1400\n0.001, 1350", "t_tot_days" => "0.001",
    "P" => "1.0e6", "D0" => "5.38e-9", "Ea" => "226000.0", "deltaV" => "7e-6",
    "D0_right" => "0.0003634023264950478", "Ea_right" => "218022.084784",
    "alpha" => "0.0", "beta" => "90.0", "gamma" => "90.0", "CompInt" => "0.5",
    "nx_left" => "50", "nx_right" => "50", "CFL" => "50.0",
)

"""Asserts every expected output file exists, is non-empty, and (for the text/bundled formats)
carries the small citation line - see `plotting.jl`'s `CITATION_LINE` and `gui.md`.

`result` is the `done` job-status payload (has `png_url`/`pdf_url`/`data_url`/`xlsx_url`/
`jld2_url`/`folder_url`) - every one of these is also fetched via `get_page` and checked for a
200 status, since these routes (results/{folder}/{filename} and results/{folder}/) resolve the
path fresh on every request rather than off a fixed set of routes registered once at server
startup (see the comment above them in app.jl for why that distinction matters here): a disk
check alone can't catch a route that 404s despite the file being right there on disk."""
function check_result_files(folder::String, result)
    for f in ("plot.png", "plot.pdf", "profile_initial.tab", "profile_final.tab",
              "inputs.toml", "profiles.xlsx", "profiles.jld2")
        path = joinpath(folder, f)
        @test isfile(path)
        @test filesize(path) > 0
    end
    @test occursin(CITATION_LINE, read(joinpath(folder, "inputs.toml"), String))
    @test occursin(CITATION_LINE, read(joinpath(folder, "profile_initial.tab"), String))
    @test occursin(CITATION_LINE, read(joinpath(folder, "profile_final.tab"), String))

    for url in (result.png_url, result.pdf_url, result.data_url, result.xlsx_url, result.jld2_url, result.folder_url)
        resp = get_page(url)
        @test resp.status == 200
        @test !isempty(resp.body)
    end
end
