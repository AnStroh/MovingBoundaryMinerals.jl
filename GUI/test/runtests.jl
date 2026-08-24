using Test

if Threads.nthreads() < 2
    error(
        "GUI tests require multiple threads, since a run is started on a background thread " *
        "via Threads.@spawn (see jobs.jl) - launch with `julia -t auto GUI/test/runtests.jl` " *
        "(or `-t 2`+), exactly like `run_gui.*` launches the app itself.",
    )
end

# Point the app at a throwaway results directory *before* app.jl is included (RESULTS_DIR is
# a `const`, read once at include time), so tests never touch - or delete! - a real GUI/results/
# folder a user might have runs saved in. Cleaned up again at the bottom of this file.
const TEST_RESULTS_DIR = mktempdir()
ENV["MBM_GUI_RESULTS_DIR"] = TEST_RESULTS_DIR
ENV["GKSwstype"] = "100"   # headless GR backend for Plots.jl, same as app.jl sets for itself

# Registers every route (single_crystal.jl, diffusion_couple.jl, thermo_growth.jl, jobs.jl,
# plotting.jl, results.jl, templates/*.jl are all `include`d from here) without starting a real
# HTTP server - `main()` at the bottom of app.jl is guarded to only run when app.jl is launched
# directly, not merely included.
include(joinpath(@__DIR__, "..", "app.jl"))
include(joinpath(@__DIR__, "testutils.jl"))

@testset "MovingBoundaryMinerals.jl GUI" begin
    include("test_pages.jl")
    include("test_progress_reporting.jl")
    include("test_single_crystal.jl")
    include("test_diffusion_couple.jl")
    include("test_thermo_growth.jl")
    include("test_validation_and_job_control.jl")
    include("test_results_management.jl")
end

rm(TEST_RESULTS_DIR; recursive = true, force = true)
