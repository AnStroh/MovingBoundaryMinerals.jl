@testset "Diffusion couple: full run via HTTP" begin
    resp, data = post_form("/diffusion-couple/run", VALID_DIFFUSION_COUPLE_PARAMS)
    @test resp.status == 200
    @test haskey(data, :id)

    result = poll_until_done(data.id)
    @test result.status == "done"
    @test occursin("mass-balance error", result.summary)

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    @test isdir(folder)
    check_result_files(folder, result)

    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["source"] == CITATION_LINE
    @test jl["mode"] == "diffusion_couple"
    # Two-phase modes have separate initial/final grid lengths (regridding, see plotting.jl) -
    # both should still be non-empty.
    @test length(jl["x_initial"]) > 0
    @test length(jl["x_final"]) > 0
end
