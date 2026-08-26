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

@testset "Diffusion couple: growth (V_ip != 0) actually moves the interface" begin
    # V_ip and t_tot_years chosen so the interface moves by ~1e-5 m (about 3% of the 3e-4 m of
    # room on the right side of Ri_interface=2e-4/Ri_total=5e-4) - small enough that
    # advect_interface_regrid! won't hit its "Interface moved too fast" guard, but large enough
    # to unambiguously distinguish a real move from numerical noise.
    params = merge(VALID_DIFFUSION_COUPLE_PARAMS, Dict("V_ip" => "3e-11", "t_tot_years" => "0.01"))
    resp, data = post_form("/diffusion-couple/run", params)
    @test resp.status == 200

    result = poll_until_done(data.id)
    @test result.status == "done"

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["parameters"]["V_ip"] == 3e-11

    # With V_ip == 0 (the other test above), advect_interface_regrid!/regrid! are no-ops and the
    # grid keeps its initial 100+150 = 250 nodes throughout. A nonzero V_ip exercises the growth
    # path restored in this change (see run_diffusion_couple) - it adds a node as the interface
    # advects into new territory each step, so the final node count should come out different
    # from the untouched initial grid's, proving that path actually ran rather than being
    # silently skipped.
    @test length(jl["x_initial"]) == 250
    @test length(jl["x_final"]) != length(jl["x_initial"])
end
