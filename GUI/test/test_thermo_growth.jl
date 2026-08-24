@testset "Thermodynamic growth: full run via HTTP (short duration)" begin
    # See VALID_THERMO_GROWTH_PARAMS in testutils.jl for why t_tot_days is cut down so far
    # from the form's own default.
    resp, data = post_form("/thermo-growth/run", VALID_THERMO_GROWTH_PARAMS)
    @test resp.status == 200
    @test haskey(data, :id)

    result = poll_until_done(data.id)
    @test result.status == "done"

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    @test isdir(folder)
    check_result_files(folder)

    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["mode"] == "thermo_growth"
    @test jl["parameters"]["CompInt"] == 0.5
end
