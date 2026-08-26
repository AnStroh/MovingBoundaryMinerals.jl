@testset "Thermodynamic growth: full run via HTTP (short duration)" begin
    # See VALID_THERMO_GROWTH_PARAMS in testutils.jl for why the path's total time is cut down
    # so far from the form's own default.
    resp, data = post_form("/thermo-growth/run", VALID_THERMO_GROWTH_PARAMS)
    @test resp.status == 200
    @test haskey(data, :id)

    result = poll_until_done(data.id)
    @test result.status == "done"

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    @test isdir(folder)
    check_result_files(folder, result)

    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["mode"] == "thermo_growth"
    @test jl["parameters"]["CompInt"] == 0.5

    # Thermo-growth-only export: the T-t path actually used, back in the form's own days/°C units.
    @test isfile(joinpath(folder, "path.tab"))
    path_lines = readlines(joinpath(folder, "path.tab"))
    @test occursin("time_days", path_lines[2])
    last_day, last_C = parse.(Float64, split(path_lines[end], '\t'))
    @test last_day ≈ 0.001
    @test last_C == 1350.0
end

@testset "Thermodynamic growth: non-monotonic (D2-style) path" begin
    # A 4-node path with a small mid-path uptick (unlike VALID_THERMO_GROWTH_PARAMS's straight
    # D1-style line), exercising the same non-monotonic-path capability as examples/D2.jl.
    # Kept to a 1 degC uptick rather than a large one: with the tiny default seed crystal (Ri =
    # 0.0001 m) and this drastically shortened total time (see VALID_THERMO_GROWTH_PARAMS's
    # docstring), the interface has barely grown before the uptick, so anything much larger
    # triggers a resorption episode the regridder can't keep up with ("Cannot proceed (Newton
    # Failure)") - a real numerical limit of the tiny/fast regime this test needs for CI speed,
    # not something specific to this uptick. examples/D2.jl's own (much larger, ±10-20 degC)
    # wiggles are fine at its full 30-day scale, where the interface has had days to grow first.
    params = merge(VALID_THERMO_GROWTH_PARAMS, Dict("path" => "0, 1400\n0.0003, 1399\n0.0007, 1401\n0.001, 1350"))
    resp, data = post_form("/thermo-growth/run", params)
    @test resp.status == 200

    result = poll_until_done(data.id)
    @test result.status == "done"

    folder = joinpath(TEST_RESULTS_DIR, result.folder)
    @test isdir(folder)
    check_result_files(folder, result)

    jl = JLD2.load(joinpath(folder, "profiles.jld2"))
    @test jl["parameters"]["t_user"] ≈ [0.0, 0.0003, 0.0007, 0.001] .* (24 * 60 * 60)
    @test jl["parameters"]["T_user"] ≈ [1400.0, 1399.0, 1401.0, 1350.0] .+ 273.0
end
