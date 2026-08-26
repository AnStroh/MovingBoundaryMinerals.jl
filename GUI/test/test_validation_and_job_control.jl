@testset "Input validation" begin
    @testset "non-numeric field -> 400 with a specific message" begin
        bad = merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("L" => "not-a-number"))
        resp, data = post_form("/single-crystal/run", bad)
        @test resp.status == 400
        @test occursin("Domain length", data.error)
    end

    @testset "diffusion couple: interface position past total length -> 400" begin
        bad = merge(VALID_DIFFUSION_COUPLE_PARAMS, Dict("Ri_interface" => "0.001", "Ri_total" => "0.0005"))
        resp, data = post_form("/diffusion-couple/run", bad)
        @test resp.status == 400
        @test occursin("Initial interface position", data.error)
    end

    @testset "thermo growth: CompInt out of (0, 1) -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("CompInt" => "1.5"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("CompInt", data.error) || occursin("composition of interest", lowercase(data.error))
    end

    @testset "thermo growth: path with only 1 point -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("path" => "0, 1400"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("at least 2 points", data.error)
    end

    @testset "thermo growth: path not starting at t=0 -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("path" => "1, 1400\n30, 1350"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("first time value must be 0", data.error)
    end

    @testset "thermo growth: path with non-increasing time -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("path" => "0, 1400\n10, 1420\n5, 1350"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("strictly increase", data.error)
    end

    @testset "thermo growth: malformed path line -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("path" => "0, 1400\nnot-a-pair"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("time, temperature", data.error)
    end

    # Found by manually bypassing the client-side HTML5 `min="0"` on these fields (form.noValidate
    # = true) and submitting anyway: with no server-side check, a negative domain length was
    # silently accepted and completed as a mirrored-but-plausible-looking domain (x from -L to 0)
    # with nothing telling the user their input was invalid. num()/intg() now take an optional
    # `min`/`strict` bound (see app.jl) so the same physically-invalid values are also rejected
    # server-side, not just nudged away from client-side - the only path a raw API request (or a
    # browser with JS disabled) actually goes through.
    @testset "single crystal: negative domain length -> 400" begin
        bad = merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("L" => "-0.001"))
        resp, data = post_form("/single-crystal/run", bad)
        @test resp.status == 400
        @test occursin("Domain length", data.error)
    end

    @testset "single crystal: CFL = 0 -> 400 (would never advance dt otherwise)" begin
        bad = merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("CFL" => "0"))
        resp, data = post_form("/single-crystal/run", bad)
        @test resp.status == 400
        @test occursin("CFL", data.error)
    end

    @testset "single crystal: nx = 1 -> 400 (below the minimum viable grid)" begin
        bad = merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("nx" => "1"))
        resp, data = post_form("/single-crystal/run", bad)
        @test resp.status == 400
        @test occursin("Number of nodes", data.error)
    end

    @testset "single crystal: invalid geometry -> 400" begin
        bad = merge(VALID_SINGLE_CRYSTAL_PARAMS, Dict("n" => "5"))
        resp, data = post_form("/single-crystal/run", bad)
        @test resp.status == 400
        @test occursin("Geometry", data.error)
    end

    @testset "diffusion couple: negative diffusion coefficient -> 400" begin
        bad = merge(VALID_DIFFUSION_COUPLE_PARAMS, Dict("D0_left" => "-2.75e-6"))
        resp, data = post_form("/diffusion-couple/run", bad)
        @test resp.status == 400
        @test occursin("Pre-exponential factor", data.error)
    end

    @testset "diffusion couple: zero distribution coefficient -> 400" begin
        bad = merge(VALID_DIFFUSION_COUPLE_PARAMS, Dict("KD_start" => "0"))
        resp, data = post_form("/diffusion-couple/run", bad)
        @test resp.status == 400
        @test occursin("Distribution coefficient", data.error)
    end

    @testset "thermo growth: negative pressure -> 400" begin
        bad = merge(VALID_THERMO_GROWTH_PARAMS, Dict("P" => "-1.0e6"))
        resp, data = post_form("/thermo-growth/run", bad)
        @test resp.status == 400
        @test occursin("Pressure", data.error)
    end
end

@testset "Concurrent run rejected (409)" begin
    resp1, data1 = post_form("/single-crystal/run", VALID_SINGLE_CRYSTAL_PARAMS)
    @test resp1.status == 200
    # Issued immediately, before the first job could possibly have finished (job creation and
    # the actual FEM solve are on different threads - see jobs.jl), so this reliably lands
    # while the first job is still :running.
    resp2, data2 = post_form("/single-crystal/run", VALID_SINGLE_CRYSTAL_PARAMS)
    @test resp2.status == 409
    poll_until_done(data1.id)   # drain it so it can't leak into a later test
end

@testset "Cancel a running job" begin
    # t_tot_years is bumped far above VALID_DIFFUSION_COUPLE_PARAMS's default (1000) so this
    # reliably takes several seconds even fully JIT-warmed, giving the cancel request a real
    # window to land before the job would otherwise finish on its own.
    slow = merge(VALID_DIFFUSION_COUPLE_PARAMS, Dict("t_tot_years" => "50000.0"))
    resp, data = post_form("/diffusion-couple/run", slow)
    @test resp.status == 200

    cresp, cdata = post_form("/jobs/$(data.id)/cancel", Dict())
    @test cresp.status == 200
    @test cdata.ok == true

    result = poll_until_done(data.id)
    @test result.status == "cancelled"
    # Cancelled runs save nothing (see save_run_outputs's docstring) - nothing new should have
    # appeared under the results dir for this job.
    @test !haskey(result, :folder)
end

@testset "Cancelling with no job running" begin
    cresp, cdata = post_form("/jobs/not-a-real-id/cancel", Dict())
    @test cresp.status == 200
    @test cdata.ok == false
end

@testset "Status of an unknown job id" begin
    resp, data = get_json("/jobs/not-a-real-id/status")
    @test resp.status == 404
end
