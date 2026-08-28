@testset "report_progress callback" begin
    # Called directly (not through HTTP) so this is fast and immune to the HTTP-polling timing
    # issue where a fast run can finish between two 1-second polls without ever being observed
    # mid-flight - see the progress-bar work in the docs history for why that matters here.
    @testset "run_single_crystal" begin
        seen = Float64[]
        run_single_crystal(; report_progress = p -> push!(seen, p))
        @test !isempty(seen)
        @test issorted(seen)
        @test seen[end] ≈ 1.0
    end

    @testset "run_diffusion_couple" begin
        seen = Float64[]
        run_diffusion_couple(; report_progress = p -> push!(seen, p))
        @test !isempty(seen)
        @test issorted(seen)
        @test seen[end] ≈ 1.0
    end

    @testset "run_thermo_growth" begin
        seen = Float64[]
        # Ri pinned to the pre-2026-08-28 default (0.0001): this test's own path is deliberately
        # extreme for CI speed (50 degC in 0.001 days) - fine at the tiny original Ri, but at the
        # new default (0.0005, chosen for robustness on realistic non-monotonic paths) it grows
        # the interface fast enough to trip the "Interface moved too fast" guard.
        run_thermo_growth(; Ri = 0.0001, t_user = [0.0, 0.001 * 24 * 60 * 60], T_user = [1400.0, 1350.0] .+ 273.0,
                           report_progress = p -> push!(seen, p))
        @test !isempty(seen)
        @test issorted(seen)
        @test seen[end] ≈ 1.0
    end
end
