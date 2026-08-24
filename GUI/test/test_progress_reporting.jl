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
        run_thermo_growth(; t_tot_days = 0.001, report_progress = p -> push!(seen, p))
        @test !isempty(seen)
        @test issorted(seen)
        @test seen[end] ≈ 1.0
    end
end
