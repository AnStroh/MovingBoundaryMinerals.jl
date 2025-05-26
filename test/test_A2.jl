using Test
using MovingBoundaryMinerals, MovingBoundaryMinerals.Benchmarks
include("../examples/A2.jl")

#Testing-----------------------------------------------------------------------
@testset "Simple diffusion - sphere" begin
    x, C, x0, C0, Di, t, t_tot  = A2()
    nterms  = 1000           #Number of terms within the analytical solution (degree of the polynomial)
    xan,Can = calc_sinus_sphere(x0,C0,Di,t_tot,nterms)
    @test Can ≈ C rtol = 1e-5
end
