using Test
using MOBILE, MOBILE.Benchmarks

include("../examples/B2.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - time transformation" begin
    x_left, x_right, x0, C_left, C_right, C0, t, D, D0, T0, T, Ea, R = B2()
    C       = [C_left; C_right]
    nterms  = 1000           #Number of terms within the analytical solution (degree of the polynomial)
    Can,xan = crank_time_transformation3(C0,x0,T0,T,Ea,R,D0[1],t,nterms)
    @test Can ≈ C rtol = 1e-2
end
