using Test
using MovingBoundaryMinerals, MovingBoundaryMinerals.Benchmarks
include("../examples/B5.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - Smith (1955)" begin
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC, Di, V_ip, t_tot = B5()
    xan, Can = smith(x_right,C_right,Ri,Di,t_tot,KD0,V_ip,n)
    @test Can ≈ C_right rtol = 1e-4
end
