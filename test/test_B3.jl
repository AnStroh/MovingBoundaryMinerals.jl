using Test
using MOBILE, MOBILE.Benchmarks
include("../examples/B3.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - Lasaga (1983)" begin
    x_left, x_right, x0, C_left, C_right, C0, Sols_left, Sols_right,Checks, CheckBC, T_pl, t_pl = B3()
    Can1 = first.(sum(Sols_left))
    Can2 = first.(sum(Sols_right))
    @test Can1 ≈ last.(sum(Sols_left)) rtol = 1e-7
    @test Can2 ≈ last.(sum(Sols_right)) rtol = 1e-7
end
