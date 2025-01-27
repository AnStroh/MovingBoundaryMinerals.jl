using Test
include("../examples/B2.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (Lasaga)" begin
    adapt_dt = true
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, Sols_left, Sols_right,Checks, CheckBC, T_pl, t_pl = main(adapt_dt,plot_sim,verbose)
    Can1 = first.(Sols_left)
    Can2 = first.(Sols_right)
    @test Can1 ≈ last.(Sols_left) rtol = 1e-6
    @test Can1 ≈ last.(Sols_right) rtol = 1e-6
end



