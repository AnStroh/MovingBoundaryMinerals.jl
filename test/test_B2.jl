using Test
include("../examples/B2.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (Lasaga)" begin
    x_left, x_right, x0, C_left, C_right, C0, Sols_left, Sols_right,Checks, CheckBC, T_pl, t_pl = main()
    Can1 = first.(Sols_left)
    Can2 = first.(Sols_right)
    #@test Can ≈ C rtol = 1e-6
end



