using Test
include("../examples/Diff_couple_no_interaction.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (no interaction)" begin
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, Cini_l, Cini_r, nmodes_l, nmodes_r,Amp_l, Amp_r, Di, Ri, t = main(plot_sim,verbose)
    xan_l = copy(x_left)
    xan_r = copy(x_right)
    Can_l = sinusoid_profile(Cini_l,nmodes_l,Ri[1],Di[1],t,Amp_l,xan_l)
    Can_r = sinusoid_profile(Cini_r,nmodes_r,Ri[2]-Ri[1],Di[2],t,Amp_r,x_right)
    @test Can_l ≈ C_left rtol = 1e-5
    @test Can_r ≈ C_right rtol = 1e-4
end



