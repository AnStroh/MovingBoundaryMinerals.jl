using Test
include("../examples/C1.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - Rayleigh fractionation" begin
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC = main(plot_sim,verbose)
    Ray_Fs, Ray_Fl, Ray_Cl, Ray_Cs, Cl_p, phi_solid = rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
    @test Ray_Fs ≈ C_left rtol = 1e-4
    @test Ray_Fl ≈ C_right rtol = 1e-4
end



