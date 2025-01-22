using Test
include("../examples/B4.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - Rayleigh fractionation" begin
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC = main()
    Ray_Fs, Ray_Fl, Ray_Cl, Ray_Cs, Cl_p, phi_solid = rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
    
    #@test Can ≈ C rtol = 1e-6
end



