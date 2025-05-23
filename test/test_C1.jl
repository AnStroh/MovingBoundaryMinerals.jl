using Test
using MOBILE, MOBILE.Benchmarks

include("../examples/C1.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (MB) - Rayleigh fractionation" begin
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, n, maxC, KD0 = C1()
    Ray_Fs, Ray_Fl, Ray_Cl, Ray_Cs, Cl_p, phi_solid = rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
    int_Ray   = trapezoidal_integration(Ray_Fs,Ray_Cs)
    int_model = trapezoidal_integration(phi_solid.^n,Cl_p)
    @test int_Ray ≈ int_model rtol = 1e-4
end
