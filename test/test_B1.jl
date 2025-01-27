using Test
include("../examples/B1.jl")
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux)" begin
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, t, D = main(plot_sim,verbose)
    C       = [C_left; C_right]
    nterms  = 1000           #Number of terms within the analytical solution (degree of the polynomial)
    xan,Can = calc_sinus_sphere(x0,C0,D[1],t,nterms)
    @test Can ≈ C rtol = 1e-2
end



