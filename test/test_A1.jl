using Test
include("../examples/A1.jl")
#Testing-----------------------------------------------------------------------
@testset "Simple diffusion planar (1D)" begin
    x, C, x0, C0, Di, t, t_tot, Cini, nmodes, Amp, L  = main()
    xan = copy(x)
    Can = sinusoid_profile(Cini,nmodes,L,Di,t,Amp,xan)
    @test Can ≈ C rtol = 1e-6
end

