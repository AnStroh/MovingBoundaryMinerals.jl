using Test
using MovingBoundaryMinerals
# using LinearAlgebra, LaTeXStrings, SparseArrays
include("../examples/A1.jl")
#Testing-----------------------------------------------------------------------
@testset "Simple diffusion - planar" begin
    x, C, x0, C0, Di, t, t_tot, Cini, nmodes, Amp, L  = A1()
    xan = copy(x)
    Can = sinusoid_profile(Cini,nmodes,L,Di,t,Amp,xan)
    @test Can ≈ C rtol = 1e-4
end
