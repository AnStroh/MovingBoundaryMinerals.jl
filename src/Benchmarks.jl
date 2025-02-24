#= Utilized functions in the framework of the Paper:

    by A.Stroh and E. Moulas
    doi:
    Version: 1.0
=#
module Benchmarks 
    
    using ..Diff_Coupled
    using LinearAlgebra, BenchmarkTools, SpecialFunctions
    
    export analytical_sol_step_function, analytical_sol_half_space, calc_sinus_sphere, crank_time_transformation1, crank_time_transformation2, lasaga, rayleigh_fractionation ,smith,crank_time_transformation3

    

    function analytical_sol_step_function(Ri,Di,x,C0,t_tot)
        #Analytical solution-------------------------------------------------------
        x0  = copy(x)
        C   = copy(C0)
        Can = (C[end] + C[1]) * inv(2.0) .+ (C[end] - C[1]) * inv(2.0) * erf.((x0 .- Ri[1]) * inv(2.0) * inv(sqrt(Di * t_tot)))       #single Step function
        return Can   
    end
    function analytical_sol_half_space(H,Di,x,C0,t_tot)
        #Analytical solution-------------------------------------------------------
        x0  = copy(x)
        C   = copy(C0)
        Can =  C[1] .+ (C[end] - C[1]) * erf.((x0 .- H) * inv(2.0) * inv(sqrt(Di * t_tot)))       #Half-space solution
        return Can   
    end

    function calc_sinus_sphere(x0,C0,D,tot,nterms)
        #Routine to test a linear diffusion code with Neumann BC
        #The solution is based on the eigenfuction expansion method
        # Crank (1975): The Mathematics of Diffusion, chapter: 6.3 
        C = copy(C0)
        x = copy(x0)
        C1   = C0[1] 
        Cs  = C0[end]
        L   = x0[end]
        for nt = 1:nterms
            C   .= C  .+ 2.0 * L * inv.(pi .* x) .* ((-1.0) .^ nt) .* inv(nt) .* sin.(nt .* pi .* x .* inv(L)) .* exp(-D * nt ^ 2 * pi ^ 2 * tot * inv(L ^ 2))
            C1   = C1 .+ 2.0 * (((-1.0) .^ nt) .* exp(-D * nt ^ 2 * pi ^ 2 * tot * inv(L ^ 2)))
        end
        C      = C  .+ 1.0
        C[1]   = C1 .+ 1.0
        C[end] = Cs
        return x,C
    end

    function crank_time_transformation1(C0,x0,T0,T,E,R,D0,t,Cl_i)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Myr2sec = 1e6 * 60 * 60 * 24 * 365.25
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        index   = findfirst(C0 -> C0 != Cl_i, C0)
        H       = x0[index]                                                             #Location of the step                                                              
        C_Crank = C0[1] .+ (C0[end] - C0[1]) * erf.((x0 .- H) * inv(2.0 * sqrt(zeta)))       #Half-space solution       
        return C_Crank
    end
    function crank_time_transformation2(C0,x0,T0,T,E,R,D0,t,Ri)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Myr2sec = 1e6 * 60 * 60 * 24 * 365.25
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        C_Crank = (C0[end] + C0[1]) * inv(2.0) .+ (C0[end] - C0[1]) * inv(2.0) * erf.((x0 .- Ri[1]) * inv(2.0 * sqrt(zeta)))       #single Step function
        return C_Crank
    end
    function crank_time_transformation3(C0,x0,T0,T,E,R,D0,t,nterms)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Myr2sec = 1e6 * 60 * 60 * 24 * 365.25
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        #Analytical solution Diffusion in a sphere, Crank (1975)
        C_Crank = copy(C0)
        xan     = copy(x0)
        C1      = C0[1] 
        Cs      = C0[end]
        L       = x0[end]
        for nt = 1:nterms
            C_Crank   .= C_Crank  .+ 2.0 * L * inv.(pi .* xan)  .* ((-1.0) .^ nt) .*inv(nt) .* sin.(nt .* pi .* xan .* inv(L)) .* exp(-zeta * nt ^ 2 * pi ^ 2 * inv(L ^ 2))
            C1   = C1 .+ 2.0 * (((-1.0) .^ nt) .* exp(-zeta * nt ^ 2 * pi ^ 2 * inv(L ^ 2)))
        end
        C_Crank      = C_Crank  .+ 1.0
        C_Crank[1]   = C1 .+ 1.0
        C_Crank[end] = Cs
        return C_Crank, xan
    end

    function rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
        #=Check Rayleigh Limits -> fractionation due to growth, no diffusion
        Albarède, F. (2003) Geochemistry: An Introduction, Cambridge Press, 262 pp.
        Cliquid/C_tot = F^(KD-1) F: melt fraction    => C_tot = 1 (closed system)
        KD = Cs/Cl =#
        #Reshape---------------------------------------------------------------
        #C_left = copy(C_left)'
        #x_left = copy(x_left)'
        #Fractions Rayleigh Limit--------------------------------------------------
        Ray_Fs = LinRange(Ri0[1] *inv(Ri0[2]),Ri[1] * inv(Ri0[2]),1000) .^ n        #Solid fraction; solid = solid/melt fraction (liquid)
        Ray_Fl = 1.0 .- Ray_Fs                                             #Melt fraction
        Ray_Cl = Ray_Fl .^ (KD0 - 1.0) .* C0_r[end]                        #Liquid concentration
        Ray_Cs = Ray_Cl .* KD0                                             #Solid concentration
        #Post processing----------------------------------------------------------
        indi   = findall(x -> x > Ri0[1], x_left)
        x_l_p  = zeros(size(x_left))
        x_l_p  = copy(x_left[LinearIndices(indi)])
        x_l_p  = [Ri0[1] x_l_p']
        phi_solid = x_l_p * inv(Ri0[2])
        C_l_p  = zeros(size(C_left))
        C_l_p  = copy(C_left[LinearIndices(indi)])
        C_bc_p = linear_interpolation_1D(x_left,C_left,Ri0[1])
        C_l_p  = [C_bc_p C_l_p'] 
        return Ray_Fs, Ray_Fl, Ray_Cl, Ray_Cs, C_l_p, phi_solid
    end

    function smith(x_right,C_right,Ri,Di,t_tot,KD,V_ip,n)
        #Smith (1955): analytical solution for the concentration of a liquid in growth medium against an advancing planar crystal interface
        if n == 1
            println("dffg")
            xan  = x_right .- Ri[1]
            @show xan
            Dan  = Di[2]
            tan  = t_tot
            Kan  = KD
            R    = V_ip
            ksi  = 0.5 .* (Dan .* tan) .^ (-0.5)
            Can0 = C_right[end]
            Can  = 1.0 .+ ((1.0 .- Kan) .* inv(2.0 * Kan)) .* exp.(- R .* inv(Dan) .* xan) .* erfc.(ksi .* (xan .- R .* tan)) .- 
                    0.5 .* erfc.(ksi .* (xan .+ R .* tan)) .+ ((1.0 .- Kan) .* inv(2.0)) .* (1.0 .* inv(1.0 .- Kan) .- 1.0 .* inv(Kan)) .* 
                    exp.(- (1.0 .- Kan) .* (R .* inv(Dan)) .* (xan .+ Kan .* R .* tan)) .* 
                    erfc.(ksi .* (xan .+ (2.0 .* Kan .- 1.0) .* R .* tan))
            Can  = Can * Can0           #solution for liquid
        end
        return xan, Can
    end
end