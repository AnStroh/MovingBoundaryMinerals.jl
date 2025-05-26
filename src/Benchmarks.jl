module Benchmarks

    using ..MovingBoundaryMinerals
    using LinearAlgebra, SpecialFunctions

    export analytical_sol_step_function, analytical_sol_half_space, calc_sinus_sphere, crank_time_transformation1, crank_time_transformation2, lasaga, rayleigh_fractionation ,smith,crank_time_transformation3

    """
        analytical_sol_step_function(Ri, Di, x, C0, t_tot)

    Calculate the analytical solution for a initiala step function diffusion problem.

    # Arguments
    - `Ri::Vector{Float64}`: Vector containing the initial position of the step.
    - `Di::Float64`: Diffusion coefficient.
    - `x::Vector{Float64}`: Spatial coordinate vector.
    - `C0::Vector{Float64}`: Initial composition vector.
    - `t_tot::Float64`: Total time for which the diffusion is calculated.

    # Returns
    - `Can::Vector{Float64}`: Analytical solution of the composition profile at time `t_tot`.

    # Description
    This function computes the analytical solution for the composition profile of a step function diffusion problem using the error function (`erf`). The solution is based on the initial composition values and the diffusion coefficient over a specified time period.
    """

    function analytical_sol_step_function(Ri,Di,x,C0,t_tot)
        #Analytical solution-----------------------------------------
        x0  = copy(x)
        C   = copy(C0)
        Can = (C[end] + C[1]) * inv(2.0) .+ (C[end] - C[1]) * inv(2.0) * erf.((x0 .- Ri[1]) * inv(2.0) * inv(sqrt(Di * t_tot)))       #Single step function
        return Can
    end

    """
        analytical_sol_half_space(H, Di, x, C0, t_tot)

    Compute the analytical solution for diffusion in a half-space.

    # Arguments
    - `H::Float64`: The position of the half-space boundary.
    - `Di::Float64`: The diffusion coefficient.
    - `x::Vector{Float64}`: The spatial coordinates where the solution is evaluated.
    - `C0::Vector{Float64}`: The initial composition profile.
    - `t_tot::Float64`: The total time for which the solution is computed.

    # Returns
    - `Can::Vector{Float64}`: The composition profile at time `t_tot`.

    # Description
    This function computes the analytical solution for the composition profile in a half-space using the error function (`erf`). The composition profile is calculated based on the initial composition `C0`, the diffusion coefficient `Di`, and the total time `t_tot`.
    """

    function analytical_sol_half_space(H,Di,x,C0,t_tot)
        #Analytical solution-----------------------------------------
        x0  = copy(x)
        C   = copy(C0)
        Can =  C[1] .+ (C[end] - C[1]) * erf.((x0 .- H) * inv(2.0) * inv(sqrt(Di * t_tot)))         #Half-space solution
        return Can
    end

    """
        calc_sinus_sphere(x0, C0, D, tot, nterms)

    Calculate the composition profile on a sphere using the eigenfunction expansion method for a linear diffusion problem with Neumann boundary conditions.

    # Arguments
    - `x0::Vector{Float64}`: Initial positions on the sphere.
    - `C0::Vector{Float64}`: Initial composition profile.
    - `D::Float64`: Diffusion coefficient.
    - `tot::Float64`: Total time for diffusion.
    - `nterms::Int`: Number of terms in the eigenfunction expansion.

    # Returns
    - `x::Vector{Float64}`: Spatial vector.
    - `C::Vector{Float64}`: Composition profile after diffusion.

    # Description
    This function computes the composition profile on a sphere over time using the eigenfunction expansion method. The method is based on the solution provided in Crank (1975): The Mathematics of Diffusion, chapter 6.3. The function iteratively updates the composition profile by summing the contributions from each term in the eigenfunction expansion.
    """

    function calc_sinus_sphere(x0,C0,D,tot,nterms)
        #Routine to test a linear diffusion code with Neumann BC
        #The solution is based on the eigenfuction expansion method
        #Crank (1975): The Mathematics of Diffusion, chapter: 6.3
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

    """
        crank_time_transformation1(C0, x0, T0, T, E, R, D0, t, Cl_i)

    Calculate diffusion with cooling and without growth following Crank (1956). Initial composition profile is a half space solution.

    # Arguments
    - `C0::Vector{Float64}`: Initial composition profile.
    - `x0::Vector{Float64}`: Spatial coordinates.
    - `T0::Float64`: Initial temperature.
    - `T::Float64`: Final temperature.
    - `E::Float64`: Activation energy.
    - `R::Float64`: Universal gas constant.
    - `D0::Float64`: Pre-exponential factor for diffusion coefficient.
    - `t::Float64`: Time duration for the transformation.
    - `Cl_i::Float64`: Initial composition level.

    # Returns
    - `C_Crank::Vector{Float64}`: Composition profile after diffusion.

    # Description
    This function calculates the diffusion profile considering cooling and no growth, based on the method described by Crank (1956). It uses the trapezoidal integration method to compute the integral of the diffusion coefficient over time and applies the error function to determine the composition profile.
    """

    function crank_time_transformation1(C0,x0,T0,T,E,R,D0,t,Cl_i)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        index   = findfirst(C0 -> C0 != Cl_i, C0)
        H       = x0[index]                                                                         #Location of the step
        C_Crank = C0[1] .+ (C0[end] - C0[1]) * erf.((x0 .- H) * inv(2.0 * sqrt(zeta)))              #Half-space solution
        return C_Crank
    end

    """
        crank_time_transformation2(C0, x0, T0, T, E, R, D0, t, Ri)

    Calculate diffusion with cooling and without growth following Crank (1956). Initial Profile: Single step.

    # Arguments
    - `C0::Vector{Float64}`: Initial composition profile.
    - `x0::Vector{Float64}`: Spatial coordinates.
    - `T0::Float64`: Initial temperature.
    - `T::Float64`: Final temperature.
    - `E::Float64`: Activation energy.
    - `R::Float64`: Gas constant.
    - `D0::Float64`: Pre-exponential factor for diffusion coefficient.
    - `t::Float64`: Time duration.
    - `Ri::Vector{Float64}`: Initial radius.

    # Returns
    - `C_Crank::Vector{Float64}`: Composition profile after diffusion.

    # Description
    This function calculates the composition profile after diffusion with cooling and without growth, following the method described by Crank (1956). It uses the trapezoidal integration method to compute the integral of the diffusion coefficient over time and applies the error function to determine the composition profile.
    """
    function crank_time_transformation2(C0,x0,T0,T,E,R,D0,t,Ri)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        C_Crank = (C0[end] + C0[1]) * inv(2.0) .+ (C0[end] - C0[1]) * inv(2.0) * erf.((x0 .- Ri[1]) * inv(2.0 * sqrt(zeta)))       #Single step function
        return C_Crank
    end


    """
        crank_time_transformation3(C0, x0, T0, T, E, R, D0, t, nterms)

    Calculate diffusion with cooling and without growth following Crank (1956). Geometry: Sphere.

    # Arguments
    - `C0::Vector{Float64}`: Initial composition profile.
    - `x0::Vector{Float64}`: Initial spatial coordinates.
    - `T0::Float64`: Initial temperature.
    - `T::Float64`: Final temperature.
    - `E::Float64`: Activation energy.
    - `R::Float64`: Gas constant.
    - `D0::Float64`: Pre-exponential factor for diffusion coefficient.
    - `t::Float64`: Time duration for diffusion.
    - `nterms::Int`: Number of terms in the series expansion.

    # Returns
    - `C_Crank::Vector{Float64}`: Composition profile after diffusion.
    - `xan::Vector{Float64}`: Spatial coordinates.

    # Description
    This function calculates the diffusion in a sphere with cooling and without growth using the analytical solution provided by Crank (1975). The diffusion coefficient is temperature-dependent and is integrated over time using trapezoidal integration. The composition profile is updated iteratively using a series expansion.
    """
    function crank_time_transformation3(C0,x0,T0,T,E,R,D0,t,nterms)
        #Calculate diffusion with cooling and without growth following Crank (1956)
        Tan     = LinRange(T0,T,1000)
        tan     = LinRange(0,t,1000)
        Dan     = D0.*exp.(-E.*inv.(R .* Tan))
        zeta    = trapezoidal_integration(tan,Dan)
        #Analytical solution diffusion in a sphere, Crank (1975)
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

    """
        rayleigh_fractionation(x_left, C_left, Ri0, Ri, C0_r, KD0, n)

    Compute the Rayleigh fractionation due to growth without diffusion.

    # Arguments
    - `x_left::Vector{Float64}`: Positions on the left side.
    - `C_left::Vector{Float64}`: Concentrations on the left side.
    - `Ri0::Vector{Float64}`: Initial radii.
    - `Ri::Vector{Float64}`: Actuel radii.
    - `C0_r::Vector{Float64}`: Initial concentrations.
    - `KD0::Float64`: Distribution coefficient.
    - `n::Int`: Geometry factor.

    # Returns
    - `Ray_Fs::Vector{Float64}`: Rayleigh's solid fraction.
    - `Ray_Fl::Vector{Float64}`: Rayleigh's liquid fraction.
    - `Ray_Cl::Vector{Float64}`: Rayleigh's solution for the liquid composition.
    - `Ray_Cs::Vector{Float64}`: Rayleigh's solution for the solid composition.
    - `C_l_p::Vector{Float64}`: Post-processed concentrations on the left side.
    - `phi_solid::Vector{Float64}`: Solid fraction within the model.

    # References
    - Albarède, F. (2003) Geochemistry: An Introduction, Cambridge Press, 262 pp.
    """
    function rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
        #=Check Rayleigh Limits -> fractionation due to growth, no diffusion
        Albarède, F. (2003) Geochemistry: An Introduction, Cambridge Press, 262 pp.
        Cliquid/C_tot = F^(KD-1) F: melt fraction    => C_tot = 1 (closed system)
        KD = Cs/Cl =#
        #Fractions Rayleigh Limit------------------------------------
        Ray_Fs = LinRange(Ri0[1] *inv(Ri0[2]),Ri[1] * inv(Ri0[2]),1000) .^ n                        #Solid fraction; solid = solid/liquid fraction
        Ray_Fl = 1.0 .- Ray_Fs                                                                      #Liquid fraction
        Ray_Cl = Ray_Fl .^ (KD0 - 1.0) .* C0_r[end]                                                 #Liquid composition
        Ray_Cs = Ray_Cl .* KD0                                                                      #Solid composition
        #Post processing---------------------------------------------
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

    """
        smith(x_right, C_right, Ri, Di, t_tot, KD, V_ip, n)

    Analytical solution for the composition of a liquid in growth medium against an advancing planar crystal interface based on Smith et al. (1955).

    # Arguments
    - `x_right::Array{Float64}`: Array of positions on the right side.
    - `C_right::Array{Float64}`: Array of concentrations on the right side.
    - `Ri::Array{Float64}`: Actuel radii.
    - `Di::Array{Float64}`: Array of diffusion coefficients.
    - `t_tot::Float64`: Total time.
    - `KD::Float64`: Distribution coefficient.
    - `V_ip::Float64`: Interface velocity.
    - `n::Int`: Geometry factor.

    # Returns
    - `xan::Array{Float64}`: Array of positions (analytical solution).
    - `Can::Array{Float64}`: Array of calculated concentrations (analytical solution).

    # Notes
    This function calculates the composition of a liquid in a growth medium against an advancing planar crystal interface using an analytical solution derived by Smith (1955). The calculation is performed only if `n == 1`.
    """
    function smith(x_right,C_right,Ri,Di,t_tot,KD,V_ip,n)
        #Smith et al. (1955): analytical solution for the composition of a liquid in growth medium against an advancing planar crystal interface
        if n == 1
            xan  = x_right .- Ri[1]
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
            Can  = Can * Can0                                                                       #Solution for liquid
        end
        return xan, Can
    end
end
