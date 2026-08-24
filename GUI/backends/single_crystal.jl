using MovingBoundaryMinerals, MovingBoundaryMinerals.Benchmarks
using LinearAlgebra, SparseArrays

"""
    run_single_crystal(; kwargs...)

GUI backend for the single-crystal diffusion mode, adapted from `examples/Simple_Diff.jl`.
Every hardcoded physics/numerics value in the original example is exposed here as a keyword
argument with the same default, so calling with no arguments reproduces the example's result.
Never plots internally (unlike the original); always returns the raw arrays for the caller to plot.

# Returns
`(; x, C, x0, C0, D, t, t_tot, MB_Error)`
"""
function run_single_crystal(;
        Di          = -1.0,       #Diffusion coefficient in [m^2/s]; -1.0 = calculate D via Arrhenius from D0/Ea1
        D0          = 2.75e-6,    #Pre-exponential factor in [m^2/s]
        L           = 0.001,      #Length of the domain in [m]
        Cinf        = 0.0,        #Composition at infinity in [-]
        Cstart      = 4.0,        #Initial composition in [-]
        rho         = 2700.0,     #Density in [kg/m^3]
        Ea1         = 292880.0,   #Activation energy in [J/mol]
        Tstart_C    = 1000.0,     #Starting temperature in [degC]
        Tstop_C     = 650.0,      #End temperature in [degC]
        t_tot_years = 100.0,      #Total time in [years]
        n           = 1,          #Geometry; 1: planar, 2: cylindrical, 3: spherical
        nx          = 100,        #Number of nodes
        CFL         = 0.99,       #CFL number for time step calculation
        should_stop::Function = () -> false,   #Checked every iteration; throw(SimulationCancelled()) if true
        report_progress::Function = (_) -> nothing,   #Called every iteration with t/t_tot in [0, 1]
    )
    R       = 8.314472
    t_tot   = t_tot_years * 365.25 * 24 * 60 * 60
    BCout   = n == 1 ? [1, 1] : [0, 1]   #Left BC must be Neumann for cylindrical/spherical geometry

    #Domain-------------------------------------------------------------
    dx    = L * inv(nx - 1)
    x     = [0:dx:L;]

    #Initial condition----------------------------------------------------
    t     = 0.0
    it    = 0
    C     = Cinf * ones(nx, 1)
    C[1]  = Cstart
    C0    = copy(C)
    x0    = copy(x)

    #History dependent parameters-------------------------------------------
    T_ar    = LinRange(Tstart_C + 273.15, Tstop_C + 273.15, 1000)
    t_ar    = LinRange(0.0, t_tot, 1000)

    dt_diff = t_ar[2:end] .- t_ar[1:end-1]

    #Preallocate variables------------------------------------------------
    Co       = zeros(size(C))
    dt       = 0.0
    L_g      = spzeros(length(x), length(x))
    Mass     = Float64[]
    nels     = length(x) - 1
    R_g      = zeros(length(x), 1)
    MB_Error = Float64[]

    #Calculate initial Ds, T-----------------------------------------------
    D, T  = update_t_dependent_param_simple!(D0, Di, Ea1, R, T_ar, t_ar, t, t_tot)
    Mass0 = calc_mass_vol_simple_diff(x, C, n, rho)

    #Input validation-------------------------------------------------------
    if BCout[1] != 0 && (n == 3 || n == 2)
        error("The code is only valid for cylindrical/spherical geometry, where the left outer BC has Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T != T_ar[1]
        error("Initial temperature must be equal to the first value in the temperature array.")
    end

    #Time loop--------------------------------------------------------------
    while t < t_tot
        should_stop() && throw(SimulationCancelled())
        Co = copy(C)
        dt = calculate_dt(D, dx, CFL)
        t, dt, it = update_time!(t, dt, it, t_tot)
        report_progress(t * inv(t_tot))
        D, T = update_t_dependent_param_simple!(D0, Di, Ea1, R, T_ar, t_ar, t, t_tot)
        L_g, R_g = fill_matrix!(C, x, D, dt, n, nels)
        ScF = sum(diag(L_g)) * inv(length(diag(L_g)))
        L_g, R_g = set_outer_bc!(BCout, L_g, R_g, Co[1], Co[end], ScF)
        C = L_g \ R_g
        redirect_stdout(devnull) do
            ErrM = calc_mass_vol_simple_diff(x, C, n, rho)
            push!(MB_Error, ErrM)
        end
    end
    Mass = calc_mass_vol_simple_diff(x, C, n, rho)
    redirect_stdout(devnull) do
        calc_mass_err(Mass, Mass0)
    end
    return (; x, C, x0, C0, D, t, t_tot, MB_Error)
end
