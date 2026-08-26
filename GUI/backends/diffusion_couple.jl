using MovingBoundaryMinerals
using LinearAlgebra, SparseArrays

"""
    run_diffusion_couple(; kwargs...)

GUI backend for the diffusion-couple (moving interface, flux-balance) mode, adapted from
`examples/Diff_couple_Flux.jl` (`V_ip = 0`, the default) and `examples/Diff_couple_Flux_growth.jl`
(`V_ip != 0`) - the same function covers both, exactly as the two example scripts differ only in
that one value and the resulting grid needing to be re-advected/regridded every step once the
interface actually moves. Every hardcoded physics/numerics value that a non-programmer would
plausibly want to change is exposed here as a keyword argument with the same default, so calling
with no arguments reproduces `Diff_couple_Flux.jl`'s result. Advanced/rarely changed settings (the
grid refinement strategy/grading, as opposed to its resolution) are kept fixed at the example's
defaults. Never plots internally; always returns the raw arrays for the caller to plot.

# Returns
`(; x_left, x_right, dx1, dx2, x0, res, Ri, C_left, C_right, C0, MB_Error, t_tot)`
"""
function run_diffusion_couple(;
        D0_left   = 2.75e-6,      #Pre-exponential factor, left side, in [m^2/s]
        D0_right  = 3.9e-7,       #Pre-exponential factor, right side, in [m^2/s]
        Ea_left   = 292879.6767,  #Activation energy, left side, in [J/mol]
        Ea_right  = 360660.4018,  #Activation energy, right side, in [J/mol]
        Ri_interface = 0.0002,    #Initial interface position in [m]
        Ri_total     = 0.0005,    #Total domain length in [m]
        Cl_i      = 0.6,          #Initial composition, left side, in [-]
        Cr_i      = 0.3,          #Initial composition, right side, in [-]
        KD_start  = 1.0,          #Distribution coefficient at t=0
        KD_end    = 0.7,          #Distribution coefficient at t=t_tot
        Tstart_C  = 1000.0,       #Starting temperature in [degC]
        Tstop_C   = 700.0,        #End temperature in [degC]
        t_tot_years = 1000.0,     #Total time in [years]
        V_ip      = 0.0,          #Interface velocity in [m/s]; 0 = stationary interface, >0 = growth (left phase consumes right), <0 = resorption
        n         = 3,            #Geometry; 1: planar, 2: cylindrical, 3: spherical
        nx_left   = 100,          #Number of nodes, left side
        nx_right  = 150,          #Number of nodes, right side
        CFL       = 0.5,          #CFL number for time step calculation
        should_stop::Function = () -> false,   #Checked every iteration; throw(SimulationCancelled()) if true
        report_progress::Function = (_) -> nothing,   #Called every iteration with t/t_tot in [0, 1]
    )
    Di          = [-1.0 -1.0;]
    D0          = [D0_left D0_right;]
    rho         = [1.0  1.0;]
    Ri          = [Ri_interface Ri_total;]
    R           = 8.314472
    t_tot       = t_tot_years * 365.25 * 24 * 60 * 60
    RefineMethod = 1
    verbose     = false

    #History dependent parameters---------------------------------
    KD_ar   = LinRange(KD_start, KD_end, 1000)
    t_ar    = LinRange(0.0, t_tot, 1000)
    T_ar    = LinRange(Tstart_C + 273.15, Tstop_C + 273.15, 1000)

    #Numerics-------------------------------------------------------
    res         = [nx_left nx_right;]
    resmin      = copy(res)
    MRefin      = 2.0
    RefineLevel = 7
    RefineCond  = 0.001
    nPoints     = 40
    BCout       = [0 0]

    #Check, if t_ar is valid (increasing in time)-----------------
    dt_diff = t_ar[2:end] .- t_ar[1:end-1]
    if any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    end

    #Create mesh, discretization and mapping----------------------
    if Ri[1] >= Ri[2]
        error("Please change the size of the system. Increase the total domain length.")
    elseif res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    x_left, x_right, dx1, dx2, x0 = create_grid!(Ri, res, MRefin, verbose)

    #Preprocess and initial condition-----------------------------
    L       = Ri[end]
    t       = 0.0
    it      = 0
    C_left  = Cl_i * ones(res[1], 1)
    C_right = Cr_i * C_left[end] * ones(res[2], 1) * inv(KD_ar[1])
    C0      = [copy(C_left); copy(C_right)]
    C       = copy(C0)
    x       = copy(x0)
    KD      = copy(KD_ar[1])

    Mass0   = calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)

    #Preallocate variables----------------------------------------
    dt      = 0.0
    L_g     = spzeros(length(x), length(x))
    Mass    = Float64[]
    R_g     = zeros(length(x), 1)
    MB_Error = Float64[]

    D_l, D_r, KD, T = update_t_dependent_param!(D0, Di, Ea_left, Ea_right, KD_ar, R, T_ar, t_ar, t, t_tot)

    #First check for correct setup--------------------------------
    if BCout[1] != 0 && (n == 3 || n == 2)
        error("The code is only valid for cylindrical/spherical geometry, where the left outer BC has Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T != T_ar[1]
        error("Initial temperature must be equal to the first value in the temperature array.")
    end

    #Time loop----------------------------------------------------
    while t < t_tot
        should_stop() && throw(SimulationCancelled())
        dt = find_dt(dx1, dx2, V_ip, D_l, D_r, CFL)
        t, dt, it = update_time!(t, dt, it, t_tot)
        report_progress(t * inv(t_tot))
        D_l, D_r, KD, T = update_t_dependent_param!(D0, Di, Ea_left, Ea_right, KD_ar, R, T_ar, t_ar, t, t_tot)
        # A no-op (Fl_regrid stays 0, x_left/x_right/res unchanged) when V_ip == 0, exactly as in
        # Diff_couple_Flux.jl; only actually advects/regrids once the interface can move.
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri, V_ip, dt, x_left, x_right, C_left, C_right, res)
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left, x_right, C_left, C_right, D_l, D_r, dt, n, res)
        L_g, R_g, ScF = set_inner_bc_flux!(L_g, R_g, KD, D_l, D_r, x_left, x_right, V_ip, rho, res)
        L_g, R_g = set_outer_bc!(BCout, L_g, R_g, Co_l[1], Co_r[end], ScF)
        C_left, C_right = solve_soe(L_g, R_g, res)
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin, RefineCond, RefineLevel, nPoints, RefineMethod, verbose)
        Massnow = calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)
        push!(Mass, Massnow)
        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error, ErrM)
        end
    end
    redirect_stdout(devnull) do
        calc_mass_err(Mass, Mass0)
    end
    return (; x_left, x_right, dx1, dx2, x0, res, Ri, C_left, C_right, C0, MB_Error, t_tot)
end
