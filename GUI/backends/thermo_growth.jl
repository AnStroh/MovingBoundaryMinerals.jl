using MovingBoundaryMinerals
using LinearAlgebra, SparseArrays, DelimitedFiles, Statistics

"""
    run_thermo_growth(; kwargs...)

GUI backend for the thermodynamically-constrained growth/resorption mode, adapted from
`examples/D2.jl`. `t_user`/`T_user` are the same user-defined temperature-time path D2.jl takes
(in seconds/Kelvin here rather than days/degC, since this is the programmatic entry point - the
GUI form layer converts from days/degC and validates the path before calling this); a 2-point
path reproduces `examples/D1.jl`'s linear-cooling behaviour, so no separate D1 wrapper is needed,
and a longer, non-monotonic path reproduces D2.jl's. Every hardcoded physics/numerics value that
a non-programmer would plausibly want to change is exposed as a keyword argument with the same
default. Advanced/rarely changed settings (crystallographic angles, resolution/refinement) are
kept fixed at the example's defaults. Never plots internally; always returns the raw arrays for
the caller to plot.

Unlike the original example, the phase-diagram CSV files are located relative to the installed
package directory (`pkgdir(MovingBoundaryMinerals)`), not the process's working directory, so this
works regardless of where the GUI server is launched from.

# Returns
`(; x_left, x_right, x0, C_left, C_right, C0, Tlin, XC_left, XC_right, Tstart, Tstop, KDlin, KD_sim, T_sim, Mass0, Mass, MB_Error, t_tot)`
"""
function run_thermo_growth(;
        Ri          = 0.0001,                          #Initial interface radius in [m]
        t_user      = [0.0, 30.0 * 24 * 60 * 60],       #Time nodes of the T-t path in [s]; must be sorted ascending and start at 0.0
        T_user      = [1400.0, 1350.0] .+ 273.0,        #Temperature nodes of the T-t path in [K]; same length as t_user, can go up and down (non-monotonic paths supported)
        P           = 1.0e6,      #Pressure in [Pa]
        D0          = 5.38e-9,    #Pre-exponential factor in [m^2/s]
        Ea          = 226000.0,   #Activation energy in [J/mol]
        CompInt     = 0.50,       #Composition of interest of the solid solution (Fe number)
        n           = 1,          #Geometry; 1: planar, 2: cylindrical, 3: spherical
        should_stop::Function = () -> false,   #Checked every iteration; throw(SimulationCancelled()) if true
        report_progress::Function = (_) -> nothing,   #Called every iteration with t/t_tot in [0, 1]
    )
    alpha       = 0.0
    beta        = 90.0
    gamma       = 90.0
    deltaV      = 7e-6
    R           = 8.314472
    verbose     = false
    RefineMethod = 1

    if length(t_user) != length(T_user)
        error("t_user and T_user must have the same length.")
    end
    if length(t_user) < 2
        error("The T-t path needs at least 2 points.")
    end
    if !issorted(t_user)
        error("t_user must be sorted in ascending order (time has to increase monotonically along the path).")
    end
    if t_user[1] != 0.0
        error("t_user[1] must be 0.0 (the T-t path has to start at the beginning of the simulation).")
    end
    t_tot = t_user[end]
    if t_tot <= 0.0
        error("The total time (t_user[end]) must be positive.")
    end

    pkgroot = pkgdir(MovingBoundaryMinerals)
    coeff       = readdlm(joinpath(pkgroot, "examples", "Examples_phase_diagram", "Coefficients_Reaction_lines.csv"))
    eq_values   = [coeff[1,1]  coeff[2,1]  coeff[3,1];
                   coeff[1,2]  coeff[2,2]  coeff[3,2]]
    rho_phases  = readdlm(joinpath(pkgroot, "examples", "Examples_phase_diagram", "density_phases copy.tab"))

    #T-t path-----------------------------------------------------------------
    tpath   = copy(t_user)
    Tpath   = copy(T_user)
    Tstart  = Tpath[1]
    Tstop   = Tpath[end]

    #Numerics---------------------------------------------------------------
    CFL         = 50.0
    res         = [50 50;]
    resmin      = copy(res)
    MRefin      = 5.0
    RefineLevel = 3
    RefineCond  = 0.01
    nPoints     = 20
    BCout       = [0 0]

    #Create data set--------------------------------------------------------
    coeff_up, coeff_do  = coeff_trans_line(eq_values)
    TMAX    = maximum(Tpath) + 1000.0
    TMIN    = minimum(Tpath) - 1000.0
    Tlin    = LinRange(TMAX, TMIN, 10000)
    XC_left, XC_right   = composition(coeff_up, coeff_do, Tlin)
    C_leftlin           = copy(XC_left)
    C_rightlin          = copy(XC_right)

    nd1         = Int(round(sqrt.(length(rho_phases[:,1]))))
    Xwm         = reshape(copy(rho_phases[:,1]), nd1, nd1)
    Twm         = reshape(copy(rho_phases[:,2]), nd1, nd1)
    rho_left    = ones(size(Xwm))
    rho_right   = ones(size(Xwm))

    R_left      = C_leftlin .* inv.((1.0 .- C_leftlin))
    R_right     = C_rightlin .* inv.((1.0 .- C_rightlin))
    KDlin       = R_left .* inv.(R_right)

    #Preprocess and initial condition---------------------------------------
    t   = 0.0
    it  = 0
    T   = copy(Tstart)
    C_leftB, C_rightB   = composition(coeff_up, coeff_do, T)
    Xc  = (1.0 - CompInt) * C_rightB + CompInt * C_leftB
    log10D_001    = log10(D0) - (Ea + (P - 1e5) * deltaV) / (2.303 * R * Tstart) + 3 * (C_leftB - 0.14)
    log10D_others = log10D_001 - log10(6.0)
    D_001 = 10^log10D_001
    D_010 = 10^log10D_others
    D_100 = 10^log10D_others
    D_l = D_001 * (cos(alpha))^2 + D_010 * (cos(beta))^2 + D_100 * (cos(gamma))^2
    D_r = exp(-7.92 - 26222 / Tstart)
    L   = (Ri[1]^n * (Xc - C_leftB) * inv(C_rightB - Xc))^(1 * inv(n)) + Ri[1]
    Ri  = [Ri L]

    if Ri[1] >= Ri[2]
        error("The interface radius is too large for this composition/geometry - please decrease it.")
    end
    if res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    x0_left, x0_right, dx1, dx2, x0 = create_grid!(Ri, res, MRefin, verbose)

    rho     = calculate_density(Xwm[:,1], Twm[1,:], rho_left, rho_right, C_leftB, C_rightB, T)
    x_left  = copy(x0_left)
    x_right = copy(x0_right)
    C_left  = C_leftB * ones(1, res[1])
    C_right = C_rightB * ones(1, res[2])
    C0      = [copy(C_left) copy(C_right)]
    dt      = minimum([dx1, dx2])^2 * inv(maximum([D_l, D_r]))

    Mass0   = calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)

    L_g     = spzeros(length(x0), length(x0))
    Mass    = Float64[]
    KD_sim  = Float64[]
    T_sim   = Float64[]
    R_g     = zeros(length(x0), 1)
    Residual = Float64[]
    MB_Error = Float64[]

    while t < t_tot
        should_stop() && throw(SimulationCancelled())
        t, dt, it = update_time!(t, dt, it, t_tot)
        report_progress(t * inv(t_tot))
        if t <= t_tot
            T = linear_interpolation_1D(tpath, Tpath, t)
        end
        C_left[end], C_right[1] = composition(coeff_up, coeff_do, T)
        dC  = C_right[1] - C_left[end]
        rho = calculate_density(Xwm[:,1], Twm[1,:], rho_left, rho_right, C_leftB, C_rightB, T)

        log10D_001    = log10(D0) - (Ea + (P - 1e5) * deltaV) / (2.303 * R * T) + 3 * (mean(C_left) - 0.14)
        log10D_others = log10D_001 - log10(6.0)
        D_001 = 10^log10D_001
        D_010 = 10^log10D_others
        D_100 = 10^log10D_others
        D_l = D_001 * (cos(alpha))^2 + D_010 * (cos(beta))^2 + D_100 * (cos(gamma))^2
        D_r = exp(-7.92 - 26222 / T)

        JL   = -D_l * rho[1] * (C_left[end] - C_left[end-1]) * inv(dx1)
        JR   = -D_r * rho[2] * (C_right[2] - C_right[1]) * inv(dx2)
        V_ip = (JR - JL) * inv(dC)
        push!(Residual, (V_ip * dC - (JR - JL))^2)

        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri, V_ip, dt, x_left, x_right, vec(C_left), vec(C_right), res)
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left, x_right, C_left, C_right, D_l, D_r, dt, n, res)
        L_g, R_g, ScF = set_inner_bc_stefan!(L_g, R_g, C_left, C_right, res)
        L_g, R_g = set_outer_bc!(BCout, L_g, R_g, Co_l[1], Co_r[end], ScF)
        C_left, C_right = solve_soe(L_g, R_g, res)

        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin, RefineCond, RefineLevel, nPoints, RefineMethod, verbose)

        Massnow = calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)
        R_left  = C_left[end] .* inv(1.0 - C_left[end])
        R_right = C_right[1] .* inv(1.0 - C_right[1])
        KDt     = R_left .* inv(R_right)
        push!(Mass, Massnow)
        push!(KD_sim, KDt)
        push!(T_sim, copy(T))

        dtV = minimum([dx1, dx2])^1 * inv(abs(V_ip))
        dtD = minimum([dx1, dx2])^2 * inv(maximum([D_l, D_r]))
        dt  = minimum([dtD, dtV]) * CFL

        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error, ErrM)
        end
    end
    redirect_stdout(devnull) do
        calc_mass_err(Mass, Mass0)
    end
    return (; x_left, x_right, x0, C_left = vec(C_left), C_right = vec(C_right), C0 = vec(C0),
             Tlin, XC_left, XC_right, Tstart, Tstop, KDlin, KD_sim, T_sim, Mass0, Mass, MB_Error, t_tot)
end
