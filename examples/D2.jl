using MovingBoundaryMinerals
using Plots, LinearAlgebra, DelimitedFiles, SparseArrays, LaTeXStrings, Statistics, Dates
#Main function--------------------------------------------------------------
#This example is a variant of D1.jl that shows how to prescribe an arbitrary,
#user-defined temperature-time (T-t) path instead of a linearly decreasing
#one. The path does not need to be monotonic: temperature can go up and down
#between nodes. The default path below falls, plateaus, then rises gently -
#genuinely non-monotonic, but with modest slope changes at each corner. At
#every time step the actual temperature is obtained by linear interpolation
#between the user-supplied (t_user, T_user) nodes with `linear_interpolation_1D`,
#which only requires `t_user` to be sorted -> `T_user` itself can move freely
#up and down.
#Fixed 2026-08-28: the previous default path (1400->1390->1390->1380->1385->
#1385->1375->1365->1365->1370->1350°C over 30 days, at Ri=0.0001, matching
#D1.jl's Tstart/Tstop/Ri) has sharp corners - large slope changes packed into
#3-day intervals - that drive the interface velocity into a stiff regime the
#regridder cannot keep up with. Under the corrected cosd() physics (see the
#cos()/cosd() fix above), a full run eventually hits "Cannot proceed (Newton
#Failure)" after several hours. This default instead uses a larger seed
#(Ri=0.0005, which also makes every run faster: bigger Ri means bigger dx,
#hence bigger dt at the same CFL) and much gentler corners, which complete the
#full 30 days in well under a minute. See CHANGELOG.md for details.
function D2(;RefineMethod = 1, plot_sim = false, verbose= false, animate_sim = false, frame_every = 150,
             t_user = [0.0, 10.0, 20.0, 30.0] .* (60*60*24),                #Time nodes of the T-t path in [s]; USER INPUT - must be sorted in ascending order, does not need to be regularly spaced. Default: 30 days total, falls then plateaus then rises gently
             T_user = [1400.0, 1385.0, 1385.0, 1390.0] .+ 273.0)      #Temperature nodes of the T-t path in [K]; USER INPUT - can increase or decrease between nodes (non-monotonic paths are supported). Default: falls 1400->1385°C (day 0-10), plateaus (day 10-20), rises 1385->1390°C (day 20-30)
    @info "D2.jl: as of 2026-08-26, the anisotropic diffusion weighting (D_l, below) uses cosd() " *
          "(degrees) instead of the previous cos() (radians) bug - results differ by ~6.7% from " *
          "versions of this package before that fix. See CHANGELOG.md and D1.jl for details."
    #If you find a [] with two entries this belong to the respective side of
    #the diffusion couple ([left right])
    #Check user-defined T-t path----------------------------------------------
    if length(t_user) != length(T_user)
        error("t_user and T_user must have the same length. Please check your inputs.")
    end
    if !issorted(t_user)
        error("t_user must be sorted in ascending order (time has to increase monotonically along the path).")
    end
    if t_user[1] != 0.0
        error("t_user[1] must be 0.0 (the T-t path has to start at the beginning of the simulation).")
    end
    #Physics-----------------------------------------------------------------
    D0          = 5.38*1e-9                                                                             #Pre-exponential factor in [m^2/s]
    Ea          = 226000                                                                                #Activation energy for the left side in [J/mol]
    alpha       = 0.0                                                                                   #Angle between [001] & profile
    beta        = 90.0                                                                                  #Angle between [010] & profile
    gamma       = 90.0                                                                                  #Angle between [100] & profile
    deltaV      = 7*10^-6                                                                               #Volume change in [m^3/mol]
    Ri          = 0.0005                                                                                #Position of the interface -> initial radius of the left phase in [m]. 5x D1.jl's 0.0001: see the note above the function signature.
    P           = 10^6                                                                                  #Pressure in [Pa]
    R           = 8.314472                                                                              #Universal gas constant in [J/(mol*K)]
    Myr2Sec     = 60*60*24*365.25*1e6                                                                   #Conversion factor from Myr to s
    n           = 1                                                                                     #Geometry; 1: planar, 2: cylindrical, 3: spherical
    CompInt     = 0.50                                                                                  #Composition of interest of the solid solution (Fe number)
    coeff       = readdlm("examples/Examples_phase_diagram/Coefficients_Reaction_lines.csv")            #Reads the coefficients for the linear least squares
    eq_values   = [coeff[1,1]  coeff[2,1]  coeff[3,1];	                                                #Coefficients for composition calculation of component B (stable at higher T) X2 = a2 + b2*T + c2*T²
                   coeff[1,2]  coeff[2,2]  coeff[3,2]]                                                  #Coefficients for composition calculation of component A (stable at lower T) X1 = a1 + b1*T + c1*T²
    rho_phases  = readdlm("./examples/Examples_phase_diagram/density_phases copy.tab")                  #Reads the density values for the phases
    #T-t path-----------------------------------------------------------------
    tpath               = copy(t_user)                                                                  #Time path in [s], taken directly from the user-defined array
    Tpath               = copy(T_user)                                                                  #Temperature path in [K], taken directly from the user-defined array (can go up and down)
    t_tot               = tpath[end]                                                                    #Total time [s] -> last node of the user-defined T-t path
    Tstart              = Tpath[1]                                                                      #Starting temperature in [K]
    Tstop               = Tpath[end]                                                                    #End temperature in [K]
    #Numerics---------------------------------------------------------------
    CFL                 = 50.0                                                                          #CFL condition
    res                 = [50 50;]                                                                      #Number of nodes
    resmin              = copy(res)                                                                     #Minimum number of nodes
    MRefin              = 5.0                                                                           #Refinement factor
    RefineLevel         = 3                                                                            #Refinement level; how many times should the grid be refined
    RefineCond          = 0.01                                                                           #Refinement condition; refine until last dx on the left side <= RefineCond * Ri[1]
    nPoints             = 20                                                                             #Number of points for initial grid (h-refinement); It can be a number, if both sides should have the same initial number of points or it can be an array ([left right]), if the resolution for the initial number of points should differ.
    BCout               = [0 0]                                                                         #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann;
                                                                                                        #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Create data set--------------------------------------------------------
    #Create arrays X(T) using linear least squares
    coeff_up, coeff_do  = coeff_trans_line(eq_values)                                                   #Extract coefficients
    TMAX                = maximum(Tpath) + 1000.0                                                       #Max temperature for T-array -> based on the highest point of the user-defined T-t path
    TMIN                = minimum(Tpath) - 1000.0                                                       #Min temperature for T-array -> based on the lowest point of the user-defined T-t path
    Tlin                = LinRange(TMAX,TMIN,10000)                                                     #Temperature values
    XC_left, XC_right   = composition(coeff_up,coeff_do,Tlin)                                           #e.g. liquidus/solidus/solvus
    XC_left_Cel, XC_right_Cel = composition(coeff_up,coeff_do,Tlin .- 273.0)                            #e.g. liquidus/solidus/solvus
    C_leftlin           = copy(XC_left)                                                                 #Store the composition of A for later calculations
    C_rightlin          = copy(XC_right)                                                                #Store the composition of B for later calculations
    #Create density plots---------------------------------------------------
    nd1                 = Int(round(sqrt.(length(rho_phases[:,1]))))
    Xwm                 = copy(rho_phases[:,1])
    Xwm                 = reshape(Xwm,nd1,nd1)                                                          #X(C1)
    Twm                 = copy(rho_phases[:,2])
    Twm                 = reshape(Twm,nd1,nd1)                                                          #Temperature in [K]
    rho_left            = ones(size(Xwm))
    rho_right           = ones(size(Xwm))
    #Create other arrays----------------------------------------------------
    R_left              = C_leftlin .* inv.((1.0 .- C_leftlin))                                         #Composition rate in phase A
    R_right             = C_rightlin .* inv.((1.0 .- C_rightlin))                                       #Composition rate in phase B
    KDlin               = R_left .* inv.(R_right)                                                       #KD coefficient
    #Preprocess and initial condition---------------------------------------
    t                   = 0.0                                                                           #Time
    it                  = 0                                                                             #Iterations
    T                   = copy(Tstart)                                                                  #Initial temperature
    C_leftB, C_rightB   = composition(coeff_up,coeff_do,T)                                              #Initial composition of the phases
    Xc                  = (1.0 - CompInt) * C_rightB + CompInt * C_leftB                                #Actual total composition (Xc intersection with black dashed line, e.g. Xc = C(Ol)*V(Ol)+C(melt)*V(melt)
    log10D_001          = log10(D0) - (Ea + (P-10^5)* deltaV)/(2.303*R*Tstart) + 3*(C_leftB - 0.14)     #log10(diffusion coefficient Olivine) following Dohmen & Chakraborty (2007)
    log10D_others       = log10D_001 - log10(6.0)                                                       #log10(diffusion coefficient Olivine) following Dohmen & Chakraborty (2007)
    D_001               = 10^log10D_001                                                                 #Diffusion coefficient in direction 001
    D_010               = 10^log10D_others                                                              #Diffusion coefficient in direction 010
    D_100               = 10^log10D_others                                                              #Diffusion coefficient in direction 100
    #Fixed 2026-08-26: alpha/beta/gamma are angles in degrees (see their definitions above), but this
    #previously called cos() (radians) instead of cosd() (degrees) - cosd(90) = 0 as intended, whereas
    #cos(90) ≈ -0.448, so the two axes meant to be perpendicular to the profile (and so contribute
    #nothing) were picking up a spurious cos(90)^2 ≈ 0.2 weight each. For the default alpha=0, beta=90,
    #gamma=90 this overestimated D_l by about 6.7%. See D1.jl for the same fix and CHANGELOG.md.
    D_l                 = D_001*(cosd(alpha))^2 + D_010*(cosd(beta))^2 + D_100*(cosd(gamma))^2             #Effective diffusion coefficient Olivine after Crank (1975), p. 7
    D_r                 = exp(-7.92-26222/(Tstart))                                                     #Diffusivity melt: Approximated after Zhang & Cherniak (2010) p. 332 EQ. 19
    L                   = (Ri[1] ^ n * (Xc - C_leftB) * inv((C_rightB - Xc)) ^ (1 * inv(n))) + Ri[1]    #Total length of the modelling domain
    Ri                  = [Ri L]                                                                        #Radii of the 2 phases
    #Create mesh, discretization and mapping--------------------------------
    if Ri[1] >= Ri[2]
        error("Please change the size of the system. Increase Ri[2].")
    end
    if res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    if RefineMethod == 1
        x0_left, x0_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    elseif RefineMethod == 2
        x0_left, x0_right, dx1, dx2, x0 = h_refinement1(Ri,RefineLevel,nPoints)
        res = [length(x0_left) length(x0_right)]
        resmin = copy(res)
    elseif RefineMethod == 3
        x0_left, x0_right, dx1, dx2, x0 = h_refinement2(Ri,RefineCond,nPoints)
        res = [length(x0_left) length(x0_right)]
        resmin = copy(res)
    else
        error("RefineMethod not valid. Please choose 1, 2 or 3.")
    end
    #Calculate densities----------------------------------------------------
    rho                 = calculate_density(Xwm[:,1],Twm[1,:],rho_left,rho_right,C_leftB,C_rightB,T)    #Initial normalized densities in [-]
    #More initial conditions------------------------------------------------
    x_left              = copy(x0_left)                                                                 #Initial x_left
    x_right             = copy(x0_right)                                                                #Initial x_right
    C_left              = C_leftB * ones(1,res[1])                                                      #Composition of component B in phase A
    C_right             = C_rightB * ones(1,res[2])                                                     #Composition of component B in phase B
    C0                  = [copy(C_left) copy(C_right)]                                                  #Store initial composition
    dt                  = minimum([dx1,dx2]) ^ 2 .* inv((maximum([D_l,D_r])))                           #Initial dt
    #Total mass-------------------------------------------------------------
    Mass0               = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)                            #Initial mass
    Mass01              = (trapezoidal_integration(x_left.^n*rho[1],C_left)+ trapezoidal_integration(x_right.^n*rho[2],C_right))/Ri[2]    #Initial composition (plot)
    #Preallocate variables--------------------------------------------------
    L_g                 = spzeros(length(x0),length(x0))                                                #Global LHS matrix
    Mass                = Float64[]                                                                     #Array to store the mass of the system
    Mass2               = Float64[]                                                                     #Array to store the mass of the system (plot)
    KD_sim              = Float64[]                                                                     #Array to store the distribution coefficient
    T_sim               = Float64[]                                                                     #Array to store the temperature
    R_left_sim          = Float64[]                                                                     #Array to store the composition ratio left side (interface)
    R_right_sim         = Float64[]                                                                     #Array to store the composition ratio right side (interface)
    R_g                 = zeros(length(x0),1)                                                           #Global RHS vector
    #Checks-----------------------------------------------------------------
    C_left_check        = [C_left[end]]                                                                 #Check composition left side
    C_right_check       = [C_right[1]]                                                                  #Check composition right side
    T_check             = [Tstart]                                                                      #Check temperature
    Residual            = Float64[]                                                                     #Residual of the velocity
    MB_Error            = Float64[]                                                                     #Mass error
    anim                = animate_sim ? Animation() : nothing                                            #Collects one frame every `frame_every` iterations if animate_sim; movie of the composition profile/phase-diagram evolution
    #-----------------------------------------------------------------------
    #Solving the moving boundary problem------------------------------------
    while t < t_tot
        #Update time--------------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        if t <= t_tot
            T  = linear_interpolation_1D(tpath,Tpath,t)                                                #Interpolate T at t from the user-defined (tpath,Tpath) nodes; supports non-monotonic paths
        end
        #Calculate Equilibrium compositions at actual T---------------------
        C_left[end], C_right[1] = composition(coeff_up,coeff_do,T)
        dC  = C_right[1] - C_left[end]                                                                  #Composition difference
        rho = calculate_density(Xwm[:,1],Twm[1,:],rho_left,rho_right,C_leftB,C_rightB,T)
        #Calculate diffusivities--------------------------------------------
        log10D_001    = log10(D0) - (Ea + (P-10^5)* deltaV)/(2.303*R*T) + 3*(mean(C_left) - 0.14)       #log10(diffusion coefficient Olivine) following Dohmen & Chakraborty (2007)
        log10D_others = log10D_001 - log10(6.0)                                                         #log10(diffusion coefficient Olivine) following Dohmen & Chakraborty (2007)
        D_001         = 10^log10D_001                                                                   #Diffusion coefficient in direction 001
        D_010         = 10^log10D_others                                                                #Diffusion coefficient in direction 010
        D_100         = 10^log10D_others                                                                #Diffusion coefficient in direction 100
        D_l           = D_001*(cosd(alpha))^2 + D_010*(cosd(beta))^2 + D_100*(cosd(gamma))^2               #Effective diffusion coefficient after Crank (1975), p. 7
        D_r           = exp(-7.92-26222/(T))                                                            #Diffusivity melt: Approximated after Zhang & Cherniak (2010) p. 332 EQ. 19
        #Stefan condition -> Composition difference-------------------------
        JL   = - D_l * rho[1] * (C_left[end] - C_left[end-1]) * inv(dx1)                                #Flux of the left side to the right side
        JR   = - D_r * rho[2] * (C_right[2]  - C_right[1])    * inv(dx2)                                #Flux of the right side to the left side
        V_ip = (JR - JL) * inv(dC)                                                                      #Velocity in x direction
        #Debugging velocity-------------------------------------------------
        push!(Residual,(V_ip*dC-(JR-JL))^2)
        #Advect interface & regrid------------------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,vec(C_left),vec(C_right),res)
        #FEM SOLVER---------------------------------------------------------
        #Construct global matrices------------------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,res)
        #Set inner boundary conditions--------------------------------------
        L_g, R_g, ScF = set_inner_bc_stefan!(L_g,R_g,C_left,C_right,res)
        #Set outer boundary conditions and scale matrices-------------------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #Solve system-------------------------------------------------------
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Debugging composition-----------------------------------------------
        push!(C_left_check,C_left[end])
        push!(C_right_check,C_right[1])
        push!(T_check, T)
        #Regrid-------------------------------------------------------------
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,RefineCond,RefineLevel,nPoints,RefineMethod,verbose)
        #Post-Preprocessing-------------------------------------------------
        for iit in enumerate(1)
            Massnow     = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            Massnow2    = (trapezoidal_integration(x_left.^n*rho[1],C_left)+ trapezoidal_integration(x_right.^n*rho[2],C_right))/Ri[2]
            R_left      = C_left[end] .* inv((1.0 - C_left[end]))
            R_right     = C_right[1] .* inv((1.0 - C_right[1]))
            KDt         = R_left .* inv(R_right)
            T_pl        = copy(T)
            push!(Mass, Massnow)                                                                        #Stores the mass of the system
            push!(Mass2, Massnow2)                                                                      #Stores the mass of the system (plot)
            push!(KD_sim, KDt)                                                                          #Stores the distribution coefficient
            push!(T_sim, T_pl)                                                                          #Stores the temperature
            push!(R_left_sim, R_left)                                                                   #Composition ratio left side (interface)
            push!(R_right_sim, R_right)                                                                 #Composition ratio right side (interface)
        end
        #Find the important dt----------------------------------------------
        dtV   = minimum([dx1,dx2]) ^1 .* inv((abs(V_ip)))
        dtD   = minimum([dx1,dx2]) ^2 .* inv((maximum([D_l,D_r])))
        dt    = minimum([dtD,dtV]) * CFL
        #Plotting-----------------------------------------------------------
        if (plot_sim && it % 200 == 0) || (animate_sim && it % frame_every == 0)
            #Plotting-------------------------------------------------------
            fs = 12.0
            maxC = maximum([maximum(C_left),maximum(C_right)])
            Tp_min = (minimum(Tpath) - 273.0) * 0.95                                                    #Based on the min/max of the user-defined T-t path (not just start/stop, since the path can go up and down)
            Tp_max = (maximum(Tpath) - 273.0) * 1.05
            #Dense tick marks, but only every 2nd one labelled (and never the very first,
            #which would sit on top of the y-axis "0.0" label) to avoid crowding on this
            #axis' range, which varies with the user's T-t path.
            Tp_ticks = range(Tp_min, Tp_max, length=11)[2:2:end]                                            #Every 2nd of 11 evenly-spaced points, skipping the very first (which would sit on the y-axis "0.0" label) - only labelled ticks are drawn, no unlabelled ones in between
            Tp_labels = string.(round.(Int, Tp_ticks))
            first_val, last_val = values_between_known_indices!(Tlin.-273.0,KDlin,maximum(Tpath)-273.0,minimum(Tpath)-273.0)

            #Composition profile
            p1 = plot(x_left*1e3,C_left, lw=2, label=L"\mathrm{Left\ side}")
            p1 = plot!(x_right*1e3,C_right, lw=2, label=L"\mathrm{Right\ side}")
            p1 = plot!(x0*1e3,vec(C0), label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
                  ylabel = L"X\mathrm{_{Fe}}", lw=1.5, grid=:on, legend = :right)
            p1 = plot!([x_left[end]; x_left[end];]*1e3, [0; 1*(maxC + 0.01)], color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}",ylim=[C0[1]-0.05; 1*(maxC + 0.01)])
            #(a)'s y-position is a fraction of panel 1's own y-range (set two lines up) rather than a fixed
            #data value: unlike D1.jl's monotonic cooling (where maxC only grows), D2's non-monotonic T-t
            #path can transiently pull maxC below a hardcoded 0.75, pushing the label above the panel.
            ylo_a, yhi_a = C0[1] - 0.05, 1*(maxC + 0.01)
            #(a)'s x-position is a fraction of the panel's own plotted x-range (0 to the domain's
            #right edge) rather than a fixed 0.015 mm: that collided with the y-axis tick labels
            #whenever the domain is small enough (e.g. a larger Ri) that 0.015 mm sits right on top
            #of them instead of comfortably inside the panel.
            p1 = annotate!(0.075*x_right[end]*1e3, ylo_a + 0.95*(yhi_a - ylo_a), L"\mathrm{(a)}")
            #Phase diagram
            p2 = plot(Tlin .- 273.0,XC_left, lw=2, label=L"\mathrm{Left\ side}")
            p2 = plot!(Tlin .- 273.0,XC_right, lw=2, label=L"\mathrm{Right\ side}")
            p2 = scatter!([T-273.0],[C_left[end]],marker=:circle, markersize=2, markercolor=:grey68,
                          markerstrokecolor=:grey68,label = "",ylabel = L"X\mathrm{_{Fe}}",xlabel=L"T\ \mathrm{[°C]}")
            p2 = scatter!([T-273.0],[C_right[1]],marker=:circle, markersize=2, markercolor=:grey68,
                          markerstrokecolor=:grey68,label = "")
            p2 = scatter!([Tstart-273.0],[C0[51]],marker=:circle, markersize=2, markercolor=:black,
                          markerstrokecolor=:black,label = "",ylabel = L"X\mathrm{_{Fe}}",xlabel=L"T\ \mathrm{[°C]}")
            p2 = scatter!([Tstart-273.0],[C0[50]],marker=:circle, markersize=2, markercolor=:black,
                          markerstrokecolor=:black,label = "")
            p2 = plot!([Tstart-273.0; Tstart-273.0],[0; maximum(C0[end])],lw=1.5,color=:black,linestyle=:dash,label=L"T(t=0.0)")
            p2 = plot!([T-273.0; T-273.0],[0; maximum([C_left[end],C_right[1]])],lw=1.5,color=:grey68,linestyle=:dashdot,label=L"T(t)")
            p2 = plot!([T-273.0; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:royalblue,linestyle =:dot)
            p2 = plot!([T-273.0; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), xticks=(Tp_ticks, Tp_labels), ylims=(0, 1),color=:crimson,linestyle =:dot)
            #(b)'s x-position is a fraction of this panel's own (Tp_min, Tp_max) range rather than a
            #fixed 1300 degC: for a narrower/different T-t path than the default, 1300 can sit right
            #on the y-axis tick labels instead of comfortably inside the panel.
            p2 = annotate!(Tp_min + 0.1*(Tp_max - Tp_min), 0.95, L"\mathrm{(b)}")
            #Evolution of KD(T)
            p3 = plot(Tlin .- 273.0, KDlin, lw=1.5, label=L"\mathrm{Thermodynamic\ data}", color=:black)
            p3 = scatter!([T_sim[end]-273.0],[KD_sim[end]],marker=:circle, markersize=3.0, markercolor=:black,markerstrokecolor=:black,
                          xlabel = L"T\ \mathrm{[°C]}", ylabel = L"K_{D}", lw=1.5,legend = :bottomleft,
                          grid=:on, label=L"\mathrm{Model}",xlims=(Tp_min, Tp_max), xticks=(Tp_ticks, Tp_labels), ylims=(first_val-0.01,last_val+0.01))
            #ln(KD) vs 1/T
            p4 = plot(10000.0 ./ (T_sim),log.(KD_sim),xlabel = L"10,000/T\ \mathrm{[K^{-1}]}", ylabel = L"ln(K_{D})", lw=1.5,
                    grid=:on, label="", color=:black, ticks=:auto, xrotation=0)
            #Figure 1
            p = plot(p1,p2,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                    legend_foreground_color = :transparent)
            plot_sim && it % 200 == 0 && display(p)
            animate_sim && it % frame_every == 0 && frame(anim, p)
        end
        # Suppress output of calc_mass_err
        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error,ErrM)
        end
    end
    #Post-process-----------------------------------------------------------
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, vec(C_left), vec(C_right), vec(C0),maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, Tpath, KDlin, KD_sim,T_sim, Mass0, Mass, Mass01, Mass2, C_left_check, C_right_check, T_check, Residual, MB_Error, anim
end
#Run calculation------------------------------------------------------------
# Refinement method: 1 = m-refinement, 2 = h-refinement based on number of refinement levels, 3 = h-refinement based on refinement condition (first/last dx on the left side)
run_and_plot = true
run_and_plot == false ? printstyled("You have disabled the simulation, change the variable run_and_plot == true", bold=true) : nothing
if run_and_plot
    plot_sim    = false
    plot_end    = true
    verbose     = false
    save_file   = false
    animate_sim = true                                                                                          #Set to true to save a movie of the composition profile/phase-diagram evolution (see save_path/save_name below)
    frame_every = 150                                                                                            #Capture a movie frame every this many iterations; smaller = smoother but slower to render and a bigger file
    #User-defined temperature-time path -> can be edited freely, temperature can go up and down between nodes.
    #Falls, plateaus, then rises gently over 30 days - see the note above the D2() function signature for why.
    t_user = [0.0, 10.0, 20.0, 30.0] .* (60*60*24)                                                               #Time nodes in [s]
    T_user = [1400.0, 1385.0, 1385.0, 1390.0] .+ 273.0                                                           #Temperature nodes in [K]: falls 1400->1385°C, plateaus, rises 1385->1390°C
    x_left, x_right, x0, C_left, C_right, C0, maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, Tpath, KDlin, KD_sim,T_sim, Mass0, Mass, Mass01, Mass2, C_left_check, C_right_check, T_check,Residual, MB_Error, anim = D2(RefineMethod = 1, plot_sim=plot_sim, verbose=verbose, animate_sim=animate_sim, frame_every=frame_every, t_user=t_user, T_user=T_user)
    save_path = "figures"
    if animate_sim
        #Movie of the composition profile/phase-diagram evolution over the whole T-t path
        !isdir(save_path) && mkdir(save_path)
        current_time = Dates.format(Dates.now(), "dd_mm_yy_HHMMSS")
        gif(anim, joinpath(save_path, "D2_$(current_time).gif"), fps = 10)
    end
    if plot_end
        #Title: Thermodynamical constrained Stefan condition with a user-defined, non-monotonic T-t path
        #Plotting in K-----------------------------------------------------------
        Tstart = Tstart - 273.0
        Tstop  = Tstop  - 273.0
        Tp_min = (minimum(Tpath) - 273.0) * 0.95                                                        #Based on the min/max of the user-defined T-t path (not just start/stop, since the path can go up and down)
        Tp_max = (maximum(Tpath) - 273.0) * 1.05
        #Dense tick marks, but only every 2nd one labelled (and never the very first,
        #which would sit on top of the y-axis "0.0" label) to avoid crowding on this
        #axis' range, which varies with the user's T-t path.
        Tp_ticks = range(Tp_min, Tp_max, length=11)[2:2:end]                                            #Every 2nd of 11 evenly-spaced points, skipping the very first (which would sit on the y-axis "0.0" label) - only labelled ticks are drawn, no unlabelled ones in between
        Tp_labels = string.(round.(Int, Tp_ticks))
        first_val, last_val = values_between_known_indices!(Tlin.-273.0,KDlin,maximum(Tpath)-273.0,minimum(Tpath)-273.0)
        fs = 12.0
        #Composition profile
        p1 = plot(x_left*1e3,C_left, lw=2, label=L"\mathrm{Left\ side}")
        p1 = plot!(x_right*1e3,C_right, lw=2, label=L"\mathrm{Right\ side}")
        p1 = plot!(x0*1e3,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
              ylabel = L"X\mathrm{_{Fe}}", lw=1.5, grid=:on, legend = :right)
        p1 = plot!([x_left[end]; x_left[end];]*1e3, [0; 1*(maxC + 0.01)], color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}",ylim=[C0[1]-0.05; 1*(maxC + 0.01)])
        #(a)'s y-position is a fraction of panel 1's own y-range (set two lines up) rather than a fixed
        #data value: unlike D1.jl's monotonic cooling (where maxC only grows), D2's non-monotonic T-t
        #path can transiently pull maxC below a hardcoded 0.75, pushing the label above the panel.
        ylo_a, yhi_a = C0[1] - 0.05, 1*(maxC + 0.01)
        p1 = annotate!(0.075*x_right[end]*1e3, ylo_a + 0.95*(yhi_a - ylo_a), L"\mathrm{(a)}")
        #Phase diagram
        p2 = plot(Tlin .- 273.0,XC_left, lw=2, label=L"\mathrm{Left\ side}")
        p2 = plot!(Tlin .- 273.0,XC_right, lw=2, label=L"\mathrm{Right\ side}")
        p2 = scatter!([T-273.0],[C_left[end]],marker=:circle, markersize=2, markercolor=:grey68,
                      markerstrokecolor=:grey68,label = "",ylabel = L"X\mathrm{_{Fe}}",xlabel=L"T\ \mathrm{[°C]}")
        p2 = scatter!([T-273.0],[C_right[1]],marker=:circle, markersize=2, markercolor=:grey68,
                      markerstrokecolor=:grey68,label = "")
        p2 = scatter!([Tstart],[C0[length(x_left)+1]],marker=:circle, markersize=2, markercolor=:black,
                      markerstrokecolor=:black,label = "",ylabel = L"X\mathrm{_{Fe}}",xlabel=L"T\ \mathrm{[°C]}")
        p2 = scatter!([Tstart],[C0[length(x_left)]],marker=:circle, markersize=2, markercolor=:black,
                      markerstrokecolor=:black,label = "")
        p2 = plot!([Tstart; Tstart],[0; maximum(C0[end])],lw=1.5,color=:black,linestyle=:dash,label=L"T(t=0.0)")
        p2 = plot!([T-273.0; T-273.0],[0; maximum([C_left[end],C_right[1]])],lw=1.5,color=:grey68,linestyle=:dashdot,label=L"T(t\mathrm{_{tot})}")
        p2 = plot!([T-273.0; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:royalblue,linestyle =:dot)
        p2 = plot!([T-273.0; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), xticks=(Tp_ticks, Tp_labels), ylims=(0, 1),color=:crimson,linestyle =:dot)
        p2 = annotate!(Tp_min + 0.1*(Tp_max - Tp_min), 0.95, L"\mathrm{(b)}")
        #Evolution of KD(T)
        p3 = plot(Tlin .- 273.0, KDlin, lw=1.5, label=L"\mathrm{Thermodynamic\ data}", color=:black)
        p3 = scatter!([T_sim[end]-273.0],[KD_sim[end]],marker=:circle, markersize=3.0, markercolor=:black,markerstrokecolor=:black,
                      xlabel = L"T\ \mathrm{[°C]}", ylabel = L"K_{D}", lw=1.5,legend = :bottomleft,
                      grid=:on, label=L"\mathrm{Model}",xlims=(Tp_min, Tp_max), xticks=(Tp_ticks, Tp_labels), ylims=(first_val-0.01,last_val+0.01))
        #T-t path -> shows the user-defined, non-monotonic T-t path used for the linear interpolation in the simulation
        p5 = plot(t_user ./ (60*60*24), Tpath .- 273.0, lw=1.5, color=:black, label=L"\mathrm{Interpolated\ path}",
                  xlabel=L"t\ \mathrm{[days]}", ylabel=L"T\ \mathrm{[°C]}", grid=:on)
        p5 = scatter!(t_user ./ (60*60*24), Tpath .- 273.0, marker=:circle, markersize=4, markercolor=:black,
                      markerstrokecolor=:black, label=L"\mathrm{User\ nodes}")
        #ln(KD) vs 1/T
        p4 = plot(10000.0 ./ (T_sim),log.(KD_sim),xlabel = L"10,000/T\ \mathrm{[K^{-1}]}", ylabel = L"ln(K_{D})", lw=1.5,
                grid=:on, label="", color=:black, ticks=:auto, xrotation=0)
        #Figure 1
        plot(p1,p2,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                legend_foreground_color = :transparent)
        display(current())
        save_path = "figures"
        save_name = "D2"
        save_figure(save_name,save_path,save_file)
        #Figure 2: the user-defined, non-monotonic T-t path
        plot(p5,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                legend_foreground_color = :transparent)
        display(current())
        save_name = "D2"
        save_figure(save_name,save_path,save_file)
    end
end
