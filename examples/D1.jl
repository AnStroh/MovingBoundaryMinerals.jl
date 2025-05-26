using MOBILE
using Plots, LinearAlgebra, DelimitedFiles, SparseArrays, LaTeXStrings, Statistics
#Main function--------------------------------------------------------------
function D1(; plot_sim = false, verbose = false)
    #If you find a [] with two entries this belong to the respective side of
    #the diffusion couple ([left right])
    #Physics-----------------------------------------------------------------
    D0          = 5.38*1e-9                                                                             #Pre-exponential factor in [m^2/s]
    Ea          = 226000                                                                                #Activation energy for the left side in [J/mol]
    alpha       = 90.0                                                                                  #Ideal angle between [001] and [010] (b-c-plane)
    beta        = 90.0                                                                                  #Ideal angle between [001] and [100] (a-c-plane)
    gamma       = 90.0                                                                                  #Ideal angle between [010] and [100] (a-b-plane)
    deltaV      = 7*10^-6                                                                               #Volume change in [m^3/mol]
    Ri          = 0.0001                                                                                #Position of the interface -> initial radius of the left phase in [m]
    Tstart      = 1400.0 + 273.0                                                                        #Starting temperature in [K]
    Tstop       = 1350.0 + 273.0                                                                        #End temperature in [K]
    P           = 10^6                                                                                  #Pressure in [Pa]
    R           = 8.314472                                                                              #Universal gas constant in [J/(mol*K)]
    Myr2Sec     = 60*60*24*365.25*1e6                                                                   #Conversion factor from Myr to s
    t_tot       = 60*60*24*30                                                                           #Total time [s]
    n           = 1                                                                                     #Geometry; 1: planar, 2: cylindrical, 3: spherical
    CompInt     = 0.50                                                                                  #Composition of interest of the solid solution (Fe number)
    coeff       = readdlm("examples/Examples_phase_diagram/Coefficients_Reaction_lines.csv")            #Reads the coefficients for the linear least squares
    eq_values   = [coeff[1,1]  coeff[2,1]  coeff[3,1];	                                                #Coefficients for composition calculation of component B (stable at higher T) X2 = a2 + b2*T + c2*T²
                   coeff[1,2]  coeff[2,2]  coeff[3,2]]                                                  #Coefficients for composition calculation of component A (stable at lower T) X1 = a1 + b1*T + c1*T²
    rho_phases  = readdlm("./examples/Examples_phase_diagram/density_phases copy.tab")                  #Reads the density values for the phases
    #Numerics---------------------------------------------------------------
    CFL                 = 50.0                                                                          #CFL condition
    res                 = [50 50;]                                                                      #Number of nodes
    resmin              = copy(res)                                                                     #Minimum number of nodes
    MRefin              = 2.0                                                                           #Refinement factor
    BCout               = [0 0]                                                                         #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann;
                                                                                                        #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Create data set--------------------------------------------------------
    #Create arrays X(T) using linear least squares
    coeff_up, coeff_do  = coeff_trans_line(eq_values)                                                   #Extract coefficients
    TMAX                = Tstart + 1000.0                                                               #Max temperature for T-array
    TMIN                = Tstop  - 1000.0                                                               #Min temperature for T-array
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
    #rho_left            = copy(rho_phases[:,3])
    #rho_left            = reshape(rho_left,nd1,nd1)                                                     #Density component A in [kg/m^3]
    #rho_right           = copy(rho_phases[:,4])
    #rho_right           = reshape(rho_right,nd1,nd1)                                                    #Density component B in [kg/m^3]
    rho_left            = ones(size(Xwm))
    rho_right           = ones(size(Xwm))
    #Create other arrays----------------------------------------------------
    R_left              = C_leftlin .* inv.((1.0 .- C_leftlin))                                         #Composition rate in phase A
    R_right             = C_rightlin .* inv.((1.0 .- C_rightlin))                                       #Composition rate in phase B
    KDlin               = R_left .* inv.(R_right)                                                       #Partition coefficient KD
    Tpath               = LinRange(Tstart,Tstop,10000)                                                  #Temperature path in [K]
    tpath               = LinRange(0,t_tot+1e-10,10000)                                                 #Time path in [s]
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
    D_l                 = D_001*(cos(alpha))^2 + D_010*(cos(beta))^2 + D_100*(cos(gamma))^2             #Effective diffusion coefficient Olivine after Crank (1975), p. 7
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
    x0_left, x0_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
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
    #-----------------------------------------------------------------------
    #Solving the moving boundary problem------------------------------------
    #anim = Animation()
    while t < t_tot
        #Update time--------------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        if t <= t_tot
            T  = linear_interpolation_1D(tpath,Tpath,t)
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
        D_l           = D_001*(cos(alpha))^2 + D_010*(cos(beta))^2 + D_100*(cos(gamma))^2               #Effective diffusion coefficient after Crank (1975), p. 7
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
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,verbose)
        #Post-Preprocessing-------------------------------------------------
        for iit in enumerate(1)
            Massnow     = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            Massnow2    = (trapezoidal_integration(x_left.^n*rho[1],C_left)+ trapezoidal_integration(x_right.^n*rho[2],C_right))/Ri[2]
            R_left     = C_left[end] .* inv((1.0 - C_left[end]))
            R_right    = C_right[1] .* inv((1.0 - C_right[1]))
            KDt        = R_left .* inv(R_right)
            T_pl       = copy(T)
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
        if plot_sim && it % 200 == 0
            #Plotting-------------------------------------------------------
            maxC = maximum([maximum(C_left),maximum(C_right)])
            Tp_min = (Tstop - 273.0) * 0.95
            Tp_max = (Tstart - 273.0)* 1.05
            first_val, last_val = values_between_known_indices!(Tlin.-273.0,KDlin,Tstart,Tstop)                    #CAUTION: works just for constantly dropping temperature

            #Composition profile
            p1 = plot(x_left*1e3,C_left, lw=2, label=L"\mathrm{Left\ side}")
            p1 = plot!(x_right*1e3,C_right, lw=2, label=L"\mathrm{Right\ side}")
            p1 = plot!(x0*1e3,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
                  ylabel = L"X_{Fe}", lw=1.5, grid=:on, legend = :right)
            p1 = plot!([x_left[end]; x_left[end];]*1e3, [0; 1*(maxC + 0.01)], color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}",ylim=[C0[1]-0.05; 1*(maxC + 0.01)])
            #Phase diagram
            p2 = plot(Tlin .- 273.0,XC_left, lw=2, label=L"\mathrm{Left\ side}")
            p2 = plot!(Tlin .- 273.0,XC_right, lw=2, label=L"\mathrm{Right\ side}")
            p2 = scatter!([T-273.0],[C_left[end]],marker=:circle, markersize=2, markercolor=:grey68,
                          markerstrokecolor=:grey68,label = "",ylabel = L"X_{Fe}",xlabel=L"T\ \mathrm{[°C]}")
            p2 = scatter!([T-273.0],[C_right[1]],marker=:circle, markersize=2, markercolor=:grey68,
                          markerstrokecolor=:grey68,label = "")
            p2 = scatter!([Tstart],[C0[51]],marker=:circle, markersize=2, markercolor=:black,
                          markerstrokecolor=:black,label = "",ylabel = L"X_{Fe}",xlabel=L"T\ \mathrm{[°C]}")
            p2 = scatter!([Tstart],[C0[50]],marker=:circle, markersize=2, markercolor=:black,
                          markerstrokecolor=:black,label = "")
            p2 = plot!([Tstart; Tstart],[0; maximum(C0[end])],lw=1.5,color=:black,linestyle=:dash,label=L"\mathrm{T(t=0.0)}")
            p2 = plot!([T-273.0; T-273.0],[0; maximum([C_left[end],C_right[1]])],lw=1.5,color=:grey68,linestyle=:dashdot,label=L"\mathrm{T(t_{tot})}")
            p2 = plot!([T-273.0; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:royalblue,linestyle =:dot)
            p2 = plot!([T-273.0; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), ylims=(0, 1),color=:crimson,linestyle =:dot)
            #p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass2[end]; Mass2[end]],color=:dimgrey,linestyle=:dashdot,lw=1.5, label=L"\mathrm{Final\ mass}")
            #p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass01; Mass01],color=:grey,linestyle=:dashdot,lw=1.5, label=L"\mathrm{Initial\ mass}",
            #            xlabel = L"T\ \mathrm{[°C]}", ylabel = L"X_{Mg}",grid=:on,legend = :topright)
            #p2 = scatter!([T_check],[C_left_check],marker=:circle, markersize=2, markercolor=:black,
            #              markerstrokecolor=:black,label = "Model")
            #p2 = scatter!([T_check],[C_right_check],marker=:circle, markersize=2, markercolor=:black,
            #              markerstrokecolor=:black,label = "mODEL")
            #Evolution of KD(T)
            p3 = plot(Tlin .- 273.0, KDlin, lw=1.5, label=L"\mathrm{Thermodynamic\ data}", color=:black)
            p3 = scatter!([T_sim[end]-273.0],[KD_sim[end]],marker=:circle, markersize=3.0, markercolor=:black,markerstrokecolor=:black,
                          xlabel = L"T\ \mathrm{[°C]}", ylabel = L"K_{D}", lw=1.5,legend = :bottomleft,
                          grid=:on, label=L"\mathrm{Model}",xlims=(Tp_min, Tp_max), ylims=(first_val-0.01,last_val+0.01))
            #ln(KD) vs 1/T
            p4 = plot(10000.0 ./ (T_sim),log.(KD_sim),xlabel = L"10,000/T\ \mathrm{[K^{-1}]}", ylabel = L"ln(K_{D})", lw=1.5,
                    grid=:on, label="", color=:black, ticks=:auto, xrotation=0)
            #Figure 1
            p = plot(p1,p2,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                    legend_foreground_color = :transparent)
            display(p)
            #frame(anim)
            #Figure 2
            #plot(p3,p4,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
            #        legend_foreground_color = :transparent)
        end
        # Suppress output of calc_mass_err
        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error,ErrM)
        end
    end
    #Post-process-----------------------------------------------------------
    #gif(anim, "figures/D1.gif", fps=5)  # Save with 10 frames per second
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, vec(C_left), vec(C_right), vec(C0),maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, KDlin, KD_sim,T_sim, Mass0, Mass, Mass01, Mass2, C_left_check, C_right_check, T_check, Residual, MB_Error
end
#Run calculation------------------------------------------------------------
run_and_plot = true
run_and_plot == false ? printstyled("You have disabled the simulation, change the variable run_and_plot == true", bold=true) : nothing
if run_and_plot
    plot_sim  = false
    plot_end  = true
    verbose   = false
    save_file = false
    x_left, x_right, x0, C_left, C_right, C0, maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, KDlin, KD_sim,T_sim, Mass0, Mass, Mass01, Mass2, C_left_check, C_right_check, T_check,Residual, MB_Error = D1(; plot_sim = plot_sim, verbose = verbose)
    if plot_end
        #Title: Thermodynamical constrained Stefan condition
        #Plotting in K-----------------------------------------------------------
        Tstart = Tstart - 273.0
        Tstop  = Tstop - 273.0
        Tp_min = Tstop * 0.95
        Tp_max = Tstart * 1.05
        first_val, last_val = values_between_known_indices!(Tlin.-273.0,KDlin,Tstart,Tstop)                    #CAUTION: works just for constantly dropping temperature
        fs = 12.0
        #Composition profile
        p1 = plot(x_left*1e3,C_left, lw=2, label=L"\mathrm{Left\ side}")
        p1 = plot!(x_right*1e3,C_right, lw=2, label=L"\mathrm{Right\ side}")
        p1 = plot!(x0*1e3,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
              ylabel = L"X_{Fe}", lw=1.5, grid=:on, legend = :right)
        p1 = plot!([x_left[end]; x_left[end];]*1e3, [0; 1*(maxC + 0.01)], color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}",ylim=[C0[1]-0.05; 1*(maxC + 0.01)])
        #Phase diagram
        p2 = plot(Tlin .- 273.0,XC_left, lw=2, label=L"\mathrm{Left\ side}")
        p2 = plot!(Tlin .- 273.0,XC_right, lw=2, label=L"\mathrm{Right\ side}")
        p2 = scatter!([T-273.0],[C_left[end]],marker=:circle, markersize=2, markercolor=:grey68,
                      markerstrokecolor=:grey68,label = "",ylabel = L"X_{Fe}",xlabel=L"T\ \mathrm{[°C]}")
        p2 = scatter!([T-273.0],[C_right[1]],marker=:circle, markersize=2, markercolor=:grey68,
                      markerstrokecolor=:grey68,label = "")
        p2 = scatter!([Tstart],[C0[51]],marker=:circle, markersize=2, markercolor=:black,
                      markerstrokecolor=:black,label = "",ylabel = L"X_{Fe}",xlabel=L"T\ \mathrm{[°C]}")
        p2 = scatter!([Tstart],[C0[50]],marker=:circle, markersize=2, markercolor=:black,
                      markerstrokecolor=:black,label = "")
        p2 = plot!([Tstart; Tstart],[0; maximum(C0[end])],lw=1.5,color=:black,linestyle=:dash,label=L"\mathrm{T(t=0.0)}")
        p2 = plot!([T-273.0; T-273.0],[0; maximum([C_left[end],C_right[1]])],lw=1.5,color=:grey68,linestyle=:dashdot,label=L"\mathrm{T(t_{tot})}")
        p2 = plot!([T-273.0; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:royalblue,linestyle =:dot)
        p2 = plot!([T-273.0; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), ylims=(0, 1),color=:crimson,linestyle =:dot)
        #p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass2[end]; Mass2[end]],color=:dimgrey,linestyle=:dashdot,lw=1.5, label=L"\mathrm{Final\ mass}")
        #p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass01; Mass01],color=:grey,linestyle=:dashdot,lw=1.5, label=L"\mathrm{Initial\ mass}",
        #            xlabel = L"T\ \mathrm{[°C]}", ylabel = L"X_{Mg}",grid=:on,legend = :topright)
        #p2 = scatter!([T_check],[C_left_check],marker=:circle, markersize=2, markercolor=:black,
        #              markerstrokecolor=:black,label = "Model")
        #p2 = scatter!([T_check],[C_right_check],marker=:circle, markersize=2, markercolor=:black,
        #              markerstrokecolor=:black,label = "mODEL")
        #Evolution of KD(T)
        p3 = plot(Tlin .- 273.0, KDlin, lw=1.5, label=L"\mathrm{Thermodynamic\ data}", color=:black)
        p3 = scatter!([T_sim[end]-273.0],[KD_sim[end]],marker=:circle, markersize=3.0, markercolor=:black,markerstrokecolor=:black,
                      xlabel = L"T\ \mathrm{[°C]}", ylabel = L"K_{D}", lw=1.5,legend = :bottomleft,
                      grid=:on, label=L"\mathrm{Model}",xlims=(Tp_min, Tp_max), ylims=(first_val-0.01,last_val+0.01))
        #ln(KD) vs 1/T
        p4 = plot(10000.0 ./ (T_sim),log.(KD_sim),xlabel = L"10,000/T\ \mathrm{[K^{-1}]}", ylabel = L"ln(K_{D})", lw=1.5,
                grid=:on, label="", color=:black, ticks=:auto, xrotation=0)
        #Figure 1
        plot(p1,p2,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                legend_foreground_color = :transparent)
        display(current())
        save_path = "figures"
        save_name = "D1"
        save_figure(save_name,save_path,save_file)
        #Figure 2
        #plot(p3,p4,dpi = 300,legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
        #        legend_foreground_color = :transparent)
        #display(current())
        #save_name = "D1_KD"
        #save_figure(save_name,save_path,save_file)
    end
end
