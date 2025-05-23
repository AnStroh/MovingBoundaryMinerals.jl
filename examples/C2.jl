using MOBILE, MOBILE.Benchmarks
using Plots, LinearAlgebra, LaTeXStrings,SparseArrays
#Main function----------------------------------------------------
function main(plot_sim,verbose)
    #If you find a [] with two entries this belong to the respective side of the diffusion couple ([left right])
    #Physics-------------------------------------------------------
    #Note: this example uses non-dimensional numbers! Given units might not be true!
    Di      = [0.001   0.005]                                   #Initial diffusion coefficient in [m^2/s]           -> in [L*V]
                                                                #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [NaN     NaN;]                                    #Pre-exponential factor in [m^2/s]                  -> not used in this example
    rho     = [1.0     1.0;]                                    #Normalized densities in [-]                        -> not used in this example
    Ri      = [5.0     10;]                                     #Initial radii [interface    total length] in [m]   -> in [L]
    Cl_i    = 0.1                                               #Initial composition left side in [-]               -> in [C]
    Cr_i    = Cl_i/50.0                                         #Initial composition right side in [-]              -> -//-
    V_ip    = -5e-3                                             #Interface velocity in [m/s]                        -> in [V]
    R       = NaN                                               #Universal gas constant in [J/(mol*K)]              -> not used in this example
    Ea1     = NaN                                               #Activation energy for the left side in [J/mol]     -> not used in this example
    Ea2     = NaN                                               #Activation energy for the right side in [J/mol]    -> not used in this example
    Myr2Sec = 60*60*24*365.25*1e6                               #Conversion factor from Myr to s                    -> not used in this example
    t_tot   = 3e2                                               #Total time [s]                                     -> in [L]/[V]
    n       = 3                                                 #Geometry; 1: planar, 2: cylindrical, 3: spherical
    #History dependent parameters---------------------------------
    KD_ar   = LinRange(Cl_i/Cr_i,Cl_i/Cr_i*0.5,1000)            #Partition coefficient array to calculate distribution coefficient history; KD changes with respect to time;
                                                                #The last value must be equal to the distribution coefficient at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)                          #Time array (in [s]) to calculate history over time. The last value must be equal to t_tot.
                                                                #The user is prompted to specify suitable time intervals in relation to the respective destination.
    T_ar    = LinRange(1273.15,1073.15,1000)                    #Temperature array in [K] to calculate temperature history; T changes with respect to time;
                                                                #The last value must be equal to the temperature at t = t_tot.
    #Numerics-----------------------------------------------------
    CFL    = 0.4                                                #CFL condition
    res    = [80 120;]                                          #Number of nodes
    resmin = copy(res)                                          #Minimum number of nodes
    MRefin = 50.0                                               #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
    BCout  = [0 0]                                              #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann;
                                                                #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Non-dimensionslization---------------------------------------
    V_ip, t_tot, t_ar, Di, D0, Ri, Lsc, Dsc, Vsc, tsc = scaling(Ri, Di, D0, V_ip, t_tot, t_ar)
    #Check, if t_ar is valid (increasing in time)-----------------
    dt_diff = zeros(length(t_ar)-1)
    dt_diff = t_ar[2:end] .- t_ar[1:end-1]
    if any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    end
    #Create mesh, discretization and mapping----------------------
    if Ri[1] >= Ri[2]
        error("Please change the size of the system. Increase Ri[2].")
    elseif res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    x_left, x_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    #Preprocess and initial condition-----------------------------
    L       = Ri[end]                                           #Length of the domain in [m]
    t       = 0.0                                               #Initial time in [s]
    it      = 0                                                 #Initial number of time iterations
    C_left  = Cl_i*ones(res[1],1)                               #Initial composition left side in [-]
    C_right = Cr_i*ones(res[2],1)                               #Initial composition right side in [-]
    C0      = [copy(C_left); copy(C_right)]                     #Store initial composition
    C       = copy(C0)                                          #Create 1 array with all concentrations
    C0_l    = copy(C_left)                                      #Store initial composition left side
    C0_r    = copy(C_right)                                     #Store initial composition right side
    x       = copy(x0)                                          #Create 1 array containing all x-values
    Ri0     = copy(Ri)                                          #Store initial radii
    KD      = copy(KD_ar[1])                                    #Initial distribution coefficient, just for pre-processing
    KD0     = copy(KD)                                          #Store initial distribution coefficient
    #Total mass---------------------------------------------------
    Mass0   = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
    #Preallocate variables----------------------------------------
    Co_l    = zeros(size(C_left))                               #Matrix to store old concentrations of left side
    Co_r    = zeros(size(C_right))                              #Matrix to store old concentrations of right side
    dt      = 0.0                                               #Initial time step
    L_g     = spzeros(length(x),length(x))                      #Global left hand side matrix
    Mass    = Float64[]                                         #Array to store the mass of the system
    R_g     = zeros(length(x),1)                                #Global right hand side vector
    #Checks-------------------------------------------------------
    MB_Error = Float64[]                                        #Array to store the mass error
    #Calculate initial Ds, KD, T----------------------------------
    D_l, D_r, KD, T = update_t_dependent_param!(D0,Di,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
    #First check for correct setup--------------------------------
    if BCout[1] != 0 || BCout[2] != 0
        error("The code is only valid for a closed system. Please set outer BC to Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T  != T_ar[1]
        error("Initial temperature must be equal to the first value in the temperature array.")
    end
    #Time loop----------------------------------------------------
    #anim = Animation()
    while t < t_tot
        #Calculate dt---------------------------------------------
        dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
        #Update time----------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters-------------------------
        D_l, D_r, KD,T = update_t_dependent_param!(D0,Di,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
        #Advect interface & regrid--------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,res)
        #Calculate current volume---------------------------------
        Vol_left, Vol_right, dVolC = calc_volume(x_left,x_right,n)
        #FEM SOLVER-----------------------------------------------
        #Construct global matrices--------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,res)
        #Set inner boundary conditions----------------------------
        L_g, R_g, ScF = set_inner_bc_mb!(L_g,R_g,dVolC,Mass0,KD,res)
        #Set outer boundary conditions and scale matrices---------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #Solve system---------------------------------------------
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Regrid---------------------------------------------------
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,verbose)
        #Post-Preprocessing---------------------------------------
        for iit in enumerate(1)
            Massnow = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            push!(Mass, Massnow)                                #Stores the mass of the system
        end
        if plot_sim && it % 20 == 0
            #Plotting---------------------------------------------
            maxC = maximum([maximum(C_left),maximum(C_right)])
            fs = 12.0
            p = plot(x_left,C_left, lw=2, label=L"\mathrm{Left\ side}")
            p = plot!(x_right,C_right, lw=2, label=L"\mathrm{Right\ side}")
            p = plot!(x0,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
                    ylabel = L"C\ \mathrm{[-]}", lw=1.5, grid=:on,dpi = 300,
                        legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                        legend_foreground_color = :transparent)
            p = plot!([Ri[1]; Ri[1]], [0; 1]*maxC, color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}")
            display(p)
            #frame(anim)
        end
        # Suppress output of calc_mass_err
        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error,ErrM)
        end
    end
    #Rescaling---------------------------------------------------
    Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar = rescale(Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar, Lsc, Dsc, Vsc, tsc)
    #Post-process------------------------------------------------
    #gif(anim, "figures/C2.gif", fps=5)  # Save with 10 frames per second
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, n, maxC, KD0
end
#Call main function-----------------------------------------------
run_and_plot = false
run_and_plot == false ? printstyled("You have disabled the simulation, change the variable run_and_plot == true", bold=true) : nothing
if run_and_plot
    plot_sim  = false
    plot_end  = true
    verbose   = false
    save_file = false
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, n, maxC, KD0 = main(plot_sim,verbose)
    if plot_end
        #Title: Diffusion couple (MB) - resorption + diffusion in a sphere
        #Plotting-------------------------------------------------
        fs = 12.0
        plot(x_left,C_left, lw=2, label=L"\mathrm{Left\ side}")
        plot!(x_right,C_right, lw=2, label=L"\mathrm{Right\ side}")
        plot!(x0,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[m]}",
              ylabel = L"C\ \mathrm{[-]}", lw=1.5, grid=:on,dpi = 300,
                    legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                    legend_foreground_color = :transparent)
        plot!([Ri[1]; Ri[1]], [0; 1]*maxC .+ 0.01, color=:grey68,linestyle=:dashdot, lw=2,label=L"\mathrm{Interface}")
        save_path = "figures"
        save_name = "C2"
        save_figure(save_name,save_path,save_file)
    end
end
