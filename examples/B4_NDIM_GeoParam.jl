using Diff_Coupled, Diff_Coupled.Benchmarks
using Plots, LinearAlgebra, Revise, LaTeXStrings, GeoParams
#Main function----------------------------------------------------
function main(plot_sim,verbose)
    CharUnits = GEO_units(;length=1e-2m, temperature = 1273.15K)
    #If you find a [] with two entires this belong to the respective side of the diffusion couple ([left right])
    #Physics-------------------------------------------------------
    Di_dim      = [1e-4   1e4;]                                 #Initial diffusion coefficient in [m^2/s]           -> in [L*V]
                                                                #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0_dim      = [9999   99999;]                               #Pre-exponential factor in [m^2/s]                  -> NOT USED
    rho_dim     = [1.0      1.0;]                               #Normalized densities in [-]                   -> NOT USED
    Ri_dim      = [1e-2       1;]                               #Initial radii [interface    total length] in [m]   -> in [L]
    Cl_i_dim    = 0.5                                           #Initial concentration left side in [mol]           -> in [C]
    Cr_i_dim    = Cl_i_dim/100                                  #Initial concentration right side in [mol]          -> -//-
    V_ip_dim    = 1.0                                           #Interface velocity in [m/s]                        -> in [V]
    R_dim       = 8.314472                                      #Universal gas constant in [J/(mol*K)]              -> NOT USED
    Ea1_dim     = 292879.6767                                   #Activation energy for the left side in [J/mol]     -> NOT USED
    Ea2_dim     = 360660.4018                                   #Activation energy for the right side in [J/mol]    -> NOT USED
    Myr2Sec     = 60*60*24*365.25*1e6                           #Conversion factor from Myr to s                    -> NOT USED
    t_tot_dim   = 0.35 * Myr2Sec                                #Total time [s]                                     -> in [L]/[V]
    n           = 3                                             #Geometry; 1: planar, 2: cylindrical, 3: spherical
    
    #History dependent parameters---------------------------------
    KD_ar   = LinRange(Cl_i_dim/Cr_i_dim,Cl_i_dim/Cr_i_dim,1000)#Partition coefficient array to calculate partition coefficient history; KD changes with respect to time;
                                                                #The last value must be equal to the partition coefficient at t = t_tot.
    t_ar    = LinRange(0.0,t_tot_dim,1000)                      #Time array (in s) to calculate history over time. The last value must be equal to t_tot.
                                                                #The user is prompted to specify suitable time intervals in relation to the respective destination.               
    T_ar    = LinRange(1273.15,1273.15,1000)                    #Temperature arrray in [K] to calculate temperature history; T changes with respect to time; 
                                                                #The last value must be equal to the temperature at t = t_tot.
    #Phyics nondimensionalized -----------------------------------
    Di      = nondimensionalize(Di_dim.*(m^2/s), CharUnits)     
    D0      = nondimensionalize(D0_dim.*m^2/s, CharUnits)     
    rho     = nondimensionalize(rho_dim.*kg/m^3, CharUnits)                                  
    Ri      = nondimensionalize(Ri_dim.*m, CharUnits)                                  
    Cl_i    = nondimensionalize(Cl_i_dim*mol,CharUnits)                                               
    Cr_i    = nondimensionalize(Cr_i_dim*mol,CharUnits)                                               
    V_ip    = nondimensionalize(V_ip_dim*m/s, CharUnits)                                          
    R       = nondimensionalize(R_dim * J/(mol*K), CharUnits)                                         
    Ea1     = nondimensionalize(Ea1_dim*J/mol, CharUnits)                                       
    Ea2     = nondimensionalize(Ea2_dim*J/mol, CharUnits)                                       
    t_tot   = nondimensionalize(t_tot_dim*s, CharUnits)                               
    t_ar    = LinRange(nondimensionalize(0.0s,CharUnits),t_tot,1000)                        
    T_ar    = LinRange(nondimensionalize(1273.15K, CharUnits),nondimensionalize(923.15K,CharUnits),1000)
    #Numerics-----------------------------------------------------
    CFL    = 0.4                                                #CFL condition
    res    = [80 120;]                                          #Number of grid points
    resmin = copy(res)                                          #Minimum number of grid points
    MRefin = 50.0                                               #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
    BCout  = [0 0]                                              #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann; 
                                                                #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Non-dimensionslization---------------------------------------
    #V_ip, t_tot, t_ar, Di, D0, Ri, Lsc, Dsc, Vsc, tsc = scaling(Ri, Di, D0, V_ip, t_tot, t_ar)
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
    C_left  = Cl_i*ones(res[1],1)                               #Initial concentration left side in [mol]
    C_right = Cr_i*ones(res[2],1)                               #Initial concentration right side in [mol]
    C0      = [copy(C_left); copy(C_right)]                     #Store initial concentration 
    C       = copy(C0)                                          #Create 1 array with all concentrations  
    C0_l    = copy(C_left)                                      #Store initial concentration left side
    C0_r    = copy(C_right)                                     #Store initial concentration right side
    x       = copy(x0)                                          #Create 1 array containing all x-values  
    Ri0     = copy(Ri)                                          #Store initial radii
    KD      = copy(KD_ar[1])                                    #Initial partition coefficient, just for pre-processing
    #Total mass---------------------------------------------------
    Mass0   = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)        
    #Preallocate variables----------------------------------------
    Co, Co_l, Co_r, dt, dx, Kloc, Lloc, L_g, Mass, Mloc, nels, nels_l, nels_r, R_g, x_1, x_2, y_interp = preallocations(x, C, C_left, C_right,res)
    #Calculate initial Ds, KD, T----------------------------------
    D_l, D_r, KD, T = update_t_dependent_param!(D0,Di,dt,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)

    KD0 = copy(KD)
    #First check for correct setup--------------------------------
    if BCout[1] != 0 && (n == 3 || n == 2)
        error("The code is only valid for cylindrical/spherical geometry, where the left outer BC has Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T  != T_ar[1]   
        error("Initial temperature must be equal to the first value in the temperature array.")
    end
    #Time loop----------------------------------------------------
    #for i in 1:1
    while t < t_tot
        #Calculate dt---------------------------------------------
        #dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
        dt  = minimum([dx1,dx2])/V_ip*CFL
        #Update time----------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters-------------------------
        D_l, D_r, KD,T = update_t_dependent_param!(D0,Di,dt,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
        #Advect interface & regrid--------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,res)
        #FEM SOLVER-----------------------------------------------
        #Construct global matrices--------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,nels_l,nels_r,Mloc,Kloc,Lloc,res)
        #Set inner boundary conditions----------------------------
        L_g, R_g, ScF = set_inner_bc_flux!(L_g,R_g,KD,D_l,D_r,x_left,x_right,V_ip,rho,res)
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
        if plot_sim
            #Plotting---------------------------------------------
            maxC = maximum([maximum(C_left),maximum(C_right)])
            p1 = plot(x_left,C_left, lw=2, label=L"Left\ side")
            p1 = plot!(x_right,C_right, lw=2, label=L"Right\ side")
            p1 = plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance", ylabel = L"Concentration", lw=1.5,
                       grid=:on, label=L"Initial\ condition",legendfontsize = 6,title = L"Diffusion\ couple\ (flux)\ -\ Rayleigh\ fractionation")
            p1 = plot!([Ri[1]; Ri[1]], [0; 1]*maxC,title = L"Concentration\ profile", color=:grey68,linestyle=:dashdot, lw=2,label=L"Interface")
            display(p1)
        end
    end
    #Rescaling---------------------------------------------------
    #Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar = rescale(Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar, Lsc, Dsc, Vsc, tsc)    
    x_left  = ustrip(dimensionalize(x_left, m, CharUnits))
    x_right = ustrip(dimensionalize(x_right, m, CharUnits))
    dx1     = ustrip(dimensionalize(dx1, m, CharUnits))
    dx2     = ustrip(dimensionalize(dx2, m, CharUnits))
    x0      = ustrip(dimensionalize(x0, m, CharUnits))
    Ri      = ustrip(dimensionalize(Ri, m, CharUnits))
    C_left  = ustrip(dimensionalize(C_left, mol, CharUnits))
    C_right = ustrip(dimensionalize(C_right, mol, CharUnits))
    C0      = ustrip(dimensionalize(C0, mol, CharUnits))
    #Post-process------------------------------------------------
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC
end
#Call main function-----------------------------------------------
run_and_plot = true
if run_and_plot
    plot_sim  = false
    plot_end  = true
    verbose   = false
    save_file = false
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC = main(plot_sim,verbose)
    @show x_left,C_left,Ri0,Ri
    Ray_Fs, Ray_Fl, Ray_Cl, Ray_Cs, Cl_p, phi_solid = rayleigh_fractionation(x_left,C_left,Ri0,Ri,C0_r,KD0,n)
    if plot_end
        #Plotting-------------------------------------------------
        p1 = plot(x_left,C_left, lw=2, label=L"Left\ side")
        p1 = plot!(x_right,C_right, lw=2, label=L"Right\ side")
        p1 = plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance\ [m]", ylabel = L"Concentration", lw=1.5,
                   grid=:on, label=L"Initial\ condition",legendfontsize = 6)
        p1 = plot!([Ri[1]; Ri[1]], [0; 1]*maxC,title = L"Concentration\ profile", color=:grey68,linestyle=:dashdot, lw=2,label=L"Interface")
        p2 = plot((x_left/Ri0[2]).^(n),C_left, lw=2, label=L"Solid")
        p2 = plot!((x_right/Ri0[2]).^(n),C_right, lw=2, label=L"Liquid")
        p2 = scatter!([Ray_Fs[1:100:end]],[Ray_Cs[1:100:end]], marker=:circle, markersize=2, markercolor=:midnightblue, markerstrokecolor=:midnightblue,label=L"Solid\ Rayleigh", 
                      xlabel = L"Fraction", ylabel = L"Concentration", grid=:on,legendfontsize = 6,
                      title = L"Fractions")
        p2 = scatter!([Ray_Fs[end]],[Ray_Cs[end]], marker=:circle, markersize=2, markercolor=:midnightblue, markerstrokecolor=:midnightblue, label="",xlim=(x_left[1]-0.0001*Ri[2], Ri[1]+0.0001*Ri[2]))
        p3 = plot((phi_solid.^n)', Cl_p', lw=2, label=L"Mass\ fraction\ solid")
        p3 = scatter!([Ray_Fs[1:100:end]],[Ray_Cs[1:100:end]], marker=:circle, markersize=2, markercolor=:midnightblue, markerstrokecolor=:midnightblue,label=L"Solid\ Rayleigh",
                      xlabel = L"Fraction", ylabel = L"Concentration", grid=:on,legendfontsize = 6,
                      title = L"Concentration\ profile\ solid")
        p3 = scatter!([Ray_Fs[end]],[Ray_Cs[end]], marker=:circle, markersize=2, markercolor=:midnightblue, markerstrokecolor=:midnightblue, label="")
        plot(p1,p2,p3,suptitle = L"Diffusion\ couple\ (flux)\ -\ Rayleigh\ fractionation", dpi = 300)
        #save_path = "figures"
        #save_name = "B4"
        #save_figure(save_name,save_path,save_file)
    end
end