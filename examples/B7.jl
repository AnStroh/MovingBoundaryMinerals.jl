using Diff_Coupled, Diff_Coupled.Benchmarks
using Plots, LinearAlgebra, Revise, LaTeXStrings
#Main function----------------------------------------------------
function main(plot_sim,verbose)
    #If you find a [] with two entires this belong to the respective side of the diffusion couple ([left right])
    #Physics-------------------------------------------------------
    #Di      = [-1.0   -1.0;]
    Di      = [1.0   1.0]                                    #Initial diffusion coefficient in [m^2/s]           -> in [L*V]
                                                                #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [1e-4   5e-4;]                                    #Pre-exponential factor in [m^2/s]                  -> NOT USED
    rho     = [1.0      1.0;]                                   #Normalized densities in [-]                   -> NOT USED
    Ri      = [0.5       1.0;]                                   #Initial radii [interface    total length] in [m]   -> in [L]
    Cl_i    = 0.5                                               #Initial concentration left side in [mol]           -> in [C]
    Cr_i    = Cl_i/0.1                                         #Initial concentration right side in [mol]          -> -//-
    V_ip    = -0.01                                                 #Interface velocity in [m/s]                        -> in [V]
    R       = 8.314472                                          #Universal gas constant in [J/(mol*K)]              -> NOT USED
    Ea1     = 292879.6767                                       #Activation energy for the left side in [J/mol]     -> NOT USED
    Ea2     = 300660.4018                                       #Activation energy for the right side in [J/mol]    -> NOT USED
    Myr2Sec = 60*60*24*365.25*1e6                               #Conversion factor from Myr to s                    -> NOT USED
    t_tot   = 1.0                                               #Total time [s]                                     -> in [L]/[V]
    n       = 1                                                 #Geometry; 1: planar, 2: cylindrical, 3: spherical
    #History dependent parameters---------------------------------
    KD_ar   = LinRange(Cl_i/Cr_i,Cl_i/Cr_i*0.5,1000)            #Partition coefficient array to calculate partition coefficient history; KD changes with respect to time;
                                                                #The last value must be equal to the partition coefficient at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)                          #Time array (in s) to calculate history over time. The last value must be equal to t_tot.
                                                                #The user is prompted to specify suitable time intervals in relation to the respective destination.               
    T_ar    = LinRange(1273.15,1073.15,1000)                    #Temperature arrray in [K] to calculate temperature history; T changes with respect to time; 
                                                                #The last value must be equal to the temperature at t = t_tot.
    #Numerics-----------------------------------------------------
    CFL    = 0.4                                                #CFL condition
    res    = [5 6;]                                         #Number of grid points
    resmin = copy(res)                                          #Minimum number of grid points
    MRefin = 1.0                                               #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
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
    println("size of x_left: ", size(x_left), " size of x_right: ", size(x_right))
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
    for i in 1:1
    #while t < t_tot
        #Calculate dt---------------------------------------------
        dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
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
        #@show L_g[5,:]
        println("det(L_g): ", det(L_g))
        println("i: ", i)
        @show dt*V_ip
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
    Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar = rescale(Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar, Lsc, Dsc, Vsc, tsc)    
    #Post-process------------------------------------------------
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC
end
#Call main function-----------------------------------------------
run_and_plot = true
if run_and_plot
    plot_sim  = true
    plot_end  = true
    verbose   = false
    save_file = false
    x_left, x_right, x0, Ri, Ri0, C_left, C_right, C0, C0_r, KD0, n, maxC = main(plot_sim,verbose)
        if plot_end
        #Plotting-------------------------------------------------
        p1 = plot(x_left,C_left, lw=2, label=L"Left\ side")
        p1 = plot!(x_right,C_right, lw=2, label=L"Right\ side")
        p1 = plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance\ [m]", ylabel = L"Concentration", lw=1.5,
                   grid=:on, label=L"Initial\ condition",legendfontsize = 6, dpi = 300)
        p1 = plot!([Ri[1]; Ri[1]], [0; 1]*maxC,title = L"Diffusion\ couple\ (flux)\ -\ Concentration\ profile", color=:grey68,linestyle=:dashdot, lw=2,label=L"Interface")
        #save_path = "figures"
        #save_name = "B7"
        #save_figure(save_name,save_path,save_file)
    end
end