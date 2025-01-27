#= Utilized functions in the framework of the Paper:

    by A.Stroh and E. Moulas
    doi:
    Version: 1.0
=#
using Diff_Coupled
using Plots, LinearAlgebra, Revise, LaTeXStrings
#Main function-------------------------------------------------
function main(plot_sim,verbose)
    #If you find a [] with two entires this belong to the respective side of the diffusion couple ([left right])
    #Phyics-------------------------------------------------------
    #Di      = [1.0    1.0]
    Di      = [2.65*1e-18   6.23*1e-20;]#Initial diffusion coefficient in [m^2/s]; 
                                    #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [2.75*1e-6    3.9*1e-7;]#Pre-exponential factor in [m^2/s]
    rho     = [1.0      1.0;]       #Normalized densities in [g/mol]
    Ri      = [0.0002    0.1;]      #Initial radii [interface    total length] in [m]
    Cl_i    = 0.6                   #Initial concentration left side in [mol]
    Cr_i    = 0.3                   #Initial concentration right side in [mol]
    V_ip    = 3.17e-14              #Interface velocity in [m/s]
    R       = 8.314472              #Universal gas constant in [J/(mol*K)]
    Ea1     = 292879.6767           #Activation energy for the left side in [J/mol]
    Ea2     = 360660.4018           #Activation energy for the right side in [J/mol]
    Myr2Sec = 60*60*24*365.25*1e6   #Conversion factor from Myr to s
    t_tot   = 1.0e-2 * Myr2Sec      #Total time [s]
    n       = 1                     #Geometry; 1: planar, 2: cylindric, 3: spherical
    #History dependent parameters-----------------------------------
    KD_ar   = LinRange(0.9,0.9,1000)        #Partition coefficient array to calculate partition coefficient history; KD changes with respect to time;
                                            #The last value must be equal to the partition coefficient at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)      #Time array (in s) to calculate history over time. The last value must be equal to t_tot.
                                            #The user is prompted to specify suitable time intervals in relation to the respective destination.               
    T_ar    = LinRange(1273.15,923.15,1000) #Temperature arrray in [K] to calculate temperature history; T changes with respect to time; 
                                            #The last value must be equal to the temperature at t = t_tot.
    #Numerics-----------------------------------------------------
    CFL    = 0.5                    #CFL condition
    res    = [100 150;]             #Number of grid points
    resmin = copy(res)              #Minimum number of grid points
    MRefin = 5.0                    #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
    BCout  = [0 0]                  #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann; 
                                    #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Check, if t_ar is valid (increasing in time)----------------------------------------
    dt_diff = zeros(length(t_ar)-1)
    dt_diff = t_ar[2:end] .- t_ar[1:end-1]
    if any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    end
    #Create mesh, discretization and mapping------------------------
    if Ri[1] >= Ri[2]
        error("Please change the size of the system. Increase Ri[2].")
    elseif res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    x_left, x_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    #Preprocess and initial condition--------------------------------------------
    L       = Ri[end]                                   #Length of the domain in [m]
    t       = 0.0                                       #Initial time in [s]
    it      = 0                                         #Initial number of time iterations                       
    C_left  = Cl_i*ones(res[1],1)                       #Initial concentration left side in [mol]
    C_right = Cr_i*C_left[end]*ones(res[2],1)*inv(KD_ar[1])  #Initial concentration right side in [mol]
    C0      = [copy(C_left); copy(C_right)]             #Store initial concentration 
    C       = copy(C0)                                  #Create 1 array with all concentrations  
    C0_l    = copy(C_left)                              #Store initial concentration left side
    C0_r    = copy(C_right)                             #Store initial concentration right side
    x       = copy(x0)                                  #Create 1 array containing all x-values  
    Ri0     = copy(Ri)                                  #Store initial radii
    KD      = copy(KD_ar[1])                            #Initial partition coefficient, just for pre-processing
    #Total mass---------------------------------------------------
    Mass0   = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)        #Initial mass
    #Preallocate variables----------------------------------------
    Co, Co_l, Co_r, dt, dx, Kloc, Lloc, L_g, Mass, Mloc, nels, nels_l, nels_r, R_g, x_1, x_2, y_interp = preallocations(x, C, C_left, C_right,res)
    #Calculate initial Ds, KD, T---------------------------------------------
    D_l, D_r, KD, T = update_t_dependent_param!(D0,Di,dt,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
    #First check for correct setup-----------------------------------------
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
    while t < t_tot
        #Calculate dt-------------------------------------------------
        dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
        #Update time---------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters-------------------------------
        D_l, D_r, KD,T = update_t_dependent_param!(D0,Di,dt,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
        #Advect interface & regrid--------------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,res)
        #FEM SOLVER-----------------------------------------------------
        #Construct global matrices--------------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,nels_l,nels_r,Mloc,Kloc,Lloc,res)
        #Set inner boundary conditions----------------------------------
        L_g, R_g, ScF = set_inner_bc_flux!(L_g,R_g,KD,D_l,D_r,x_left,x_right,V_ip,rho,res)
        #Set outer boundary conditions and scale matrices---------------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #Solve system---------------------------------------------------
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Regrid---------------------------------------------------------
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,verbose)
        #Post-Preprocessing---------------------------------------------
        for iit in enumerate(1)       
            Massnow = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            push!(Mass, Massnow)  #Stores the mass of the system
        end
        if plot_sim
            #Plotting------------------------------------------------------
            p = plot(x_left,C_left, lw=2, label=L"Left\ side")
            p = plot!(x_right,C_right, lw=2, label=L"Right\ side")
            p = plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance", ylabel = L"Concentration", title = L"Diffusion\ couple\ -\ growth\ (flux)", lw=1.5,
                    grid=:on, label=L"Initial\ condition")
            #plot!([Ri[1]; Ri[1]], [0; 1]*maxC, color=:grey68, lw=2,label=L"Interface", linestyle=:dashdot)
            display(p)
        end
        println("it $it")
        println("")
    end
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, dx1, dx2, x0, res, Ri, C_left, C_right, C0, maxC
end
#Call main function-------------------------------------------------------------
run_and_plot = false
if run_and_plot
    plot_sim = false
    plot_end = true
    verbose  = false
    x_left, x_right, dx1, dx2, x0, res, Ri, C_left, C_right, C0, maxC = main(plot_sim,verbose)
    if plot_end    
        #Plotting------------------------------------------------------
        plot(x_left,C_left, lw=2, label=L"Left\ side")
        plot!(x_right,C_right, lw=2, label=L"Right\ side")
        plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance", ylabel = L"Concentration", title = L"Diffusion\ couple\ -\ growth\ (flux)", lw=1.5,
              grid=:on, label=L"Initial\ condition")
        #plot!([Ri[1]; Ri[1]], [0; 1]*maxC, color=:grey68, lw=2,label=L"Interface", linestyle=:dashdot)
    end
end