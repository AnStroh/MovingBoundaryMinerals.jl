#= Utilized functions in the framework of the Paper:

    by A.Stroh and E. Moulas
    doi:
    Version: 1.0
=#
using Diff_Coupled, Diff_Coupled.Benchmarks
using Plots, LinearAlgebra, LaTeXStrings
# Main function -------------------------------------------------------
function main()
    #------------------------------------------------------------------
    verbose  = false
    plot_sim = true
    global plot_end = true
    # Physics ---------------------------------------------------------
    Di      = -1.0              #Diffusion coefficient in [m^2/s]
                                #If you want to calculate D with the Arrhenius equation, set Di = [-1.0]
    D0      = 2.75*1e-6         #Pre-exponential factor in [m^2/s]
    L       = 0.001             #Length of the domain in [m]
    Cinf    = 0.0               #Concentration at infinity in [mol]
    Cstart  = 4.0               #Initial concentration in [mol]
    rho     = 1.0               #Density in [kg/m^3]
    R       = 8.314472          #Universal gas constant in [J/(mol*K)]
    Ea1     = 292880.0          #Activation energy for the left side in [J/mol]
    Myr2Sec = 60*60*24*365.25*1e6 #Conversion factor from Myr to s
    t_tot   = 1e-4 * Myr2Sec    #Total time [s]
    n       = 1                 #Geometry; 1: planar, 2: cylindric, 3: spherical
    # Numerics --------------------------------------------------------
    nx     = 100                #Number of grid points
    CFL    = 0.99               #CFL number for time step calculation
    # Domain ----------------------------------------------------------
    dx    = L*inv(nx-1)            #Grid spacing
    x     = [0:dx:L;]           #Grid points
    BCout = [1, 1]              #Boundary condition; 0: Neumann, 1: Dirichlet 
    # Initial condition -----------------------------------------------
    t     = 0.0                 #Initial time in [s]
    it    = 0                   #Time iterations
    C     = Cinf * ones(nx,1)   #Concentration array in [mol]
    C[1]  = Cstart              #Set initial concentration at the first grid point
    C0    = copy(C)             #Store initial concentration 
    x0    = copy(x)             #Store initial grid points
    #History dependent parameters--------------------------------------
    T_ar    = LinRange(1273.15,923.15,1000) #Temperature arrray in [K] to calculate temperature history; T changes with respect to time; 
                                            #The last value must be equal to the temperature at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)      #Time array (in s) to calculate history over time. The last value must be equal to t_tot.
                                            #The user is prompted to specify suitable time intervals in relation to the respective destination. 
    #Calculate values for t check---------------------------------------
    dt_diff = zeros(length(t_ar)-1)
    dt_diff = t_ar[2:end] .- t_ar[1:end-1]            
    # Preallocate variables --------------------------------------------
    Co, Co_l, Co_r, dt, dx, Kloc, Lloc, L_g, Mass, Mloc, nels, nels_l, nels_r, R_g, x_1, x_2, y_interp = preallocations(x, C, 0, 0,[0 0])
    # Calculate grid -------------------------------------------------
    dx    = L * inv(nx - 1.0)
    #Calculate initial Ds, KD, T----------------------------------------
    D, T  = update_t_dependent_param_simple!(D0,Di,dt,Ea1,R,T_ar,t_ar,t,t_tot)
    #Initial mass calculation-------------------------------------------
    Mass0 = calc_mass_vol_simple_diff(x,C,n,rho)
    #First check for correct setup--------------------------------------
    if BCout[1] != 0 && (n == 3 || n == 2)
        error("The code is only valid for cylindrical/spherical geometry, where the left outer BC has Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T  != T_ar[1]   
        error("Initial temperature must be equal to the first value in the temperature array.")
    end
    #Time loop ---------------------------------------------------------
    while t < t_tot
        #Store old values
        Co = copy(C)
        #Calculate dt --------------------------------------------------
        dt = calculate_dt(D,dx,CFL) 
        #Update time ---------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot) 
        #Update time-dependent parameters-------------------------------
        D, T = update_t_dependent_param_simple!(D0,Di,dt,Ea1,R,T_ar,t_ar,t,t_tot)
        #FEM SOLVER ----------------------------------------------------
        #Fill matrix ---------------------------------------------------
        L_g, R_g = fill_matrix!(C,x,D,dt,n,nels)
        #Reduce the condition Number -----------------------------------
        ScF      = sum(diag(L_g)) * inv(length(diag(L_g)))
        #Set boundary conditions and scale matrices --------------------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co[1],Co[end],ScF)
        #Solve system --------------------------------------------------
        C = L_g \ R_g
        if plot_sim
            # Plotting -------------------------------------------------------------
            p = plot(x,C, lw=2, label=L"Current\ concentration")
            p = plot!(x0,C0, label=L"Initial\ concentration",color=:black,linestyle=:dash,xlabel = L"Distance",
                    ylabel = L"Concentration", title = L"Simple\ diffusion\ (planar)", lw=1.5, grid=:on) 
            display(p)
        end
    end
    Mass = calc_mass_vol_simple_diff(x,C,n,rho)
    calc_mass_err(Mass,Mass0)
    return x, C, x0, C0, D, t, t_tot
end
#Run main function------------------------------------------------------
x, C, x0, C0, Di, t, t_tot  = main()
if plot_end
    # Plotting -------------------------------------------------------------
    plot(x,C, lw=2, label=L"Current\ concentration")
    plot!(x0,C0, label=L"Initial\ concentration",color=:black,linestyle=:dash,xlabel = L"Distance",
          ylabel = L"Concentration", title = L"Simple\ diffusion\ (planar)", lw=1.5, grid=:on) 
end  


