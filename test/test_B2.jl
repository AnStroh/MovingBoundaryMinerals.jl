using Test
using Diff_Coupled, Diff_Coupled.Benchmarks
using LinearAlgebra, Revise, LaTeXStrings, SparseArrays
#Main function----------------------------------------------------
function main(adapt_dt,plot_sim,verbose)
    #If you find a [] with two entires this belong to the respective side of the diffusion couple ([left right])
    #Phyics-------------------------------------------------------
    Di      = [-1.0    -1.0]                                                #Initial diffusion coefficient in [m^2/s]
                                                                            #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [2.75*1e-6    3.9*1e-7;]                                      #Pre-exponential factor in [m^2/s]
    rho     = [1.0      1.0;]                                               #Normalized densities in [-]
    Ri      = [0.0001    0.0002;]                                           #Initial radii [interface    total length] in [m]
    Cl_i    = 0.6                                                           #Initial concentration left side in [mol fraction]
    Cr_i    = 0.3                                                           #Initial concentration right side in [mol fraction]
    V_ip    = 0.0                                                           #Interface velocity in [m/s]
    R       = 8.3144626182                                                  #Universal gas constant in [J/(mol*K)]
    Ea1     = 70 * 4184.0                                                   #Activation energy for the left side in [J/mol]
    Ea2     = 86.2 * 4184.0                                                 #Activation energy for the right side in [J/mol]
    dH0     = 5 * 4184.0                                                    #Standard enthalpy of reaction in [J/mol]
    Myr2Sec = 60*60*24*365.25*1e6                                           #Conversion factor from Myr to s
    t_tot   = 10.0 * Myr2Sec                                                #Total time [s]
    T0      = 1200.00                                                       #Initial maximal temperature in [K]
    s       = 10.0 * inv(Myr2Sec)                                           #Cooling rate in [K/s]
    n       = 1                                                             #Geometry; 1: planar, 2: cylindrical, 3: spherical
    #History dependent parameters---------------------------------
    t_ar    = LinRange(0.0,t_tot,10000)                                     #Time array (in s) to calculate history over time. The last value must be equal to t_tot.
                                                                            #The user is prompted to specify suitable time intervals in relation to the respective destination.
    T_ar    = T0 .* inv.( 1.0 .+ (s .* t_ar .* inv(T0)))                    #Temperature arrray in [K] to calculate temperature history; T changes with respect to time;
                                                                            #The last value must be equal to the temperature at t = t_tot.
    #Numerics-----------------------------------------------------
    CFL    = 0.10                                                           #CFL condition
    res    = [80 100;]                                                      #Number of grid points
    resmin = copy(res)                                                      #Minimum number of grid points
    MRefin = 50.0                                                           #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
    BCout  = [0 0]                                                          #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann;
                                                                            #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
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
    elseif V_ip != 0.0
        error("Please change V_ip to 0.0. This code deals with a fixed interface.")
    end
    x_left, x_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    #Preprocess and initial condition-----------------------------
    L       = Ri[end]                                                       #Length of the domain in [m]
    t       = 0.0                                                           #Initial time in [s]
    it      = 0                                                             #Initial number of time iterations
    C_left  = Cl_i*ones(res[1],1)                                           #Initial concentration left side in [mol fraction]
    C_right = Cr_i*ones(res[2],1)                                           #Initial concentration right side in [mol fraction]
    C0      = [copy(C_left); copy(C_right)]                                 #Store initial concentration
    C       = copy(C0)                                                      #Create 1 array with all concentrations
    C0_l    = copy(C_left)                                                  #Store initial concentration left side
    C0_r    = copy(C_right)                                                 #Store initial concentration right side
    x       = copy(x0)                                                      #Create 1 array containing all x-values
    Ri0     = copy(Ri)                                                      #Store initial radii
    KD0     = (Cl_i * inv(1 - Cl_i)) * inv((Cr_i * inv(1 - Cr_i)))          #Initial partition coefficient, just for pre-processing
    #Total mass---------------------------------------------------
    Mass0   = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
    #Preallocate variables----------------------------------------
    Co_l    = zeros(size(C_left))                                           #Matrix to store old concentrations of left side
    Co_r    = zeros(size(C_right))                                          #Matrix to store old concentrations of right side
    dt      = 0.0                                                           #Initial time step
    L_g     = spzeros(length(x),length(x))                                  #Global left hand side matrix
    Mass    = Float64[]                                                     #Array to store the mass of the system
    R_g     = zeros(length(x),1)                                            #Global right hand side vector
    T_pl       = []                                                         #Temperature for plotting
    t_pl       = []                                                         #Time for plotting
    Sols_left  = []                                                         #Array to store solutions for the left side
    Sols_right = []                                                         #Array to store solutions for the right side
    Checks     = []                                                         #Array to store checks
    CheckBC    = []                                                         #Array to store checks for boundary conditions
  #Calculate initial D, KD, T-----------------------------------
    KD  = copy(KD0)
    D_l = D0[1] * exp(-Ea1 * inv(R) * inv(T0))
    D_r = D0[2] * exp(-Ea2 * inv(R) * inv(T0))
    T   = T0
    #Parameter for the semi-analytical solution of Lasaga (1983)
    #doi: 10.1007/978-1-4612-5587-1-------------------------------
    betap = dH0 * s * inv(R) .* inv.(T0) .^ 2
    beta  = (dH0 * s * sqrt.(D_r .* inv.(D_l))) .* inv.((R .* T0 .^ 2 .* (sqrt.(D_r .* inv.(D_l)) .* (1 + Cl_i .* inv(1 - Cl_i)) + (Cl_i .* inv(Cr_i) + Cl_i * inv(1 - Cr_i)))))
    KD_ar = KD0 * exp.(-betap .* t_ar)                                      #Partition coefficient array to calculate partition coefficient history; KD changes with respect to time;
                                                                            #The last value must be equal to the partition coefficient at t = t_tot.
    #First check for correct setup--------------------------------
    if n != 1
        error("For this setup, you have to choose a planar geometry. Please set n = 1.")
    elseif BCout[1] != 0 && (n == 3 || n == 2)
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
        #Calculate dt---------------------------------------------
        if adapt_dt
            dt1 = dx1 .^ 2 .* inv.(D0[1] * exp(-Ea1 * inv(R) * inv(T0))) .* inv(3.0)
            dt2 = dx2 .^ 2 .* inv.(D0[2] * exp(-Ea2 * inv(R) * inv(T0))) .* inv(3.0)
            dt  = minimum([dt1 dt2]) * 1e6
        else
            dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
        end
        #Update time----------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters-------------------------
        KD  = KD0 * exp.(-betap .* t)                                       #KD value
        T   = T0 .* inv((1.0 .+ (s .* t .* inv(T0))))                       #Temperature
        D_l = D0[1] * exp(-Ea1 * inv(R) * inv(T))                           #Diffusion coefficient left
        D_r = D0[2] * exp(-Ea2 * inv(R) * inv(T))                           #Diffusion coefficient right
        #Advect interface & regrid--------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,res)
        #FEM SOLVER-----------------------------------------------
        #Construct global matrices--------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,res)
        #Set inner boundary conditions----------------------------
        L_g, R_g, ScF, BC_left, BC_right, BC_left_Las, BC_right_Las = set_inner_bc_Lasaga!(Cl_i,beta,t, KD,D_r,D_l,D0,C_left,C_right,dx1,dx2,rho,L_g,R_g,res)
        #Set outer boundary conditions and scale matrices---------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #Solve system---------------------------------------------
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Regrid---------------------------------------------------
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,verbose)
        #Post-Preprocessing---------------------------------------
        for iit in enumerate(1)
            Massnow = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            push!(Mass, Massnow)                                            #Stores the mass of the system
        end
        if mod(it,15000) == 0
            #-----------------------------------------------------
            Check1 = (C_left[end] * inv((1 - C_left[end]))) * inv((C_right[1] * inv((1 - C_right[1])))) - KD
            Check2 = D_l * inv(D_r) * (C_left[end] - C_left[end-1]) - dx1 * inv(dx2) * (C_right[2] - C_right[1])
            #-----------------------------------------------------
            T_pl       = push!(T_pl,T)                                      #Temperature for plotting
            t_pl       = push!(t_pl,t)                                      #Time for plotting
            Sols_left  = push!(Sols_left,[BC_left_Las, BC_left])            #[Semi-analytical solution, numerical solution]
            Sols_right = push!(Sols_right,[BC_right_Las, BC_right])         #[Semi-analytical solution, numerical solution]
            Checks     = push!(Checks,[C_left[end-2:end], C_right[1:3]])    #For benchmarking
            CheckBC    = push!(CheckBC,[Check1, Check2])                    #For benchmarking
        end
    end
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    #calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, C_left, C_right, C0, Sols_left, Sols_right,Checks, CheckBC, T_pl, (t_pl ./ Myr2Sec), Ri, maxC, minC
end
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (flux) - Lasaga (1983)" begin
    adapt_dt = true
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, Sols_left, Sols_right,Checks, CheckBC, T_pl, t_pl = main(adapt_dt,plot_sim,verbose)
    Can1 = first.(Sols_left)
    Can2 = first.(Sols_right)
    @test Can1 ≈ last.(Sols_left) rtol = 1e-3
    @test Can2 ≈ last.(Sols_right) rtol = 1e-3
end
