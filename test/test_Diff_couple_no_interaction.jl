using Test
using MOBILE
using LinearAlgebra, LaTeXStrings,SparseArrays
#Main function----------------------------------------------------
function main(plot_sim,verbose)
    #If you find a [] with two entries this belong to the respective side of the diffusion couple ([left right])
    #Physics-------------------------------------------------------
    Di      = [2.65*1e-18   2.65*1e-18;]                #Initial diffusion coefficient in [m^2/s];
                                                        #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [NaN          NaN;]                       #Pre-exponential factor in [m^2/s]
    rho     = [1.0          1.0;]                       #Normalized densities in [-]
    Ri      = [0.0002       0.0004;]                    #Initial radii [interface    total length] in [m]
    Cl_i    = 0.6                                       #Initial composition left side in [-]
    Cr_i    = 0.6                                       #Initial composition right side in [-]
    V_ip    = 0.0                                       #Interface velocity in [m/s]; CAUTION: This code only gives correct results, if V_ip is 0.
    R       = NaN                                       #Universal gas constant in [J/(mol*K)]
    Ea1     = NaN                                       #Activation energy for the left side in [J/mol]
    Ea2     = NaN                                       #Activation energy for the right side in [J/mol]
    Myr2Sec = 60*60*24*365.25*1e6                       #Conversion factor from Myr to s
    t_tot   = 1e-5 * Myr2Sec                            #Total time [s]
    n       = 1                                         #Geometry; 1: planar, 2: cylindrical, 3: spherical
    #History dependent parameters---------------------------------
    KD_ar   = LinRange(1.0,1.0,1000)                    #Partition coefficient array to calculate distribution coefficient history; KD changes with respect to time;
                                                        #The last value must be equal to the distribution coefficient at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)                  #Time array (in [s]) to calculate history over time. The last value must be equal to t_tot.
                                                        #The user is prompted to specify suitable time intervals in relation to the respective destination.
    T_ar    = LinRange(1273.15,1273.15,1000)            #Temperature array in [K] to calculate temperature history; T changes with respect to time;
                                                        #The last value must be equal to the temperature at t = t_tot.
    #Numerics-----------------------------------------------------
    CFL    = 0.5                                        #CFL condition
    res    = [500 500;]                                 #Number of nodes
    resmin = copy(res)                                  #Minimum number of nodes
    MRefin = 1.0                                        #Refinement factor; If negative, it uses MRefin = 1 on the left, and abs(MRefin) on the right
    Dirichlet = 1                                       #Dirchlet BC at the interface
    BCout  = [1 1]                                      #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann;
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
        error("Please change V_ip to 0.0. This code cannot handle moving interface.")
    elseif any(KD_ar .!= 1.0)
        error("Please change each entry of KD_ar to 1.0.")
    end
    x_left, x_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    #Create initial profile---------------------------------------
    Cini_l   = zeros(length(x_left))                    #Initial background composition left
    nmodes_l = [1; 2; 5; 7; 12] .* 1.0                  #Modes of the sinusoids left
    Amp_l    = [12; 0.5; 3; -2; 1] .* 1.0               #Initial amplitudes of the sinusoids left
    C_left   = sinusoid_profile(Cini_l,nmodes_l,Ri[1],Di[1],0.0,Amp_l,x_left)
    Cini_r   = zeros(length(x_right))                   #Initial background composition right
    nmodes_r = [1; 2; 3] .* 1.0                         #Modes of the sinusoids right
    Amp_r    = [1; 8.0; 1*inv(3.0)] .* 1.0              #Initial amplitudes of the sinusoids right
    C_right  = sinusoid_profile(Cini_r,nmodes_r,Ri[2]-Ri[1],Di[2],0.0,Amp_r,x_right)
    #Preprocess and initial condition-----------------------------
    L       = Ri[end]                                   #Length of the domain in [m]
    t       = 0.0                                       #Initial time in [s]
    it      = 0                                         #Initial number of time iterations
    C0      = [copy(C_left); copy(C_right)]             #Store initial composition
    C       = copy(C0)                                  #Create 1 array with all concentrations
    C0_l    = copy(C_left)                              #Store initial composition left side
    C0_r    = copy(C_right)                             #Store initial composition right side
    x       = copy(x0)                                  #Create 1 array containing all x-values
    Ri0     = copy(Ri)                                  #Store initial radii
    KD      = copy(KD_ar[1])                            #Initial distribution coefficient, just for pre-processing
    #Total mass---------------------------------------------------
    Mass0   = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
    #Preallocate variables----------------------------------------
    Co_l    = zeros(size(C_left))                       #Matrix to store old concentrations of left side
    Co_r    = zeros(size(C_right))                      #Matrix to store old concentrations of right side
    dt      = 0.0                                       #Initial time step
    L_g     = spzeros(length(x),length(x))              #Global left hand side matrix
    Mass    = Float64[]                                 #Array to store the mass of the system
    R_g     = zeros(length(x),1)                        #Global right hand side vector
    #Checks-------------------------------------------------------
    MB_Error = Float64[]                                #Array to store the mass error
    #Calculate initial Ds, KD, T----------------------------------
    D_l, D_r, KD, T = update_t_dependent_param!(D0,Di,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
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
    while t < t_tot
        #Calculate dt---------------------------------------------
        dt = find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
        #Update time----------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters-------------------------
        D_l, D_r, KD,T = update_t_dependent_param!(D0,Di,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
        #FEM SOLVER-----------------------------------------------
        #Construct global matrices--------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,res)
        #Set inner boundary conditions----------------------------
        if Dirichlet == 1
            ScF = 1.0
            #fill!(L_g[res[1],:],0.0)
            L_g[res[1],:]           .=  0.0
            L_g[res[1],res[1]]       =  1.0 * ScF
            R_g[res[1]]              =  C_left[end] * ScF
            #fill!(L_g[res[1]+1,:],0.0)
            L_g[res[1]+1,:]         .=  0.0
            L_g[res[1]+1,res[1]+1]   =  1.0 * ScF
            R_g[res[1]+1]            =  C_right[1] * ScF
        end
        #Set outer boundary conditions and scale matrices---------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #Solve system---------------------------------------------
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Post-Preprocessing---------------------------------------
        for iit in enumerate(1)
            Massnow = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            push!(Mass, Massnow)  #Stores the mass of the system
        end
        # Suppress output of calc_mass_err
        redirect_stdout(devnull) do
            ErrM = calc_mass_err(Mass, Mass0)
            push!(MB_Error,ErrM)
        end
    end
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    calc_mass_err(Mass,Mass0)
    return x_left, x_right, x0, C_left, C_right, C0, Cini_l, Cini_r, nmodes_l, nmodes_r,Amp_l, Amp_r, Di, Ri, t
end
#Testing-----------------------------------------------------------------------
@testset "Diffusion couple (no interaction)" begin
    plot_sim = false
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, Cini_l, Cini_r, nmodes_l, nmodes_r,Amp_l, Amp_r, Di, Ri, t = main(plot_sim,verbose)
    xan_l = copy(x_left)
    xan_r = copy(x_right)
    Can_l = sinusoid_profile(Cini_l,nmodes_l,Ri[1],Di[1],t,Amp_l,xan_l)
    Can_r = sinusoid_profile(Cini_r,nmodes_r,Ri[2]-Ri[1],Di[2],t,Amp_r,x_right)
    @test Can_l ≈ C_left rtol = 1e-5
    @test Can_r ≈ C_right rtol = 1e-4
end
