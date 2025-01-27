#= Utilized functions in the framework of the Paper:

    by A.Stroh and E. Moulas
    doi:
    Version: 1.0
    -----------------------------------------------------------------------
    The software can be used to calculate concentration profiles taking into
    account diffusion and a moving boundary. Mass balance describes
    concentration changes at the interface. 
    Phase A: left component (stable at lower T), Phase B: right component
    (stable at higher T)
==========================================================================#
using Diff_Coupled
using Plots, LinearAlgebra, DelimitedFiles, SparseArrays, LaTeXStrings
#Main function-------------------------------------------------
function main(plot_sim,verbose)
    #If you find a [] with two entires this belong to the respective side of 
    #the diffusion couple ([left right])
    #Phyics-------------------------------------------------------
    Di      = [2.47*1e-7   1*1e-7;] #Initial diffusion coefficient in [m^2/s] (Spandler & O'Neill (2009), Zhang et al. (2010))); 
                                        #If you want to calculate D with the Arrhenius equation, set Di = [-1.0 -1.0;]
    D0      = [2.75*1e-6    3.9*1e-7;]  #Pre-exponential factor in [m^2/s]
    Ri      = 0.02                      #Position of the interface -> initial radius of the left phase
    Tstart  = 1350.0 + 273.15           #Starting temperature in [K]
    Tstop   = 1350.0 + 273.15           #End temperature in [K]
    Myr2Sec = 60*60*24*365.25*1e6       #Conversion factor from Myr to s
    t_tot   = 1e-4 * Myr2Sec            #Total time [s]
    n       = 1                         #Geometry; 1: planar, 2: cylindric, 3: spherical
    CompInt    = 0.3                    #Composition of interest of the solid solution
    coeff      = readdlm("examples/Examples_phase_diagram/Coefficients_Reaction_lines.csv")
    eq_values  = [coeff[1,1]  coeff[2,1]  coeff[3,1];	    #Coefficients for composition calculation of component B (stable at higher T) X2 = a2 + b2*T + c2*T²   
                  coeff[1,2]  coeff[2,2]  coeff[3,2]]       #Coefficients for composition calculation of component A (stable at lower T) X1 = a1 + b1*T + c1*T²     
                                                            #CAUTION: The user is able to use either the variable version or Float64 values
    rho_phases = readdlm("./examples/Examples_phase_diagram/density_phases copy.tab")
    #Numerics-----------------------------------------------------
    CFL    = 0.3                    #CFL condition
    res    = [100 150;]             #Number of grid points
    resmin = copy(res)              #Minimum number of grid points
    MRefin = 15.0                     #Refinement factor
    BCout  = [0 0]                  #Outer BC at the [left right]; 1 = Dirichlet, 0 = Neumann; 
                                    #CAUTION for n = 3 the left BC must be Neumann (0)! -> right phase grows around the left phase
    #Create data set--------------------------------------------------------------
    #Create arrays X(T) using linear least squares
    coeff_up, coeff_do = coeff_trans_line(eq_values)                          #Extract coefficients
    TMAX      = Tstart + 1000.0                                             #Max temperature for T-array
    TMIN      = Tstop  - 1000.0                                             #Min temperature for T-array 
    Tlin      = LinRange(TMAX,TMIN,10000)                                       #Temperature values
    XC_left, XC_right = composition(coeff_up,coeff_do,Tlin)                 #e.g. liquidus/solidus/solvus
    C_leftlin   = copy(XC_left)                                             #Store the composition of A for later calculations
    C_rightlin  = copy(XC_right)                                            #Store the composition of B for later calculations
    #Create density plots------------------------------------------------------
    nd1       = Int(round(sqrt.(length(rho_phases[:,1]))))                                 
    Xwm       = copy(rho_phases[:,1]) 
    Xwm       = reshape(Xwm,nd1,nd1)                    #X(C1)
    Twm       = copy(rho_phases[:,2]) 
    Twm       = reshape(Twm,nd1,nd1)                    #Temperature in K
    rho_left  = copy(rho_phases[:,3]) 
    rho_left  = reshape(rho_left,nd1,nd1)               #Density component A in kg/m³
    rho_right = copy(rho_phases[:,4]) 
    rho_right = reshape(rho_right,nd1,nd1)              #Density component B in kg/m³
    #Create other arrays-------------------------------------------------------
    R_left    = C_leftlin .* inv.((1.0 .- C_leftlin))   #Concentration rate in phase A
    R_right   = C_rightlin .* inv.((1.0 .- C_rightlin)) #Concentration rate in phase B
    KDlin     = R_left .* inv.(R_right)                 #Partition coefficient KD -> see Lasaga Geospeedometry
    Tpath     = LinRange(Tstart,Tstop,10000)            #Temperature path in [K]
    tpath     = LinRange(0,t_tot+1e-10,10000)           #Time path in [s]
    #Preprocess and initial condition--------------------------------------------
    t         = 0.0                #Time
    it        = 0                  #Iterations
    D_l       = Di[1]              #Diffusion coefficient left comp.
    D_r       = Di[2]              #Diffusion coefficient right comp.
    T         = copy(Tstart)       #Initial temperature
    C_leftB, C_rightB = composition(coeff_up,coeff_do,T)                    #Initial composition of the phases
    Xc        = (1.0 - CompInt) * C_rightB + CompInt * C_leftB            #Actual total composition (green line in subplot(221)), Xc intersection with black dashed line, e.g. Xc = C(Ol)*V(Ol)+C(melt)*V(melt)
    L         = (Ri[1] ^ n * (Xc - C_leftB) * inv((C_rightB - Xc)) ^ (1 * inv(n))) + Ri[1] #Total length of the modelling domain 
    Ri        = [Ri L]                                                #Radii of the 2 phases
    #Create mesh, discretization and mapping------------------------
    if Ri[1] >= Ri[2]
        error("Please change the size of the system. Increase Ri[2].")
    end
    if res[1] > res[2]
        error("Please change the resolution of the system. res[2] >= res[1].")
    end
    x0_left, x0_right, dx1, dx2, x0 = create_grid!(Ri,res,MRefin,verbose)
    #Calculate densities-------------------------------------------------------
    rho         = calculate_density(Xwm[:,1],Twm[1,:],rho_left,rho_right,C_leftB,C_rightB,T)
    #More initial conditions---------------------------------------------------
    x_left      = copy(x0_left)                                     #Initial x_left
    x_right     = copy(x0_right)                                    #Initial x_right
    C_left      = C_leftB * ones(1,res[1])                          #Concentration of component B in phase A
    C_right     = C_rightB * ones(1,res[2])                         #Concentration of component B in phase B
    C0          = [copy(C_left) copy(C_right)]                      #Store initial concentration
    dt          = minimum([dx1,dx2]) ^ 2.0 .* inv((maximum([D_l,D_r])))  #Initial dt
    #Total mass---------------------------------------------------
    Mass0       = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)  #Initial mass
    #Preallocate variables----------------------------------------
    #Co          = zeros(size(C0))
    #Co_left     = zeros(size(C_left))
    #Co_right    = zeros(size(C_right))
    #dx          = zeros(length(x0) - 1,1)
    Kloc        = zeros(2, 2)
    Lloc        = zeros(2, 2)
    L_g         = spzeros(length(x0),length(x0))
    Mass        = Float64[]
    KD_sim      = Float64[]                    
    T_sim       = Float64[]                    
    R_left_sim  = Float64[]                    
    R_right_sim = Float64[]                    
    Mloc        = zeros(2, 2)
    nels_l      = res[1] - 1
    nels_r      = res[2] - 1
    R_g         = zeros(length(x0),1)
    #---------------------------------------------------------------------------
    #Solving the moving boundary problem---------------------------------------
    #while t < t_tot
    for iit = 1:1000
    println("time: $t")
        #Update time---------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        if t <= t_tot          
            T  = linear_interpolation_1D(tpath,Tpath,t)
        end
        #Calculate Equilibrium compositions at actual T------------------------
        C_left[end], C_right[1] = composition(coeff_up,coeff_do,T)
        dC  = C_right[1] - C_left[end]          #Composition difference
        rho = calculate_density(Xwm[:,1],Twm[1,:],rho_left,rho_right,C_leftB,C_rightB,T)
        #Stefan condition -> Composition difference----------------------------                     
        JL   = - D_l * rho[1] * (C_left[end] - C_left[end-1]) * inv(dx1)        #Flux of the left side to the right side
        JR   = - D_r * rho[2] * (C_right[2]  - C_right[1])    * inv(dx2)        #Flux of the right side to the left side
        V_ip = (JR - JL) * inv(dC)                                              #Velocity in x direction
        #Advect interface & regrid--------------------------------------
        Fl_regrid, x_left, x_right, C_left, C_right, res, Ri = advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,vec(C_left),vec(C_right),res)
        #FEM SOLVER-----------------------------------------------------
        #Construct global matrices--------------------------------------
        L_g, R_g, Co_l, Co_r = construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,nels_l,nels_r,Mloc,Kloc,Lloc,res)
        @show det(L_g)
        #if det(L_g) == 0.0
        #    error("Check1: Matrix L_g is singular")
        #end
        #Set inner boundary conditions----------------------------------
        L_g, R_g, ScF = set_inner_bc_stefan!(L_g,R_g,C_left,C_right,res)
        #if det(L_g) == 0.0
        #    error("Check2: Matrix L_g is singular ")
        #end
        #Set outer boundary conditions and scale matrices---------------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co_l[1],Co_r[end],ScF)
        #if det(L_g) == 0.0
        #    error("Check3: Matrix L_g is singular ")
        #end
        #Solve system---------------------------------------------------
        println("Temperature: $T")
        C_left, C_right = solve_soe(L_g,R_g,res)
        #Regrid---------------------------------------------------------
        x_left, x_right, C_left, C_right, dx1, dx2, res = regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, res, resmin, MRefin,verbose)
        #Post-Preprocessing---------------------------------------------
        for iit in enumerate(1)       
            Massnow = calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
            R_left     = C_left[end] .* inv((1.0 - C_left[end]))
            R_right    = C_right[1] .* inv((1.0 - C_right[1]))
            KDt        = R_left .* inv(R_right)
            T_pl       = copy(T)
            push!(Mass, Massnow)  #Stores the mass of the system
            push!(KD_sim, KDt)    #Stores the partition coefficient
            push!(T_sim, T_pl)     #Stores the temperature
            push!(R_left_sim, R_left)   #Concentration ratio left side (interface)
            push!(R_right_sim, R_right) #Concentration ratio right side (interface)
        end
        #Find the important dt-------------------------------------------------
        dtV   = minimum([dx1,dx2]) ^1 .* inv((abs(V_ip)))
        dtD   = minimum([dx1,dx2]) ^2 .* inv((maximum([D_l,D_r])))
        dt    = minimum([dtD,dtV]) * CFL
        if plot_sim
            #Plotting------------------------------------------------------
            maxC = maximum([maximum(C_left),maximum(C_right)])
            Tp_min = Tstop * 0.95
            Tp_max = Tstart * 1.05
            first_val, last_val = values_between_known_indices!(Tlin,KDlin,Tstart,Tstop)        #Just works for constantly dropping temperature
            #Concentration profile
            p1 = plot(x_left,C_left, lw=2, label=L"Left\ side")
            p1 = plot!(x_right,C_right, lw=2, label=L"Right\ side")
            p1 = plot!([x_left[end]; x_left[end]],[0; 1]*maxC,color=:grey,linestyle=:dash, lw=2, label=L"Interface")
            p1 = plot!(x0,C0',color=:black,linestyle=:dash,xlabel = L"Distance", ylabel = L"Concentration", title = L"Concentration\ profile", lw=1.5,
                  grid=:on, label=L"Initial\ condition")
            p1 = annotate!(0.05, maxC - 0.05,"(d)", transform=:axes)
            #Phase diagram
            p2 = plot(Tlin,XC_left, lw=2, label=L"Left\ side")
            p2 = plot!(Tlin,XC_right, lw=2, label=L"Right\ side")
            p2 = scatter!([T],[C_left[end]],marker=:circle, markersize=2, markercolor=:black, 
                          markerstrokecolor=:black,label = "")
            p2 = scatter!([T],[C_right[1]],marker=:circle, markersize=2, markercolor=:black, 
                          markerstrokecolor=:black,label = "")
            p2 = plot!([T; T],[0; maximum([C_left[end],C_right[1]])],lw=1.5, label="",color=:black,linestyle=:dash)
            p2 = plot!([T; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:midnightblue,linestyle =:dot)
            p2 = plot!([T; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), ylims=(0, 1),color=:crimson,linestyle =:dot)
            p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass[end]; Mass[end]],color=:purple,linestyle=:dash,lw=1.5, label=L"Final\ mass")
            p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass0; Mass0],color=:green,linestyle=:dash,lw=1.5, label=L"Initial\ mass",
                        xlabel = L"Temperature", ylabel = L"X", title = L"Phase\ diagram",grid=:on)
            p2 = annotate!(0.05, 0.95,"(a)", transform=:axes)
            #evolution of KD(T)
            p3 = plot(Tlin, KDlin, lw=2, label="")
            p3 = scatter!([T_sim[end]],[KD_sim[end]],marker=:circle, markersize=2, markercolor=:black,markerstrokecolor=:black,
                          xlabel = L"Temperature", ylabel = L"K_{D}", title = L"K_{D}(T)\ evolution", lw=1.5,
                          grid=:on, label="",xlims=(Tp_min, Tp_max), ylims=(first_val-0.01,last_val+0.01))
            p3 = annotate!(0.05, maximum(KDlin)- 0.05, "(b)", transform=:axes)
            #ln(KD) vs 1/T)
            #p4 = plot(1.0 ./ T_sim,log.(KD_sim),xlabel = L"1/T", ylabel = L"ln(K_{D})", title = L"Arrhenius plot", lw=1.5,
            #            grid=:on, label="")
            #p4 = annotate!(0.05, maximum(log.(KD_sim)) - 0.05,"(c)", transform=:axes)
            p = plot(p2,p3,p1,suptitle = L"Thermodynamical\ constrained,\ Stefan\ condition")
            display(p)
        end
    end
    maxC = maximum([maximum(C_left),maximum(C_right)])
    minC = minimum([minimum(C_left),minimum(C_right)])
    ErrM = (Mass[end] - Mass0) / Mass0
    println("The total mass error is $(ErrM[end]*100)% using fluxes.")    
    C_left_PD, C_right_PD = composition(coeff_up,coeff_do,T)
    ErrCR = (C_right[1] - C_right_PD) / C_right_PD
    ErrCL = (C_left[end] - C_left_PD) / C_left_PD
    println("The concentration error for the left phase is $(ErrCL*100)%.")
    println("The concentration error for the right phase is $(ErrCR*100)%.")
    return x_left, x_right, x0, vec(C_left), vec(C_right), vec(C0),maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, KDlin, KD_sim,T_sim, Mass0, Mass 
end
#Calculate the Stefan condition--------------------------------
run_and_plot = true
if run_and_plot
    plot_sim = true
    plot_end = true
    verbose  = false
    x_left, x_right, x0, C_left, C_right, C0, maxC, Tlin, XC_left, XC_right, T, Tstart, Tstop, KDlin, KD_sim,T_sim, Mass0, Mass = main(plot_sim,verbose)
    if plot_end
        #Plotting------------------------------------------------------
        Tp_min = Tstop * 0.95
        Tp_max = Tstart * 1.05
        first_val, last_val = values_between_known_indices!(Tlin,KDlin,Tstart,Tstop)        #Just works for constantly dropping temperature
        #Concentration profile
        p1 = plot(x_left,C_left, lw=2, label=L"Left\ side")
        p1 = plot!(x_right,C_right, lw=2, label=L"Right\ side")
        p1 = plot!([x_left[end]; x_left[end]],[0; 1]*maxC,color=:grey,linestyle=:dash, lw=2, label=L"Interface")
        p1 = plot!(x0,C0,color=:black,linestyle=:dash,xlabel = L"Distance", ylabel = L"Concentration", title = L"Concentration\ profile", lw=1.5,
              grid=:on, label=L"Initial\ condition")
        p1 = annotate!(0.05, maxC - 0.05,"(d)", transform=:axes)
        #Phase diagram
        p2 = plot(Tlin,XC_left, lw=2, label=L"Left\ side")
        p2 = plot!(Tlin,XC_right, lw=2, label=L"Right\ side")
        p2 = scatter!([T],[C_left[end]],marker=:circle, markersize=2, markercolor=:black, 
                      markerstrokecolor=:black,label = "")
        p2 = scatter!([T],[C_right[1]],marker=:circle, markersize=2, markercolor=:black, 
                      markerstrokecolor=:black,label = "")
        p2 = plot!([T; T],[0; maximum([C_left[end],C_right[1]])],lw=1.5, label="",color=:black,linestyle=:dash)
        p2 = plot!([T; 0],[C_left[end];C_left[end]],lw=1.5, label="",color=:midnightblue)
        p2 = plot!([T; 0],[C_right[1];C_right[1]],lw=1.5, label="",xlims=(Tp_min, Tp_max), ylims=(0, 1),color=:crimson)
        p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass[end]; Mass[end]],color=:purple,linestyle=:dash,lw=1.5, label="")
        p2 = plot!([Tstop*0.3; Tstart*1.5],[Mass0; Mass0],color=:green,linestyle=:dash,lw=1.5, label="",
                    xlabel = L"Temperature", ylabel = L"X", title = L"Phase\ diagram",grid=:on)
        p2 = annotate!(0.05, 0.95,"(a)", transform=:axes)
        #evolution of KD(T)
        p3 = plot(Tlin, KDlin, lw=2, label="")
        p3 = scatter!([T_sim[end]],[KD_sim[end]],marker=:circle, markersize=2, markercolor=:black,markerstrokecolor=:black,
                      xlabel = L"Temperature", ylabel = L"K_{D}", title = L"K_{D}(T)\ evolution", lw=1.5,
                      grid=:on, label="",xlims=(Tp_min, Tp_max), ylims=(first_val-0.01,last_val+0.01))
        p3 = annotate!(0.05, maximum(KDlin)- 0.05, "(b)", transform=:axes)
        #ln(KD) vs 1/T)
        #p4 = plot(1.0 ./ T_sim,log.(KD_sim),xlabel = L"1/T", ylabel = L"ln(K_{D})", title = L"Arrhenius plot", lw=1.5,
        #            grid=:on, label="")
        #p4 = annotate!(0.05, maximum(log.(KD_sim)) - 0.05,"(c)", transform=:axes)
        plot(p2,p3,p1,suptitle = L"Thermodynamical\ constrained,\ Stefan\ condition")
    end
end    