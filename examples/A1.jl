using Diff_Coupled
using Plots, LinearAlgebra, LaTeXStrings, SparseArrays
# Main function -------------------------------------------------------
function main(plot_sim)
    # Physics ---------------------------------------------------------
    Di      = 2.65*1e-18                        #Diffusion coefficient in [m^2/s]
                                                #If you want to calculate D with the Arrhenius equation, set Di = [-1.0]
    D0      = 2.75*1e-6                         #Pre-exponential factor in [m^2/s]
    L       = 0.005                             #Length of the domain in [m]
    rho     = 2700.0                            #Density in [kg/m^3]
    R       = 8.314472                          #Universal gas constant in [J/(mol*K)]
    Ea1     = 292880.0                          #Activation energy for the left side in [J/mol]
    Myr2Sec = 60*60*24*365.25*1e6               #Conversion factor from Myr to s
    t_tot   = 1e-3 * Myr2Sec                    #Total time [s]
    n       = 1                                 #Geometry; 1: planar, 2: cylindrical, 3: spherical
    # Numerics --------------------------------------------------------
    res   = 500                                 #Number of grid points
    CFL   = 0.99                                #CFL number for time step calculation
    # Domain ----------------------------------------------------------
    dx    = L*inv(res-1)                        #Grid spacing
    x     = [0:dx:L;]                           #Grid points
    BCout = [1, 1]                              #Boundary condition; 0: Neumann, 1: Dirichlet
    #Create initial profile--------------------------------------------
    Cini   = zeros(length(res))                 #Initial background concentration
    nmodes = [1; 2; 5; 7; 12] .* 1.0            #Modes of the sinusoids
    Amp    = [12; 0.5; 3; -2; 1] .* 1.0         #Initial amplitudes of the sinusoids
    C      = sinusoid_profile(Cini,nmodes,L,Di,0.0,Amp,x)
    # Initial condition -----------------------------------------------
    t      = 0.0                                #Initial time in [s]
    it     = 0                                  #Time iterations
    C0     = copy(C)                            #Store initial concentration
    x0     = copy(x)                            #Store initial grid points
    #History dependent parameters--------------------------------------
    T_ar    = LinRange(1273.15,1273.15,1000)    #Temperature array in [K] to calculate temperature history; T changes with respect to time;
                                                #The last value must be equal to the temperature at t = t_tot.
    t_ar    = LinRange(0.0,t_tot,1000)          #Time array (in [s]) to calculate history over time. The last value must be equal to t_tot.
                                                #The user is prompted to specify suitable time intervals in relation to the respective destination.
    #Calculate values for t check--------------------------------------
    dt_diff = zeros(length(t_ar)-1)
    dt_diff = t_ar[2:end] .- t_ar[1:end-1]
    #Preallocate variables --------------------------------------------
    Co      = zeros(size(C))                    #Old concentration
    dt      = 0.0                               #Initial time step
    dx      = zeros(length(x) - 1,1)            #Grid spacing
    L_g     = spzeros(length(x),length(x))      #Global matrix
    Mass    = Float64[]                         #Mass array
    nels    = length(x) - 1                     #Number of elements
    R_g     = zeros(length(x),1)                #Global vector
    #Checks------------------------------------------------------------
    MB_Error = Float64[]                        #Array to store the mass error   
    #Calculate grid ---------------------------------------------------
    dx     = L * inv(res - 1.0)
    #Calculate initial D, KD, T----------------------------------------
    D, T   = update_t_dependent_param_simple!(D0,Di,Ea1,R,T_ar,t_ar,t,t_tot)
    #Initial mass calculation------------------------------------------
    Mass0  = calc_mass_vol_simple_diff(x,C,n,rho)
    #First check for correct setup-------------------------------------
    if BCout[1] != 0 && (n == 3 || n == 2)
        error("The code is only valid for cylindrical/spherical geometry, where the left outer BC has Neumann conditions (0).")
    elseif t != 0.0
        error("Initial time must be zero.")
    elseif any(dt_diff .<= 0.0) || any(t_ar .< 0.0) || any(t_ar .> t_tot)
        error("The time array is not valid. Please check your inputs.")
    elseif T  != T_ar[1]
        error("Initial temperature must be equal to the first value in the temperature array.")
    end
    #Time loop --------------------------------------------------------
    #anim = Animation()
    while t < t_tot
        #Store old values
        Co = copy(C)
        #Calculate dt -------------------------------------------------
        dt = calculate_dt(D,dx,CFL)
        #Update time --------------------------------------------------
        t, dt, it = update_time!(t,dt,it,t_tot)
        #Update time-dependent parameters------------------------------
        D, T = update_t_dependent_param_simple!(D0,Di,Ea1,R,T_ar,t_ar,t,t_tot)
        #FEM SOLVER ---------------------------------------------------
        #Fill matrix --------------------------------------------------
        L_g, R_g = fill_matrix!(C,x,D,dt,n,nels)
        #Reduce the condition Number ----------------------------------
        ScF      = sum(diag(L_g)) * inv(length(diag(L_g)))
        #Set boundary conditions and scale matrices -------------------
        L_g, R_g = set_outer_bc!(BCout,L_g,R_g,Co[1],Co[end],ScF)
        #Solve system -------------------------------------------------
        C = L_g \ R_g
        if plot_sim && it % 15 == 0
            # Plotting ------------------------------------------------
            fs = 12.0
            p = plot(x*1000,C, lw=2, label=L"\mathrm{Current\ composition}")
            p = plot!(x0*1000,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
                    ylabel = L"C\ \mathrm{[-]}", lw=1.5, grid=:on,dpi = 300,
                        legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                        legend_foreground_color = :transparent)      
            display(p)
            #frame(anim)
        end
        # Suppress output of calc_mass_err
        redirect_stdout(devnull) do
            ErrM = calc_mass_vol_simple_diff(x,C,n,rho)
            push!(MB_Error,ErrM)
        end
    end
    #gif(anim, "figures/A1.gif", fps=10)  # Save with 10 frames per second
    Mass = calc_mass_vol_simple_diff(x,C,n,rho)
    calc_mass_err(Mass,Mass0)
    return x, C, x0, C0, D, t, t_tot, Cini, nmodes, Amp, L
end

#Run main function-----------------------------------------------------
run_and_plot = true
if run_and_plot
    plot_sim  = false
    plot_end  = true
    save_file = false
    x, C, x0, C0, Di, t, t_tot, Cini, nmodes, Amp, L  = main(plot_sim)
    xan = copy(x)
    Can = sinusoid_profile(Cini,nmodes,L,Di,t,Amp,xan)
    if plot_end
        # Title: Simple diffusion - planar
        # Plotting ----------------------------------------------------
        fs = 12.0
        plot(x*1000,C, lw=2, label=L"\mathrm{Current\ composition}")
        plot!(x0*1000,C0, label=L"\mathrm{Initial\ composition}",color=:black,linestyle=:dash,xlabel = L"x\ \mathrm{[mm]}",
              ylabel = L"C\ \mathrm{[-]}", lw=1.5, grid=:on)
        scatter!([xan[1:5:end]]*1000,[Can[1:5:end]], marker=:circle, markersize=2.0, label=L"\mathrm{Analytical\ solution}",
                    markerstrokecolor=:crimson, markercolor=:crimson)
        scatter!([xan[end]]*1000,[Can[end]], marker=:circle, markersize=2.0, label="",
                    markerstrokecolor=:crimson, markercolor=:crimson,dpi = 300,
                    legendfontsize=fs-2,guidefontsize=fs, tickfontsize=fs-1,
                    legend_foreground_color = :transparent)
        #save_path = "figures"
        #save_name = "A1"
        #save_figure(save_name,save_path,save_file)
    end
end
