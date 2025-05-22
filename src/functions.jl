using SparseArrays, LinearAlgebra, BenchmarkTools, Revise, Dates, Plots
export advect_interface_regrid!, blkdiag, calculate_dt, calc_mass_vol, calc_mass_vol_simple_diff, calc_volume, create_grid!, find_dt, fill_matrix!, linear_interpolation_1D, linspace_interface, preallocations, regrid!, set_inner_bc_mb!, set_inner_bc_flux!,set_inner_bc_Lasaga!, set_outer_bc!, trapezoidal_integration, update_time!, update_t_dependent_param!, update_t_dependent_param_simple!, construct_matrix_fem, solve_soe,calc_mass_err, make_dx_right, newton_solver, define_new_grid, sinusoid_profile,save_figure,scaling,rescale
#Functions----------------------------------------------------
"""
    advect_interface_regrid!(Ri, V_ip, dt, x_left, x_right, C_left, C_right, nr)

Update the interface position and calculate new grids based on the advection velocity. Units may differ from SI units if
non-dimensionalisation has been performed.

# Arguments
- `Ri::Float64`: Radii [interface    total length in [m].
- `V_ip::Float64`: Advection velocity in [m/s].
- `dt::Float64`: Time step in [s].
- `x_left::Vector{Float64}`: Left distance nodes in [m].
- `x_right::Vector{Float64}`: Right distance nodes in [m].
- `C_left::Vector{Float64}`: Composition values on the left nodes in [-].
- `C_right::Vector{Float64}`: Composition values on the right nodes in [-].
- `nr::Vector{Int}`: Resolution of the left and the right grid.

# Returns
- `Fl_regrid::Int`: Flag indicating if regridding was performed (1) or not (0).
- `x_left::Vector{Float64}`: Updated left distance nodes.
- `x_right::Vector{Float64}`: Updated right distance nodes.
- `C_left::Vector{Float64}`: Updated composition values on the left nodes.
- `C_right::Vector{Float64}`: Updated composition values on the right nodes.
- `nr::Vector{Int}`: Updated resolution.
- `Ri::Float64`: Updated radii [interface    total length].
"""
function advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,nr)
    Rio   = copy(Ri)
    Ri[1] = Rio[1] + V_ip * dt                                                                      #Update interface position
    if V_ip > 0                                                                                     #Calculate new grid for positive velocity
        x_left      = [x_left; x_left[end]]                                                         #Add new grid point
        C_left      = [C_left; C_left[end]]                                                         #Add new composition value
        x_left[end] = copy(Ri[1])                                                                   #Update grid point
        dx1         = x_left[end] - x_left[end-1]                                                   #Calculate new dx on the left side of interface
        if Ri[1] > x_right[2]                                                                       #Check if Ri moved to fast
            @show x_right[1] x_right[2] Ri[1]
            error("Interface moved too fast. Ri is larger than x_right[2].")
        end
        x_right[1]  = copy(Ri[1])                                                                   #Update grid point
        dx2         = x_right[2] - x_right[1]                                                       #Calculate new dx on the right side of interface
        Fl_regrid   = 1                                                                             #Flag for regridding
        nr[1]       = nr[1] + 1                                                                     #Update resolution
    elseif V_ip < 0                                                                                 #Calculate new grid for negative velocity
        x_right     = [x_right[1];x_right]                                                          #Add new grid point
        C_right     = [C_right[1];C_right]                                                          #Add new composition value
        x_right[1]  = copy(Ri[1])                                                                   #Update grid point
        dx2         = x_right[2] - x_right[1]                                                       #Calculate new dx on the right side of interface
        if Ri[1] < x_left[end-1]                                                                    #Check if Ri moved to fast
            @show x_left[end-1] x_left[end] Ri[1]
            error("Interface moved too fast. Ri is smaller than x_left[end-1].")
        end
        x_left[end] = copy(Ri[1])                                                                   #Update grid point
        dx1         = x_left[end] - x_left[end-1]                                                   #Calculate new dx on the left side of interface
        Fl_regrid   = 1                                                                             #Flag for regridding
        nr[2]       = nr[2] + 1                                                                     #Update resolution
    else                                                                                            #Leave grid as it was, if velocity = 0.0
        Fl_regrid   = 0                                                                             #Flag for regridding
        dx1         = x_left[end] - x_left[end-1]                                                   #Calculate new dx on the left side of interface
        dx2         = x_right[2] - x_right[1]                                                       #Calculate new dx on the right side of interface
    end
    return Fl_regrid, x_left, x_right, C_left, C_right, nr, Ri
end


"""
    blkdiag(matrices::AbstractMatrix...)

Constructs a block diagonal matrix from a variable number of input matrices.

# Arguments
- `matrices::AbstractMatrix...`: Input matrices to be placed on the diagonal.

# Returns
- `result::Matrix`: Block diagonal matrix constructed from the input matrices.
"""

function blkdiag(matrices::AbstractMatrix...)
    # Calculate total size for the block diagonal matrix
    total_rows = sum(size(m, 1) for m in matrices)
    total_cols = sum(size(m, 2) for m in matrices)
    # Create an empty matrix with the total size filled with zeros
    result = zeros(eltype(matrices[1]), total_rows, total_cols)
    # Place each matrix on the diagonal
    row_offset = 0
    col_offset = 0
    for mat in matrices
        rows, cols = size(mat)
        result[row_offset+1:row_offset+rows, col_offset+1:col_offset+cols] .= mat
        row_offset += rows
        col_offset += cols
    end
    return result
end

"""
    blocktest(L1, R1, L2, R2)

Constructs a block matrix and a block vector from given input matrices and vectors.

# Arguments
- `L1::Matrix`: The first input matrix.
- `R1::Vector`: The first input vector.
- `L2::Matrix`: The second input matrix.
- `R2::Vector`: The second input vector.

# Returns
- `Lblock::SparseMatrixCSC`: The block matrix constructed from `L1` and `L2`.
- `Rblock::Vector`: The block vector constructed from `R1` and `R2`.
"""
function blocktest(L1,R1,L2,R2)
    n1 = length(R1)                                                                                 #Length of R1
    n2 = length(R2)                                                                                 #Length of R2
    Lblock = spzeros(n1+n2,n1+n2)                                                                   #Preallocate Lblock
    Rblock = zeros(n1+n2)                                                                           #Preallocate Rblock
    #Set block matrix and vector------------------------------
    for i in 1:n1                                                                                   #Set block matrix part L1 and vector part R1
       if i == 1                                                                                    #First row
        Lblock[i,i] = L1[1,1]
        Lblock[i,2] = L1[1,2]
       elseif i == n1                                                                               #Last row
        Lblock[i,i]   = L1[n1,n1]
        Lblock[i,i-1] = L1[n1,n1-1]
       else                                                                                         #Middle rows
        Lblock[i,i-1] = L1[i,i-1]
        Lblock[i,i+0] = L1[i,i+0]
        Lblock[i,i+1] = L1[i,i+1]
       end
       Rblock[i] = R1[i]
    end
    for i in 1:n2                                                                                   #Set block matrix part L2 and vector part R2
        if i == 1                                                                                   #First row
            Lblock[i+n1,i+n1] = L2[1,1]
            Lblock[i+n1,2+n1] = L2[1,2]
           elseif i == n2                                                                           #Last row
            Lblock[i+n1,i+n1]   = L2[n2,n2]
            Lblock[i+n1,i-1+n1] = L2[n2,n2-1]
           else                                                                                     #Middle rows
            Lblock[i+n1,i-1+n1] = L2[i,i-1]
            Lblock[i+n1,i+0+n1] = L2[i,i+0]
            Lblock[i+n1,i+1+n1] = L2[i,i+1]
        end
        Rblock[i+n1] = R2[i]
    end
    return Lblock, Rblock
end

"""
    calculate_dt(D, dx, CFL)

Calculate the time step `dt` for a diffusion process. Units may differ from SI units if non-dimensionalisation
has been performed.

# Arguments
- `D::Float64`: The diffusion coefficient in [m^2/s].
- `dx::Float64`: The spatial step size in [m].
- `CFL::Float64`: The Courant-Friedrichs-Lewy (CFL) number.

# Returns
- `dt::Float64`: The calculated time step.
"""

function calculate_dt(D,dx,CFL)
    #Calculate dt---------------------------------------------
    dt = CFL * (dx .^ 2) *inv(D)
    return dt
end

"""
    calc_mass_err(Mass, Mass0)

Calculate the mass error between the final mass `Mass[end]` and the initial mass `Mass0`.

# Arguments
- `Mass::Vector`: A vector containing the mass values in [-].
- `Mass0::Number`: The initial mass value in [-].

# Output
- `ErrM::Number`: The calculated mass error.

"""
function calc_mass_err(Mass,Mass0)
    ErrM = (Mass[end] - Mass0) * inv(Mass0)
    println("The total mass difference is $(ErrM[end]*100)% (relevant for closed systems).")
    return ErrM
end

"""
    calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)

Calculate the total mass based on the volume of the phase.

# Arguments
- `x_left::Float64`: Left nodes in [m].
- `x_right::Float64`: Right nodes in [m].
- `C_left::Vector{Float64}`: Composition values of the left phase in [-].
- `C_right::Vector{Float64}`: Composition values of the right phase in [-].
- `n::Int`: Number which defines the geometry.
- `rho::Vector{Float64}`: Densities of the left and right phase [kg/m³].

# Returns
- `Mtot::Float64`: The total mass.
"""
function calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
    #Calculate total mass based on the volume-----------------
    V_left_ini, V_right_ini, dVC_ini = calc_volume(x_left,x_right,n)                                #Extract volumes
    M_left  = trapezoidal_integration(rho[1]*V_left_ini,C_left)                                     #Calculate mass left phase
    M_right = trapezoidal_integration(rho[2]*V_right_ini,C_right)                                   #Calculate mass right phase
    Mtot    = M_left + M_right                                                                      #Total mass
    return Mtot
end

"""
    calc_mass_vol_simple_diff(x1, C1, ndim, rho)

Calculate the total mass based on the volume of a phase. This function is written for a calculation of one phase.

# Arguments
- `x1::Vector{Float64}`: A vector containing the spatial coordinates.
- `C1::Vector{Float64}`: A vector containing the composition values at the spatial coordinates.
- `ndim::Int`: The number related to the geometry.
- `rho::Vector{Float64}`: A vector containing the density values.

# Returns
- `Mtot::Float64`: The total mass calculated using trapezoidal integration.

# Description
This function calculates the total mass based on the volume for a simple diffusion process. It first preallocates arrays
for volume (`V`) and volume change (`dV`). It then calculates the volumes for each spatial coordinate and the volume changes
between consecutive coordinates. The total volume change (`dVC`) is computed by averaging the volume changes. Finally, the
total mass (`Mtot`) is calculated using trapezoidal integration of the product of density and volume with respect to the
composition.

"""

function calc_mass_vol_simple_diff(x1,C1,ndim,rho)
    #Calculate total mass based on the volume-----------------
    #Preallocations
    V = zeros(length(x1),1)
    dV = zeros(length(x1),1)
    #Calculate volumes
    for (i,_) in enumerate(1:length(x1))
        V[i]   = x1[i] ^ ndim			                                                            #Volume
    end
    for (i,_) in enumerate(1:length(V)-1)
        dV[i] = (V[i+1] - V[i])                                                                     #Volume change
    end
    dV1 = [dV; 0] * 0.5
    dV2 = [0; dV] * 0.5
    dVC    = dV1 + dV2                                                                              #Total volume change
    Mtot  = trapezoidal_integration(rho[1]*V,C1)                                                    #Total mass calculation
    return Mtot
end

"""
    calc_volume(x1, x2, ndim)

Calculation of all volumes. The density in both phases is constant. Subsequently, shrinking and expanding volumes are not
considered.

# Arguments
- `x1`: Nodes of the left phase.
- `x2`: Nodes of the right phase.
- `ndim`: Geometry factor.

# Returns
- `V1`: Array of volumes for the left phase.
- `V2`: Array of volumes for the right phase.
- `dVC`: Array of total volume changes.

"""
function calc_volume(x1,x2,ndim)
    #=Calculation of all volumes. We assume a constant density in both phases,
    which is why the density does not need to be included in the calculations.
    Shrinking or expanding volumes are not considered.
    ==========================================================
    =#
    #Preallocations
    V1 = zeros(length(x1),1)
    V2 = zeros(length(x2),1)
    dV = zeros(length(x1)+length(x2)-1,1)
    #Calculate volumes
    for (i,_) in enumerate(1:length(x1))
        V1[i]   = x1[i] ^ ndim			                                                            #Volume left phase
    end
    for (i,_) in enumerate(1:length(x2))
        V2[i]   = x2[i] ^ ndim			                                                            #Volume right phase
    end
    V = [V1; V2]                                                                                    #Total volume
    for (i,_) in enumerate(1:length(V)-1)
        dV[i] = (V[i+1] - V[i])                                                                     #Volume change
    end
    dV1 = [dV; 0] * 0.5
    dV2 = [0; dV] * 0.5
    dVC = dV1 + dV2                                                                                 #Total volume change
    return V1, V2, dVC
end

"""
    construct_matrix_fem(x_left, x_right, C_left, C_right, D_l, D_r, dt, n, res)

Constructs the global matrix for the FEM solver in a diffusion couple problem. Units may differ from SI units if
non-dimensionalisation has been performed.

# Arguments
- `x_left::Vector{Float64}`: Left grid spatial points in [m].
- `x_right::Vector{Float64}`: Right grid spatial points in [m].
- `C_left::Vector{Float64}`: Composition values of the left phase in [-].
- `C_right::Vector{Float64}`: Composition values of the right phase in [-].
- `D_l::Float64`: Diffusion coefficient on the left side in [m²/s].
- `D_r::Float64`: Diffusion coefficient on the right side in [m²/s].
- `dt::Float64`: Time step in [s].
- `n::Int`: Geometry definition.
- `res::Vector{Float64}`: Resolution.

# Returns
- `L_g::SparseMatrixCSC{Float64, Int}`: Global stiffness matrix.
- `R_g::Vector{Float64}`: Global RHS vector.
- `Co_l::Vector{Float64}`: Stores left side composition values before the update.
- `Co_r::Vector{Float64}`: Stores right side composition values before the update.
"""
function construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,res)
    #Calculate number of elements-----------------------------
    nels_l      = res[1] - 1
    nels_r      = res[2] - 1
    #Store old values for RHS---------------------------------
    Co_l = copy(C_left)
    Co_r = copy(C_right)
    #Left side matrix-----------------------------------------
    #Fill matrix----------------------------------------------
    L1_g, R1_g = fill_matrix!(C_left,x_left,D_l,dt,n,nels_l)
    #Right side matrix----------------------------------------
    #Fill matrix----------------------------------------------
    L2_g, R2_g = fill_matrix!(C_right,x_right,D_r,dt,n,nels_r)
    #Merge systems--------------------------------------------
    L_g = blkdiag(L1_g,L2_g)
    R_g = [R1_g; R2_g]
    #L_g2,R_g2 = blocktest(L1_g,R1_g,L2_g,R2_g)
    #check_block = maximum(abs.(L_g2[:] .- L_g[:]));
    #println("Error on blk ",check_block)
    return L_g, R_g, Co_l, Co_r
end

"""
    create_grid!(Ri, nr, MRefin,verbose)

Create grid with or without variable spacing.Units may differ from SI units if non-dimensionalisation has been performed.

# Arguments
- `Ri::Vector{Float64}`: Initial radii [interface    total length] in [m].
- `nr::Vector{Int}`: Resolution vector.
- `MRefin::Int`: The refinement factor for the grid.
- `verbose::Bool`: A boolean indicating whether to print additional information.

# Returns
- `x_left::Matrix{Float64}`: Left nodes in [m].
- `x_right::Matrix{Float64}`: Right nodes in [m].
- `dx1::Float64`: The grid spacing on the left side in [m].
- `dx2::Float64`: The grid spacing on the right side in [m].
- `x0::Matrix{Float64}`: Initial grid spacing for the whole domain in [m].

"""
function create_grid!(Ri,nr,MRefin,verbose)
    Ri, nr, x_left, x_right, dx_left, dx_right = define_new_grid(Ri,nr,MRefin,verbose)              #Define new grid
    dx1 = x_left[end] - x_left[end-1]                                                               #Last dx on the left side
    dx2 = x_right[2]  - x_right[1]                                                                  #First dx on the right side
    x0  = [copy(x_left); copy(x_right)]                                                             #Initial grid spacing
    return x_left, x_right, dx1, dx2, x0
end

"""
    define_new_grid(Ri, nr, Rfact, verbose)

This function defines a new grid based on the given parameters. Units may differ from SI units if non-dimensionalisation
has been performed.

## Arguments
- `Ri`: Radii [interface    total length] in [m].
- `nr`: Resolution of nodes on the left and right side.
- `Rfact`: Grid refinement factor.
- `verbose`: A boolean indicating whether to print additional information.

## Returns
- `Ri`: Radii [interface    total length] in [m].
- `nr`: Resolution of nodes on the left and right sides.
- `x_left`:Left nodes in [m].
- `x_right`: Right nodes in [m].
- `dx_left`: Left grid spacing in [m].
- `dx_right`: Right grid spacing in [m].
- `Sc_left`: The scaling factor for the left side.
- `Sc_right`: The scaling factor for the right side.

"""
function define_new_grid(Ri,nr,Rfact,verbose)
    if Rfact == 1.0                                                                                 #Equally spaced grid
        x_left   = collect(LinRange(0.0, Ri[1], nr[1]))                                             #Create vector with nr[1] equally spaced nodes
        x_right  = collect(LinRange(Ri[1],Ri[2], nr[2]))                                            #Create vector with nr[2] equally spaced nodes
        dx_left  = diff(x_left)                                                                     #Calculate dx on the left side
        dx_right = diff(x_right)                                                                    #Calculate dx on the right side
    elseif Rfact > 1.0                                                                              #Refine both
        dx_left  = collect(LinRange(1.0, 1.0 * inv(Rfact), nr[1]-1))                                #Spacing left side
        x_left   = [0; cumsum(dx_left)] * Ri[1] * inv(sum(dx_left))                                 #Applying to left size; ratio cumsum(dx_left)/sum(dx_left) defines the spacing
        dx_left  = diff(x_left)                                                                     #Calculate dx on the left side
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]                                                                           #Length of the right domain
        d = dx_left[end]                                                                            #Last dx on the left side
        R = newton_solver(S, d, nr[2]-1, 1e-8, 200,verbose)                                         #Apply Newton solver to find the new right grid
        dx_right = make_dx_right(R, d, nr[2]-1)                                                     #Calculate new dx on the right side
        x_right  = [Ri[1]; Ri[1] .+ cumsum(dx_right)]                                               #Calculate new nodes on the right side
        dx_right = diff(x_right)                                                                    #Calculate dx on the right side
    else                                                                                            #Refine only right grid
        Rfact    = abs(Rfact)                                                                       #Make sure Rfact is positive
        x_left   = collect(LinRange(0.0, Ri[1], nr[1]-1))                                           #Create vector with nr[1] equally spaced nodes
        dx_left  = diff(x_left)                                                                     #Calculate dx on the left side
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]                                                                           #Length of the right domain
        d = dx_left[end]                                                                            #Last dx on the left side
        R = newton_solver(S, d, nr[2]-1, 1e-8, 200,verbose)                                         #Apply Newton solver to find the new right grid
        dx_right = make_dx_right(R, d, nr[2]-1)                                                     #Calculate new dx on the right side
        x_right  = [Ri[1]; Ri[1] .+ cumsum(dx_right)]                                               #Calculate new nodes on the right side
        dx_right = diff(x_right)                                                                    #Calculate dx on the right side
        if verbose == true
            println("MRefin: - $Rfact ; replaced by $Rfact on the right side, kept equal spacing on the left")
        end
    end
    Sc_left  = dx_left[1]*inv(dx_left[end])                                                         #Scaling factor for left side
    Sc_right = dx_right[end]*inv(dx_right[1])                                                       #Scaling factor for right side
    if verbose == true
        #Print results to check
        println("dx left  : ", dx_left[end])
        println("dx right : ", dx_right[1])
        println("x_left   : ", x_left[end])
        println("x_right  : ", x_right[end])
        println("Refinement factor for left is ",Sc_left)
        println("Refinement factor for right is ",Sc_right)
        println("dx1/dx2 (left)  = ",dx_left[1]*inv(dx_left[2]))
        println("dx2/dx1 (right) = ",dx_right[2]*inv(dx_right[1]))
    end
    return Ri, nr, x_left, x_right, dx_left, dx_right, Sc_left, Sc_right
end

"""
    find_dt(dx1, dx2, V_ip, D_l, D_r, CFL)

Find the important time step `dt` based on the given parameters. The function calculates the time step `dt` based on the
advection and diffusion properties of the system. Usually, advection time scale `dtV` are more dominat than diffusion time
scale `dtD`. However, we included a dumping of dt, if `dt`` > `dtD` to ensure the visibility of diffusion processes within
the code. If the advection velocity `V_ip` is zero, `dtD` is used instead. Units may differ from SI units if
non-dimensionalisation has been performed.


## Arguments
- `dx1`: Left spatial step size next to the interface.
- `dx2`: Right spatial step size next to the interface.
- `V_ip`: The advection velocity in [m/s].
- `D_l`: The diffusion coefficient on the left side in [m²/s].
- `D_r`: The diffusion coefficient on the right side in [m²/s].
- `CFL`: The Courant-Friedrichs-Lewy number.

## Returns
- `dt`: The calculated time step.

"""
function find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
    dt_drop = 0.99
    #Find the important dt------------------------------------
    dtV   = minimum([dx1,dx2]) ^1 * inv(abs(V_ip))                                                  #Advection time
    dtD   = minimum([dx1,dx2]) ^2 * inv(maximum([D_l,D_r]))*50                                      #Diffusion time
    dtV   = dtV * dt_drop                                                                           #Dropped advection time
    
    if V_ip == 0.0                                                                                  #dt for pure diffusion
        dt = dtD * CFL;
    else
        if dtV < dtD                                                                                #dt for pure diffusion
            dt = dtV * CFL;                                                                         #Advection time with CFL
        else
            dt = dtD *CFL                                                                         #Advection time with CFL
        end
    end
    
    #dtV2  = dtV * CFL * 5.0                                                                         #Advection time with CFL
    #dt    = minimum([dtV1,dtV2])                                                                    #Calculate dt
    #if V_ip == 0.0                                                                                  #dt for pure diffusion
    #    dt   = dtD * CFL
    #elseif dt > dtD                                                                                 #Dumping of dt, if advection and diffusion occur
    #    dt   = dtD * CFL *1e4
    #end
    return dt
end

"""
    fill_matrix!(C, x, D, dt, ndim, nels)

fill_matrix! function fills the global matrices L_g and R_g with the corresponding local matrices and vectors. Units may
differ from SI units if non-dimensionalisation has been performed.

# Arguments
- `C`: Composition matrix in [-].
- `x`: Spatial grid  points.
- `D`: Diffusion coefficient in [m²/s].
- `dt`: Time step in [s].
- `ndim`: Geometry factor.
- `nels`: Number of elements.

# Returns
- `L_g`: Global LHS matrix.
- `R_g`: Global RHS vector.
"""
function fill_matrix!(C,x,D,dt,ndim,nels)
    #Check dimensions-----------------------------------------
    if ndim != 1 && ndim != 2 && ndim != 3
        error("Error in geometry (n must be 1,2 or 3).")
    end
    #Reset matrices-------------------------------------------
    L_g     = spzeros(length(x),length(x))                                                          #Size changes every iteration
    R_g     = zeros(length(x),1)                                                                    #Size changes every iteration
    Co      = copy(C)                                                                               #Copy composition values
    _dt     = inv(dt)                                                                               #Inverse of dt
    #Make global matrices-------------------------------------
    for (iel,_) in enumerate(1:nels)
        x_1 = copy(x[iel])                                                                          #Select x[1:end-1]
        x_2 = copy(x[iel+1])                                                                        #Select x[2:end]
        #Define local matrices using linear shape functions---
        if ndim  == 3
            Mloc = [-((x_1 - x_2)*(6*x_1^2 + 3*x_1*x_2 + 1*x_2^2))*inv(30)    -((x_1 - x_2)*(3*x_1^2 + 4*x_1*x_2 + 3*x_2^2))*inv(60);   #Local M (sub-matrix)
                    -((x_1 - x_2)*(3*x_1^2 + 4*x_1*x_2 + 3*x_2^2))*inv(60)    -((x_1 - x_2)*(1*x_1^2 + 3*x_1*x_2 + 6*x_2^2))*inv(30)]
            Kloc = [-(D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))              (D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2));    #Local K (sub-matrix)
                     (D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))             -(D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))]
        elseif ndim == 2
            Mloc  = [-((x_1 - x_2)*(3*x_1 + x_2))*inv(12)     x_2^2*inv(12) - x_1^2*inv(12);        #Local M (sub-matrix)
                        x_2^2*inv(12) - x_1^2*inv(12)      -((x_1 - x_2)*(x_1 + 3*x_2))*inv(12)]
            Kloc  = [-(D*(x_1 + x_2))*inv(2*(x_1 - x_2))     (D*(x_1 + x_2))*inv(2*(x_1 - x_2));    #Local K (sub-matrix)
                      (D*(x_1 + x_2))*inv(2*(x_1 - x_2))    -(D*(x_1 + x_2))*inv(2*(x_1 - x_2))]
        elseif ndim == 1
            Mloc  = [x_2*inv(3) - x_1*inv(3)   x_2*inv(6) - x_1*inv(6);                             #Local M (sub-matrix)
                     x_2*inv(6) - x_1*inv(6)   x_2*inv(3) - x_1*inv(3)]
            Kloc  = [-D*inv(x_1 - x_2)   D*inv(x_1 - x_2);                                          #Local K (sub-matrix)
                      D*inv(x_1 - x_2)  -D*inv(x_1 - x_2)]
        end
        #Local matrices---------------------------------------
        Lloc = Mloc .* _dt .+ Kloc
        Rloc = Mloc .* _dt  * Co[iel:iel+1]
        #Global matrices--------------------------------------
        L_g[iel:iel+1,iel:iel+1]  .+= Lloc
        R_g[iel:iel+1]            .+= Rloc
    end
    return L_g, R_g
end

"""
    linear_interpolation_1D(x, y, x_interp)

Perform linear interpolation in 1D.

# Arguments
- `x::AbstractArray`: The x coordinates of data points.
- `y::AbstractArray`: The y coordinates of data points.
- `x_interp::Real`: The x coordinate at which interpolation is desired.

# Returns
- `y_interp::Real`: The interpolated y value at `x_interp`.

"""
function linear_interpolation_1D(x, y, x_interp)
    if length(x) != length(y)
        #Check if input arrays are of the same length
        error("x data and y data must be of the same length. Please check your inputs.")
    elseif x_interp < minimum(x) || x_interp > maximum(x)
        #Check if x_interp is within the bounds of x
        error("x_interp is out of bounds of x.")
    end
    # Find the interval in which x_interp lies and interpolate
    for (i,_) in enumerate(1:length(x)-1)
        if x_interp >= x[i] && x_interp <= x[i + 1]
            # Linear interpolation formula for a unequally spaced grid
            x0, y0   = x[i], y[i]
            x1, y1   = x[i + 1], y[i + 1]
            N1       = 1 - (x_interp - x0) * inv(x1 - x0)
            N2       =     (x_interp - x0) * inv(x1 - x0)
            y_interp =  N1 * y0 + N2 * y1
            return y_interp
        end
    end
end

"""
    linspace_interface(L1, L2, LIP, nx1, nx2, dX1_dXN)

This function calculates the adaptive grid depending on the position of the interface.
Note: This function assumes that `L1 < L2` and `nx1 <= nx2`.
Units may differ from SI units if non-dimensionalisation has been performed.

## Arguments
- `L1`: Length of the left side in [m].
- `L2`: Length of the right side in [m].
- `LIP`: Interface position in [m].
- `nx1`: Number of nodes on the left side.
- `nx2`: Number of nodes on the right side.
- `dX1_dXN`: Ratio of the first grid spacing to the last grid spacing.

## Returns
- `x_left`: Array of nodes on the left side in [m].
- `x_right`: Array of nodes on the right side in [m].

"""
function linspace_interface(L1,L2,LIP,nx1,nx2,dX1_dXN)
    #This function calculates the adaptive grid depending on the position of
    #the interface.
    #Preallocations-------------------------------------------
    mag_left    = 0.0
    c           = 0.0
    Fa          = zeros(nx2)
    Fb          = zeros(nx2)
    Fc          = zeros(nx2)
    Resa        = zeros(nx2)
    Resb        = zeros(nx2)
    Resc        = zeros(nx2)
    #=
    ==========================================================
    Create Left Grid------------------------------------------
    mag is the magnification factor between sequential dx.
    Note that dx_i+1 = mag*dx_i -> Thus dx_n = mag^(nx-1)*dx_i
    The dxs are found by constructing a system of equations as shown below
    L1: length left side, L2: length right side, LIP: interface position
    |1     1      1|   |dx1|   |L2-L1|
    |-mag  1      0| = |dx2| = |  0  |
    |0     -mag   1| = |dx3| = |  0  |
    =#
    #Solve the SoE from above
    mag_left  = (dX1_dXN) .^ (1*inv((1-nx1)))                                                       #Calculate mag factor for left side
    ndx       = nx1 - 1                                                                             #Number of dxs (left side)
    LHS       = Matrix{Float64}(I,ndx,ndx)                                                          #Preallocate LHS
    LHS[1,:] .= 1.0
    for (i,_) in enumerate(1:ndx-1)
        LHS[i+1,i] = -mag_left                                                                      #Fill LHS
    end
    RHS       = zeros(ndx,1)                                                                        #Preallocate RHS
    RHS[1]    = LIP-L1
    dx_left   = (LHS\RHS)'                                                                          #Calculate dxs
    x1        = [0 cumsum(dx_left, dims = 2);]
    x_left    = (x1 .+ L1)'                                                                         #Calculate new left grid
    dx_left_last  = x_left[end] - x_left[end-1]                                                     #Calculate last dx
    #=
    ----------------------------------------------------------
    Then calculate mag_right based on
    dx1   +    dx2     +   ...   + dxn   = L2-LIP
    dx1   + magdx1     +   mag^(n-1)*dx1 = L2-LIP
    (1+mag+mag^2.+mag^(n-1))*dx1         = L2-LIP
    mag_right is calculated using the bisection method
    =#
    a   = [1e-10]
    b   = [1e10]
    Fa  = 0.0
    Fb  = 0.0
    for (i,_) in enumerate(1:nx2)
        Fa = Fa .+ a .^ (i - 1.0)
        Fb = Fb .+ b .^ (i - 1.0)
    end
    Resa = Fa .* dx_left_last .- (L2 - LIP)
    Resb = Fb .* dx_left_last .- (L2 - LIP)
    Err  = 1e23
    while Err > Float64(1e-13*dx_left_last)
        c  = 0.5 * a + 0.5 * b
        Fc = 0.0
        for (i,_) in enumerate(1:nx2)
            Fc = Fc .+ c .^ (i - 1.0)
        end
        Resc = Fc .* dx_left_last .- (L2 - LIP)
        if dot(Resa,Resc) < 0.0
            b = copy(c)
        else
            a = copy(c)
        end
        Err = float(sum(copy(abs.(Resc))))
    end
    mag_right = float(sum(copy(c)))
    #Calculate Right Grid-------------------------------------
    #Same approach as before with the left grid
    ndx      = nx2 - 1                                                                              #Number of dxs (right side)
    LHS      = Matrix{Float64}(I,ndx,ndx)                                                           #Preallocate LHS
    LHS[1,:] .= 1.0
    for (i,_) in enumerate(1:ndx-1)
        LHS[i+1,i] = - mag_right                                                                    #Fill LHS
    end
    RHS      = zeros(ndx,1)                                                                         #Preallocate RHS
    RHS[1]   = L2 - LIP
    dx_right = (LHS\RHS)'                                                                           #Calculate dxs
    x2       = [0 cumsum(dx_right, dims=2);]                                                        #Calculate x2
    x_right  = (x2 .+ LIP)'                                                                         #Calculate new right grid
    return x_left, x_right
end

"""
    make_dx_right(R, d1, n)

Constructs an array containing the right side `dx`.

# Arguments
- `R::Number`: The common ratio to scale all `dx`.
- `d1::Number`: The initial value of `dx`.
- `n::Integer`: Number of `dx` elements.

# Returns
- `dx::Array{Float64,1}`: Spatial distancenes for the right grid.
"""
function make_dx_right(R, d1, n)
    dx    = zeros(n)
    dx[1] = d1
    for i in 2:n
        dx[i] = dx[i - 1] * R
    end
    return dx
end

"""
    newton_solver(S, d, n, tol, max_iter, verbose)

Solves the non-linear equation F(x) = (1 - x^n) / (1-x) - S / d using the Newton method.

# Arguments
- `S::Float64`: The value of S in the equation (scaling factor).
- `d::Float64`: The value of d in the equation (e.g. spatial distance).
- `n::Int`: The value of n in the equation (exponential factor).
- `tol::Float64`: The tolerance for convergence.
- `max_iter::Int`: The maximum number of iterations.
- `verbose::Bool`: Whether to print iteration information.

# Returns
- `x::Float64`: The solution to the equation.
"""
function newton_solver(S, d, n, tol, max_iter,verbose)
    #Newton Solver--------------------------------------------
    x   = 1.1                                                                                       #CAUTION: Do not change this value (initial guess)
    Res = 1e23                                                                                      #Some large number
    for i in 1:max_iter
        Fx  = (1 - x^n) * inv(1 - x) - (S *inv(d))                                                  #Function to solve
        _dFx = inv(fma(-n * x^(n - 1), (1 - x), (1 - x^n)) *inv((1 - x)^2))                         #Derivative
        Res = abs(Fx)                                                                               #Residual
        if Res < tol
            if verbose == true
                println("Newton converged in $i iterations - Res: $Res")
            end
            return x
        end
        x = x - 1.0* Fx * _dFx                                                                      #New value
        if verbose == true
            println("Newton iteration $i: x = $x with Res: $Res")
        end
    end
    if Res>1e-3
        error("Cannot proceed (Newton Failure) - increase Ri[1] or increase grid resolution.")
    else
        println("Warning Newton RES: ",Res," - increase Ri[1] or increase grid resolution.")
    end
    return x
end

"""
    pchip(x, y, X)

Shape-preserving piecewise Cubic Hermite Interpolating Polynomial.

This function computes the shape-preserving piecewise cubic Hermite interpolating polynomial for
the given data points `(x, y)`.

# Arguments
- `x::Vector`: The x-coordinates of the data points.
- `y::Vector`: The y-coordinates of the data points.
- `X::Vector`: The x-coordinates at which to evaluate the interpolating polynomial.

# Returns
- `P::Vector`: The interpolated values at the points `X`.
"""
function pchip(x,y,X)
    #Shape-preserving piecewise Cubic Hermite Interpolating Polynomial
    #Define parameters----------------------------------------
    msl = 3.0                                                                                       #Mon. Slope Limiter
    n   = length(x)
    @views h   = x[2:end] .- x[1:end-1]
    @views y1  = y[1:end-1]
    @views y2  = y[2:end]
    m   = (y2 .- y1) .* inv.(h)
    #Calculate slopes-----------------------------------------
    d   = zeros(n,1)
    for i in 1:n
        if i == 1
            d[i] = (fma(fma(2.0, h[1], h[2]), m[1], - (h[1] * m[2]))) * inv(h[1] + h[2])
            if m[i] * m[i+1] < 0.0 || m[i] == 0.0 || m[i+1] == 0.0
                d[i] = 0.0
            end
        elseif i == n
            d[i] = (fma(fma(2.0, h[n-1], h[n-2]), m[n-1], - (h[n-1] * m[n-2]))) * inv(h[n-1] + h[n-2])
            if m[i-1] * m[i-2] < 0.0 || m[i-1] == 0.0 || m[i-2] == 0.0
                d[i] = 0.0
            end
        else
            wim1    = 2.0 * inv(x[i]   - x[i-1])
            wi      = 2.0 * inv(x[i+1] - x[i])
            d[i]    = (wim1 * m[i-1] + wi * m[i]) * inv(wim1 + wi)
            if m[i-1] * m[i] < 0.0 || m[i-1] == 0.0 || m[i] == 0.0
                d[i] = 0.0
            end
            d[i] = min(abs(d[i]),msl * min(abs(m[i-1]),abs(m[i]))) * sign(d[i])                     #Monotonicity adjustment
        end
    end
    #Interpolation--------------------------------------------
    P = zeros(length(X),1)
    for (j,_) in enumerate(1:length(X))
        if !isempty(findall(x -> x == X[j],x))
            k = findall(x -> x == X[j],x)
            P[j] = y[k][1]
        elseif X[j] > maximum(x)
            k = maximum(findall(x -> x == maximum(x),x))
            P[j] = y[k][1]
        elseif X[j] < minimum(x)
            k = minimum(findall(x -> x == minimum(x),x))
            P[j] = y[k][1]
        else
            k = maximum(findall(x -> x < X[j],x))
            if k > length(x) - 1
                error("k out of bounds")
            end
            @views x1 = x[k]
            @views x2 = x[k+1]
            @views y1 = y[k]
            @views y2 = y[k+1]
            @views d1 = d[k]
            @views d2 = d[k+1]
            s   = X[j] - x1
            h   = x2   - x1
            P[j] = (      3.0 * h * s ^ 2   .- 2.0 * s ^ 3) * inv(h ^ 3) * y2 .+
                   (h^3 - 3.0 * h * s ^ 2   .+ 2.0 * s ^ 3) * inv(h ^ 3) * y1 .+
                   (s ^ 2  * (s - h) ^ 1)                   * inv(h ^ 2) * d2 .+
                   (s ^ 1  * (s - h) ^ 2)                   * inv(h ^ 2) * d1
        end
    end
    return P
end

"""
    regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, nr, nmin, MRefin, verbose)

Regrid the grid and interpolate the composition profiles. Units may differ from SI units if non-dimensionalisation has
been performed.

# Arguments
- `Fl_regrid::Int`: Flag indicating whether to regrid or not.
- `x_left::Vector`: Vector of left spatial nodes in [m].
- `x_right::Vector`: Vector of right spatial nodes in [m].
- `C_left::Vector`: Vector of left composition values in [-].
- `C_right::Vector`: Vector of right composition values in [-].
- `Ri::Vector`: Radii [interface    total length] in [m].
- `V_ip::Float64`: Velocity of the interface in [m/s].
- `nr::Vector`: Resolution.
- `nmin::Int`: Minimum grid size.
- `MRefin::Int`: Refinement factor.
- `verbose::Bool`: Whether to print additional information.

# Returns
- `x_left::Matrix`: Matrix of left nodes in [m].
- `x_right::Matrix`: Matrix of right nodes in [m].
- `C_left::Vector`: Vector of left composition values in [-].
- `C_right::Vector`: Vector of right composition values in [-].
- `dx1::Float64`: Grid spacing at the left interface in [m].
- `dx2::Float64`: Grid spacing at the right interface in [m].
- `nr::Vector`: Updated resolution.
"""
function regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, nr, nmin, MRefin,verbose)
    x_left  = copy(vec(x_left))
    x_right = copy(vec(x_right))
    if Fl_regrid == 1
        if V_ip > 0.0                                                                               #Store new grid size for positive velocities
            nr[1] = round(Ri[1] * inv(Ri[2] - Ri[1]) * nr[2])                                       #Define resolution for left side
            if nr[1] < nmin[1]                                                                      #Check if resolution is above minimum (left side)
                nr[1] = nmin[1]
            end
            if nr[2] < nmin[2]                                                                      #Check if resolution is above minimum (right side)
                nr[2] = nmin[2]
            end
        elseif V_ip < 0.0                                                                           #Store new grid size for negative velocities
            nr[2] = round(Ri[2] * inv(Ri[2] - Ri[1]) * nr[1])                                       #Define resolution for right side
            if nr[1] < nmin[1]                                                                      #Check if resolution is above minimum (left side)
                nr[1] = nmin[1]
            end
            if nr[2] < nmin[2]                                                                      #Check if resolution is above minimum (right side)
                nr[2] = nmin[2]
            end
        end
        #Calculate new grid
        #X_left, X_right = linspace_interface(0, Ri[2], Ri[1], nr[1], nr[2], MRefin)                #Use the bisection method
        Ri, nr, X_left, X_right = define_new_grid(Ri,nr,MRefin,verbose)                             #Use the Newton method (better performance!)
        dx1     = X_left[end] - X_left[end-1]                                                       #Calculate last dx on the left side
        dx2     = X_right[2] - X_right[1]                                                           #Calculate first dx on the right side
        C_left  = pchip(x_left, C_left, vec(X_left))                                                #Interpolate composition on the left side
        C_right = pchip(x_right, C_right, vec(X_right))                                             #Interpolate composition on the right side
        #Update grid
        x_left  = copy(collect(X_left))
        x_right = copy(collect(X_right))
    else
        dx1 = x_left[end] - x_left[end-1]                                                           #Calculate last dx on the left side
        dx2 = x_right[2] - x_right[1]                                                               #Calculate first dx on the right side
    end
    return x_left, x_right, C_left, C_right, dx1, dx2, nr
end

"""
    rescale(Ri0, Ri_input, x_left_input, x_right_input, x0_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input, Lsc, Dsc, Vsc, tsc)

Rescales various input parameters based on provided scaling factors. Rescaled factors are provided in Si units.

# Arguments
- `Ri0`: Non-dimensionalized initial radii.
- `Ri_input`: Non-dimensionalized input radius.
- `x_left_input`: Non-dimensionalized spatial vector left side.
- `x_right_input`: Non-dimensionalized spatial vector right side.
- `x0_input`: Non-dimensionalized initial position.
- `Di_input`: Non-dimensionalized diffusion coefficient.
- `D0_input`: Non-dimensionalized pre-exponential factor.
- `V_input`: Non-dimensionalized input velocity.
- `t_tot_input`: Non-dimensionalized total time.
- `t_ar_input`: Non-dimensionalized array of time points.
- `Lsc`: Length non-dimensionalization factor.
- `Dsc`: Diffusion non-dimensionalization factor.
- `Vsc`: Velocity non-dimensionalization factor.
- `tsc`: Time non-dimensionalization factor.

# Returns
- `Ri0`: Dimensionalized initial radii in [m].
- `Ri`: Dimensionalized radii in [m].
- `x_left`: Dimensionalized spatial vector left side in [m].
- `x_right`: Dimensionalized spatial vector right side in [m].
- `x0`: Dimensionalized initial spatial position in [m].
- `Di`: Dimensionalized input diffusion coefficient in [m²/s].
- `D0`: Dimensionalized initial diffusion coefficient in [m²/s].
- `V_ip`: Dimensionalized input velocity in [m/s].
- `t_tot`: Dimensionalized total time in [s].
- `t_ar`: Dimensionalized array of time points in [s].
"""

function rescale(Ri0, Ri_input, x_left_input, x_right_input, x0_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input, Lsc, Dsc, Vsc, tsc)
    #Rescaling utilized numbers
    #Dependent numbers
    V_ip    = V_input        * Vsc
    t_tot   = t_tot_input    * tsc
    t_ar    = t_ar_input     * tsc
    D0      = D0_input      .* Dsc
    Di      = Di_input      .* Dsc
    Ri      = Ri_input      .* Lsc
    Ri0     = Ri0           .* Lsc
    x_left  = x_left_input  .* Lsc
    x_right = x_right_input .* Lsc
    x0      = x0_input      .* Lsc
    return Ri0, Ri, x_left, x_right, x0, Di, D0, V_ip, t_tot, t_ar
end

"""
    scaling(Ri_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input)

Non-dimensionalizes the input parameters for a diffusion-coupled growth model.

# Arguments
- `Ri_input::Vector{Float64}`: Initial radii of the phases [m].
- `Di_input::Vector{Float64}`: Diffusion coefficients of the particles [m²/s]. If `Di_input` is `[-1, -1]`, the Arrhenius equation is used to calculate `D`.
- `D0_input::Vector{Float64}`: Pre-exponential factor.
- `V_input::Float64`: Initial velocity [m/s].
- `t_tot_input::Float64`: Total time [s].
- `t_ar_input::Float64`: Array of time points [s].

# Returns
- `V_ip::Float64`: Non-dimensionalized initial velocity.
- `t_tot::Float64`: Non-dimensionalized total time.
- `t_ar::Float64`: Non-dimensionalized array of time points.
- `Di::Vector{Float64}`: Non-dimensionalized diffusion coefficients.
- `D0::Vector{Float64}`: Non-dimensionalized initial diffusion coefficients.
- `Ri::Vector{Float64}`: Non-dimensionalized initial radii.
- `Lsc::Float64`: Length scale.
- `Dsc::Float64`: Diffusion scale.
- `Vsc::Float64`: Velocity scale.
- `tsc::Float64`: Time scale.

# Description
This function performs non-dimensionalization of the input parameters based on the given scales.
The length scale (`Lsc`) is fixed at `1e-3` meters. The diffusion scale is the average of `Di_input` or `D0_input`.
The function then calculates the dependent scales (`tsc`, `Vsc`) and non-dimensionalizes the input parameters accordingly.
"""

function scaling(Ri_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input)
    #Non-dimensionalization of input parameters
    #Independent scales
    Lsc = 1e-3                                                                                      #[m]
    if Di_input == [-1 -1]
        Dsc = (D0_input[1]+D0_input[2]) * inv(2.0)                                                  #[m²/s]
    else
        Dsc = (Di_input[1]+Di_input[2]) * inv(2.0)                                                  #[m²/s]
    end
    #Dependent scales
    tsc   = Lsc^2 * inv(Dsc)
    Vsc   = Dsc   * inv(Lsc)
    #Independent numbers
    V_ip  = V_input     * inv(Vsc)
    t_tot = t_tot_input * inv(tsc)
    t_ar  = t_ar_input  * inv(tsc)
    D0    = D0_input   .* inv(Dsc)
    Di    = Di_input   .* inv(Dsc)
    Ri    = Ri_input   .* inv(Lsc)
    return V_ip, t_tot, t_ar, Di, D0, Ri, Lsc, Dsc, Vsc, tsc
end

"""
    save_figure(save_path::String, save_file::Bool)

Save the current figure to a file if `save_file` is `true`. The file will be saved in the directory specified by `save_path`
with a filename `save_name` that includes the current date and time.

# Arguments
- `save_path::String`: The directory path where the figure will be saved.
- `save_file::Bool`: A flag indicating whether to save the figure or not.
"""

function save_figure(save_name::String = "My_example",save_path::String = "Diff-coupled-growth",save_file::Bool = false)
    if save_file
        current_time = Dates.format(now(), "dd_mm_yy_HHMMSS")
        savefig(joinpath(save_path,"$(save_name)_$(current_time).pdf"))
        savefig(joinpath(save_path,"$(save_name)_$(current_time).png"))
    end
end

"""
    set_inner_bc_flux!(L_g, R_g, KD, D_l, D_r, x_left, x_right, V_ip, rho, nr)

Set the inner boundary conditions at the interface using fluxes. Units may differ from SI units if non-dimensionalisation has
been performed.

# Arguments
- `L_g::Matrix`: The global matrix representing the system of equations (LHS).
- `R_g::Vector`: The global vector representing the right-hand side (RHS) of the system of equations.
- `KD::Float64`: The distribution coefficient.
- `D_l::Float64`: The diffusion coefficient on the left side in [m²/s].
- `D_r::Float64`: The diffusion coefficient on the right side in [m²/s].
- `x_left::Vector`: Left spatial nodes in [m].
- `x_right::Vector`: Right spatial nodes in [m].
- `V_ip::Float64`: The interface velocity in [m/s].
- `rho::Vector`: The density value in [kg/m³].
- `nr::Vector`: Resolution.

# Returns
- `L_g::Matrix`: Updated LHS matrix.
- `R_g::Vector`: Updated RHS vector.
- `ScF::Float64`: The scale factor used to reduce the condition number.

"""
function set_inner_bc_flux!(L_g,R_g,KD,D_l,D_r,x_left,x_right,V_ip,rho,nr)
    #Reduce the condition Number------------------------------
    ScF = 1.0
    #Inner BC1------------------------------------------------
    #KD = C_left/C_right
    L_g[nr[1],:]        .= 0.0
    L_g[nr[1],nr[1]]     = 1.0 * ScF
    L_g[nr[1],nr[1]+1]   = - KD * ScF
    R_g[nr[1]]           = 0.0
    #Inner BC2------------------------------------------------
    #J_left = J_right
    L_g[nr[1]+1,:]      .= 0.0
    L_g[nr[1]+1,nr[1]+1] = (-V_ip*rho[2] + rho[2] * D_r * inv(x_right[2]  - x_right[1]))    * ScF
    L_g[nr[1]+1,nr[1]+2] =              (- rho[2] * D_r * inv(x_right[2]  - x_right[1]))    * ScF
    L_g[nr[1]+1,nr[1]+0] = (+V_ip*rho[1] + rho[1] * D_l * inv(x_left[end] - x_left[end-1])) * ScF
    L_g[nr[1]+1,nr[1]-1] =              (- rho[1] * D_l * inv(x_left[end] - x_left[end-1])) * ScF
    R_g[nr[1]+1]         = 0.0
    return L_g, R_g, ScF
end

"""
    set_inner_bc_mb!(L_g, R_g, dVolC, Mtot, KD, nr)

Set the inner boundary conditions at the interface using mass balance (MB). Units may differ from SI units if non-dimensionalisation has
been performed.

# Arguments
- `L_g::Matrix`: The global matrix representing the system of equations (LHS).
- `R_g::Vector`: The global vector representing the right-hand side (RHS) of the system of equations.
- `dVolC::Vector`: The volume change vector in [m³].
- `Mtot::Float64`: The total mass in [-].
- `KD::Float64`: The distribution coefficient.
- `nr::Vector`: Resolution.

# Returns
- `L_g::Matrix`: Updated LHS matrix.
- `R_g::Vector`: Updated RHS vector.
- `ScF::Float64`: The scale factor used to reduce the condition number.

"""
function set_inner_bc_mb!(L_g,R_g,dVolC,Mtot,KD,nr)
    #Reduce the condition number------------------------------
    ScF = 1.0
    #Inner BC1 (MB)-------------------------------------------
    #Sum(dVolC) = Mtot
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],:]  = copy(dVolC)
    R_g[nr[1]]    = copy(Mtot)
    #Inner BC2 (KD)-------------------------------------------
    #KD = C_left/C_right
    L_g[nr[1]+1,:] .= 0.0
    L_g[nr[1]+1,nr[1]+0] = - 1.0 * ScF
    L_g[nr[1]+1,nr[1]+1] = + KD * ScF
    R_g[nr[1]+1]         = 0.0
    return L_g, R_g, ScF
end

"""
    set_inner_bc_Lasaga!(Cl_i, beta, t, KD, D_r, D_l, D0, C_left, C_right, dx1, dx2, rho, L_g, R_g, nr)

Set inner boundary conditions for the special case of major element diffusion in a diffusion couple (Lasaga, 1983). Units may
differ from SI units if non-dimensionalisation has been performed.

# Arguments
- `Cl_i::Float64`: Initial composition on the left side in [-].
- `beta::Float64`: Variable from Lasagas semi-analytical solution.
- `t::Float64`: Time in [s].
- `KD::Float64`: Partition coefficient.
- `D_r::Float64`: Diffusion coefficient on the right side in [m²/s].
- `D_l::Float64`: Diffusion coefficient on the left side in [m²/s].
- `D0::Array{Float64}`: Pre-exponential factor within the equation for the diffusion coefficient (`D`` at T0).
- `C_left::Array{Float64}`: Array of concentrations on the left side in [-].
- `C_right::Array{Float64}`: Array of concentrations on the right side in [-].
- `dx1::Float64`: Grid spacing on the left side in [m].
- `dx2::Float64`: Grid spacing on the right side in [m].
- `rho::Array{Float64}`: Array of densities in [kg/m³].
- `L_g::Array{Float64}`: Global left-hand side matrix.
- `R_g::Array{Float64}`: Global right-hand side vector.
- `nr::Array{Int64}`: Resolution.

# Returns
- `L_g::Array{Float64}`: Updated global left-hand side matrix.
- `R_g::Array{Float64}`: Updated global right-hand side vector.
- `ScF::Float64`: Scaling factor.
- `BC_left::Float64`: Modelled left inner boundary condition (interface) in [-].
- `BC_right::Float64`: Modelled right inner boundary condition (interface) in [-].
- `BC_left_Las::Float64`: Left inner boundary condition following Lasaga (1983) in [-].
- `BC_right_Las::Float64`: Right inner boundary condition following Lasaga (1983) in [-].
"""
function set_inner_bc_Lasaga!(Cl_i,beta,t, KD,D_r,D_l,D0,C_left,C_right,dx1,dx2,rho,L_g,R_g,nr)
    #Set inner boundary conditions for major elements following Lasaga (1983)
    #Semi-analytical solution
    BC_left_Las  = Cl_i * exp(- beta * t)                                                           #Conentration at the left side of the interface
    C_right_Las  = BC_left_Las * inv(1 - BC_left_Las) * inv(KD)
    BC_right_Las = C_right_Las * inv(1 + C_right_Las)                                               #Conentration at the right side of the interface
    #=Numerical solution--------------------------------------
    Eq1 = (C_left[end]/(1-C_left[end])/(C_right[1]/(1-C_right[1})) == KD)
    Eq2 = (-rho[2}*D0[2]/dx2*(C_right[2] - C_right[1]) == Jr)
    Eq3 = (-rho[1]*D0[1]/dx1*(C_left[end] - C_left[end-1]) == Jl)
    Eq4 = (-rho[1]*D0[1]/dx1*(C_left[end] - C_left[end-1])-(-rho[2]*D0[2]/dx2*(C_right[2] - C_right[1])) == 0)

    1. Solve Eq1 for C_left[end]
    2. Substitute solution of step 1 in Eq4
    3. Solve Eq4 for C_right[1]
    =#
    sol1 = (D_r * dx1 * rho[2] - (C_left[end-1] ^ 2 * D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 -
            2 * C_left[end-1] ^ 2 * D_l ^ 2 * KD * dx2 ^ 2 * rho[1] ^ 2 + C_left[end-1] ^ 2 * D_l ^ 2 * dx2 ^ 2 * rho[1] ^ 2 +
            2 * C_left[end-1] * C_right[2] * D_l * D_r * KD ^ 2 * dx1 * dx2 * rho[1] * rho[2] -
            4 * C_left[end-1] * C_right[2] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] +
            2 * C_left[end-1] * C_right[2] * D_l * D_r * dx1 * dx2 *rho[1] * rho[2] -
            2 * C_left[end-1] * D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 + 2 * C_left[end-1] * D_l ^ 2 * KD * dx2 ^ 2 * rho[1] ^ 2 +
            2 * C_left[end-1] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] - 2 * C_left[end-1] * D_l * D_r * dx1 * dx2 * rho[1] * rho[2] +
            C_right[2] ^ 2 * D_r ^ 2 * KD ^ 2 * dx1 ^ 2 * rho[2] ^ 2 - 2 * C_right[2] ^ 2 * D_r ^ 2 * KD * dx1 ^ 2 * rho[2] ^ 2 +
            C_right[2] ^ 2 * D_r ^ 2 * dx1 ^ 2 *rho[2] ^ 2 - 2 * C_right[2] * D_l * D_r * KD ^ 2 * dx1 * dx2 * rho[1] * rho[2] +
            2 * C_right[2] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] + 2 * C_right[2] * D_r ^ 2 * KD * dx1 ^ 2 * rho[2] ^ 2 -
            2 * C_right[2] * D_r ^ 2 * dx1 ^ 2 * rho[2] ^ 2 + D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 +
            2 * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] + D_r ^ 2 * dx1 ^ 2 * rho[2] ^ 2) ^ (1*inv(2)) + C_left[end-1] * D_l * dx2 * rho[1] +
            C_right[2] * D_r * dx1 * rho[2] + D_l * KD * dx2 * rho[1] - C_left[end-1] * D_l * KD * dx2 * rho[1] -
            C_right[2] * D_r * KD * dx1 * rho[2]) * inv((2 * (D_r * dx1 * rho[2] - D_r * KD * dx1 * rho[2])))
    sol2 = ((C_left[end-1] ^ 2 * D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 - 2 * C_left[end-1] ^ 2 * D_l ^ 2 * KD * dx2 ^ 2 * rho[1] ^ 2 +
            C_left[end-1] ^ 2 * D_l ^ 2 * dx2 ^ 2 * rho[1] ^ 2 +
            2 * C_left[end-1] * C_right[2] * D_l * D_r * KD ^ 2 * dx1 * dx2 * rho[1] * rho[2] -
            4 * C_left[end-1] * C_right[2] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] +
            2 * C_left[end-1] * C_right[2] * D_l * D_r * dx1 * dx2 * rho[1] * rho[2] -
            2 * C_left[end-1] * D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 + 2 * C_left[end-1] * D_l ^ 2 * KD * dx2 ^ 2 * rho[1] ^ 2 +
            2 * C_left[end-1] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] - 2 * C_left[end-1] * D_l * D_r * dx1 * dx2 * rho[1] * rho[2] +
            C_right[2] ^ 2 * D_r ^ 2 * KD ^ 2 * dx1 ^ 2 * rho[2] ^ 2 - 2 * C_right[2] ^ 2 * D_r ^ 2 * KD * dx1 ^ 2 * rho[2] ^ 2 +
            C_right[2] ^ 2 * D_r ^ 2 * dx1 ^ 2 * rho[2] ^ 2 - 2 * C_right[2] * D_l * D_r * KD ^ 2 * dx1 * dx2 * rho[1] * rho[2] +
            2 * C_right[2] * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] +
            2 * C_right[2] * D_r ^ 2 * KD * dx1 ^ 2 * rho[2] ^ 2 - 2 * C_right[2] * D_r ^ 2 * dx1 ^ 2 * rho[2] ^ 2 +
            D_l ^ 2 * KD ^ 2 * dx2 ^ 2 * rho[1] ^ 2 + 2 * D_l * D_r * KD * dx1 * dx2 * rho[1] * rho[2] +
            D_r ^ 2 * dx1 ^ 2 * rho[2] ^ 2) ^ (1*inv(2)) + D_r * dx1 * rho[2] + C_left[end-1] * D_l * dx2 * rho[1] +
            C_right[2] * D_r * dx1 * rho[2] + D_l * KD * dx2 * rho[1] - C_left[end-1] * D_l * KD * dx2 * rho[1] -
            C_right[2] * D_r * KD * dx1 * rho[2]) * ((2 * (D_r * dx1 * rho[2] - D_r * KD * dx1 * rho[2])))
    #Find realistic solution
    if sol2>1 || sol2<0
        BC_right  = copy(sol1)
    else
        BC_right  = copy(sol2)
    end
    BC_left  = C_left[end-1] - D_r * dx1 *inv(D_l * dx2) * (BC_right - C_right[2])
    #Reduce the condition Number------------------------------
    ScF = 1.0
    #Inner BC1------------------------------------------------
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],nr[1]]     = 1.0 * ScF
    R_g[nr[1]]           = BC_left * ScF
    #Inner BC2------------------------------------------------
    L_g[nr[1]+1,:] .= 0.0
    L_g[nr[1]+1,nr[1]+1] = 1.0 * ScF
    R_g[nr[1]+1]         = BC_right * ScF
    return L_g, R_g, ScF, BC_left, BC_right, BC_left_Las, BC_right_Las
end

"""
    set_outer_bc!(BCout, L_g, R_g, C_left, C_right, ScF)

Set the outer boundary conditions (Dirichlet or Neumann) for the diffusion-advection problem. Units may differ from SI units
if non-dimensionalisation has been performed.

# Arguments
- `BCout`: An array indicating the type of boundary condition at the outer boundaries.
- `L_g`: The global left-hand side matrix of the diffusion-advection problem.
- `R_g`: The global right-hand side vector of the diffusion-advection problem.
- `C_left`: Composition vector of the left side in [-].
- `C_right`: Composition vector of the right side in [-].
- `ScF`: A scaling factor.

# Details
- If `BCout[1]` is equal to 1, the left outer boundary condition is Dirichlet, otherwise it is Neumann.
- If `BCout[2]` is equal to 1, the right outer boundary condition is Dirichlet, otherwise it is Neumann.

# Returns
- `L_g`: The updated global left-hand side matrix.
- `R_g`: The updated global right-hand side vector.
"""
function set_outer_bc!(BCout,L_g,R_g,C_left,C_right,ScF)
    #Set outer boundary conditions
    if BCout[1] == 1                                                                                #Enter the loop for Dirichlet BC at the left outer BC, otherwise Neuman
        L_g[1,:]      .=  0.0
        L_g[1,1]       =  1.0 * ScF
        R_g[1]         =  C_left[1] * ScF
    end
    if BCout[2] == 1                                                                                #Enter the loop for Dirichlet BC at the right outer BC, otherwise Neuman
        L_g[end,:]    .=  0.0
        L_g[end,end]   =  1.0 * ScF
        R_g[end]       =  C_right[end] * ScF
    end
    return L_g, R_g
end

"""
    sinusoid_profile(C0, n, L, D, t, G)

Calculates a sinusoidal composition profile.

This function takes the initial composition `C0`, the number of sinusoidal modes `n`, the length of the system `L`,
the diffusion coefficient `D`, the time `t`, and the amplitude `G` as input parameters. It calculates the composition
profile at a given time `t` using the sinusoidal equation. Units may differ from SI units if non-dimensionalisation has
been performed.

# Arguments
- `C0`: Initial composition at position x at `t = 0.0 in [-].
- `n`: Number of sinusoidal modes.
- `L`: Length of the modelling domain in [m].
- `D`: Diffusion coefficient in [m²/s].
- `t`: Time in [s].
- `G`: Amplitude.

# Returns
- `C`: Composition profile at time `t` in [-].
"""
function sinusoid_profile(C0,n,L,D,t,G,x)
    #Calculates a sinusoidal composition profile
    Sin_temp = zeros(length(x))
    for (i,_) in enumerate(1:length(n))
        Sin_temp = Sin_temp .+ exp.(-((n[i] .* pi .* inv(L)) ^ 2) .* D .* t) .*G[i] .* sin.(n[i] .* pi .* x .* inv(L))
    end
    C = C0 .+ Sin_temp
    return C
end

"""
    solve_soe(L_g, R_g, res)

Solves a system of equations.

# Arguments
- `L_g`: The left-hand side matrix of the system of equations.
- `R_g`: The right-hand side vector of the system of equations.
- `res`: Resolution.

# Returns
- `C_left`: The updated composition vector of the left side in [-].
- `C_right`: The updated composition vector of the right side in [-].
"""
function solve_soe(L_g,R_g,res)
    #Solve the system of equations
    CN      = L_g \ R_g
    C_left  = CN[1:res[1]]
    C_right = CN[res[1]+1:end]
    return C_left, C_right
end

"""
    trapezoidal_integration(x, fx)

Trapezoidal integration of `fx` over `x`.

# Arguments
- `x`: Array of x values.
- `fx`: Array of f(x) values.

# Returns
- `Int`: The integrated value.
"""
function trapezoidal_integration(x,fx)
    #Trapezoidal integration of f(x) over x
    nx  = length(x) - 1                                                                             #Number of intervals
    dx  = zeros(nx,1)                                                                               #Preallocation to store dx
    Int = zeros(nx,1)                                                                               #Preallocation to store integrals
    #Calculate dx and Int
    for (i,_) in enumerate(dx)
        dx[i]  = x[i+1] - x[i]
    end
    for (i,_) in enumerate(Int)
        Int[i] = 0.5 * dx[i] * (fx[i] + fx[i+1])
    end
    return cumsum(Int, dims = 1)[end]
end

"""
    update_t_dependent_param!(D0, Di, Ea1, Ea2, KD_ar, R, T_ar, t_ar, t, t_tot)

Update the time dependent parameters `D_l`, `D_r`, `KD`, and `T` based on the given inputs.
If Di = [-1.0 -1.0], the diffusion coefficient will be calculated based on the Arrhenius relation. Units may differ from SI
units if non-dimensionalisation has been performed.

# Arguments
- `D0::Vector{Float64}`: Pre-exponential factor within the calculation of the diffusion coefficient in [m²/s].
- `Di::Vector{Float64}`: Diffusion coefficients in [m²/s].
- `Ea1::Float64`: Activation energy for the left phase in [J/mol].
- `Ea2::Float64`: Activation energy for the right phase in [J/mol].
- `KD_ar::Vector{Float64}`: Array of distributioncoefficients.
- `R::Float64`: Gas constant in [J/(mol*K)].
- `T_ar::Vector{Float64}`: Array of temperatures in [K].
- `t_ar::Vector{Float64}`: Array of time in [s].
- `t::Float64`: Current time in [s].
- `t_tot::Float64`: Total time in [s].

# Returns
- `D_l::Float64`: Updated diffusion coefficient for the left side in [m²/s].
- `D_r::Float64`: Updated diffusion coefficient for the right side in [m²/s].
- `KD::Float64`: Updated distribution coefficient.
- `T::Float64`: Updated temperature in [K].
"""
function update_t_dependent_param!(D0,Di,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
    #Interpolate KD and T-------------------------------------
    KD  = linear_interpolation_1D(t_ar,KD_ar,t)                                                     #New distribution coefficient
    T   = linear_interpolation_1D(t_ar,T_ar,t)                                                      #New temperature
    #Check for boundaries-------------------------------------
    tol = 1e-12
    if t == t_tot
        T  = copy(T_ar[end])
        KD = copy(KD_ar[end])
        println("t = t_tot -> T and KD are set to the last values of the arrays.")
    elseif minimum(KD) < minimum(KD_ar) - tol || maximum(KD) > maximum(KD_ar) + tol
        @show KD
        error("Partition coefficient exceeds the values of KD_ar.")
    elseif minimum(T) < minimum(T_ar) - tol || maximum(T) > maximum(T_ar) + tol
        @show T
        error("Temperature exceeds the values of T_ar.")
    end
    #Calculate D----------------------------------------------
    if Di == [-1 -1]                                                                                #Use Arrhenius relation
        D_l = D0[1] * exp(-Ea1 * inv(R * T))                                                        #Diffusion coefficient left
        D_r = D0[2] * exp(-Ea2 * inv(R * T))                                                        #Diffusion coefficient right
    elseif Di != [-1 -1]                                                                            #Use initial diffusivities
        D_l = Di[1]                                                                                 #Diffusion coefficient left
        D_r = Di[2]                                                                                 #Diffusion coefficient right
    else
        error("Check initial diffusivities. If you want to use Di, the values should unequal to -1. If you want to calculate D based on the Arrhenius relation both entries of Di should be -1.")
    end
    return D_l, D_r, KD, T
end

"""
    update_t_dependent_param_simple!(D0, Di, Ea1, R, T_ar, t_ar, t, t_tot)

Update the dependent parameters `D` and `T` based on the given inputs. If Di = [-1.0], the
diffusion coefficient will be calculated based on the Arrhenius relation. Units may differ from SI units if
non-dimensionalisation has been performed.

# Arguments
- `D0`: Pre-exponential factor within the calculation of the diffusion coefficient in [m²/s].
- `Di`: Initial diffusion coefficient in [m²/s].
- `Ea1`: Activation energy in [J/mol].
- `R`: Gas constant in [J/(mol*K)].
- `T_ar`: Array of temperatures in [K].
- `t_ar`: Array of time in [s].
- `t`: Current time in [s].
- `t_tot`: Total time in [s].

# Returns
- `D_l`: Updated diffusion coefficient in [m²/s].
- `T`: Updated temperature in [K].

"""
function update_t_dependent_param_simple!(D0,Di,Ea1,R,T_ar,t_ar,t,t_tot)
    #Interpolate T--------------------------------------------
    T   = linear_interpolation_1D(t_ar,T_ar,t)
    #Check for boundaries-------------------------------------
    tol = 1e-12
    if t == t_tot
        T  = copy(T_ar[end])
        println("t = t_tot -> T is set to the last values of the arrays.")
    elseif minimum(T) < minimum(T_ar) - tol || maximum(T) > maximum(T_ar) + tol
        @show T
        error("Temperature exceeds the values of T_ar.")
    end
    #Calculate D----------------------------------------------
    if Di == -1                                                                                     #Use Arrhenius relation
        D_l = D0 * exp(-Ea1 * inv(R * T))                                                           #Diffusion coefficient
    elseif Di != -1                                                                                 #Use initial diffusivity
        D_l = Di                                                                                    #Diffusion coefficient
    else
        error("Check initial diffusivity. If you want to use Di, the value should unequal to -1. If you want to calculate D based on the Arrhenius relation Di should be -1.")
    end
    return D_l, T
end

"""
    update_time!(t, dt, it, t_tot)

Update the time and time  related variables `t`, `dt`, and `it` for a given total time `t_tot`. Units may differ from SI
units if non-dimensionalisation has been performed.

# Arguments
- `t::Number`: Current time in [s].
- `dt::Number`: Time step in [s].
- `it::Integer`: Time iterations.
- `t_tot::Number`: Total time in [s].

# Returns
- `t::Number`: Updated time in [s].
- `dt::Number`: Updated time step in [s].
- `it::Integer`: Updated time iterations.
"""
function update_time!(t,dt,it,t_tot)
    t   = t + dt                                                                                    #New time
    it  = it + 1                                                                                    #Time iterations
    if t > t_tot                                                                                    #Correct for over shooting
        dt  = dt - (t - t_tot)
        t   = t_tot
    end
    return t, dt, it
end
