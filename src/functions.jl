#= Utilized functions in the framework of the Paper:

by A.Stroh and E. Moulas
doi:
Version: 1.0
=#
using SparseArrays, LinearAlgebra, BenchmarkTools, Revise, Dates, Plots
export advect_interface_regrid!, blkdiag, calculate_dt, calc_mass_vol, calc_mass_vol_simple_diff, calc_volume, create_grid!, find_dt, fill_matrix!, linear_interpolation_1D, linspace_interface, preallocations, regrid!, set_inner_bc_mb!, set_inner_bc_flux!,set_inner_bc_Lasaga!, set_outer_bc!, trapezoidal_integration, update_time!, update_t_dependent_param!, update_t_dependent_param_simple!, construct_matrix_fem, solve_soe,calc_mass_err, make_dx_right, newton_solver, define_new_grid, sinusoid_profile,save_figure,scaling,rescale
#Functions----------------------------------------------------
"""
    advect_interface_regrid!(Ri, V_ip, dt, x_left, x_right, C_left, C_right, nr)

Update the interface position and calculate new grids based on the advection velocity.

# Arguments
- `Ri::Float64`: Radii [interface    total length] in [m].
- `V_ip::Float64`: Advection velocity in [m/s].
- `dt::Float64`: Time step in [s].
- `x_left::Vector{Float64}`: Left grid points.
- `x_right::Vector{Float64}`: Right grid points.
- `C_left::Vector{Float64}`: Concentration values on the left grid points in [mol].
- `C_right::Vector{Float64}`: Concentration values on the right grid points in [mol].
- `nr::Vector{Int}`: Resolution of the left and the right grid.

# Returns
- `Fl_regrid::Int`: Flag indicating if regridding was performed (1) or not (0).
- `x_left::Vector{Float64}`: Updated left grid points.
- `x_right::Vector{Float64}`: Updated right grid points.
- `C_left::Vector{Float64}`: Updated concentration values on the left grid points.
- `C_right::Vector{Float64}`: Updated concentration values on the right grid points.
- `nr::Vector{Int}`: Updated resolution.
- `Ri::Float64`: Updated radii [interface    total length].
"""
function advect_interface_regrid!(Ri,V_ip,dt,x_left,x_right,C_left,C_right,nr)
    Rio   = copy(Ri)
    Ri[1] = Rio[1] + V_ip * dt                           #Update interface position
    if V_ip > 0                                     #Calculate new grid for positive velocity
        x_left      = [x_left; x_left[end]]
        C_left      = [C_left; C_left[end]]
        x_left[end] = Ri[1]
        dx1         = x_left[end] - x_left[end-1]
        if Ri[1] > x_right[2] - x_right[1] + Rio[1]                     #Check if Ri moved to fast
            @show x_right[1] x_right[2] Ri[1] 
            error("Interface moved too fast. Ri is larger than x_right[2].")
        end
        x_right[1]  = Ri[1]
        dx2         = x_right[2] - x_right[1]
        Fl_regrid   = 1
        nr[1]       = nr[1] + 1
    elseif V_ip < 0                                 #Calculate new grid for negative velocity
        x_right     = [x_right[1]; x_right]
        C_right     = [C_right[1]; C_right]
        if Ri[1] < -(x_left[end] - x_left[end-1]) + Rio[1]                       #Check if Ri moved to fast
            @show x_left[end-1] x_left[end] Ri[1]
            error("Interface moved too fast. Ri is smaller than x_left[end-1].")
        end
        x_left[end] = Ri[1]
        dx1         = x_left[end] - x_left[end-1]
        x_right[1]  = Ri[1]
        dx2         = x_right[2] - x_right[1]
        Fl_regrid   = 1
        nr[2]       = nr[2] + 1
    else                                            #Leave grid as it was, if velocityis 0
        Fl_regrid   = 0
        dx1         = x_left[end] - x_left[end-1]
        dx2         = x_right[2] - x_right[1]
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
    n1 = length(R1);
    n2 = length(R2);
    Lblock = spzeros(n1+n2,n1+n2)
    Rblock = zeros(n1+n2)
    for i in 1:n1
       if i == 1
        Lblock[i,i] = L1[1,1]
        Lblock[i,2] = L1[1,2]
       elseif i == n1
        Lblock[i,i]   = L1[n1,n1]
        Lblock[i,i-1] = L1[n1,n1-1]        
       else
        Lblock[i,i-1] = L1[i,i-1];
        Lblock[i,i+0] = L1[i,i+0];
        Lblock[i,i+1] = L1[i,i+1];
       end
       Rblock[i] = R1[i]
    end
    for i in 1:n2
        if i == 1
            Lblock[i+n1,i+n1] = L2[1,1]
            Lblock[i+n1,2+n1] = L2[1,2]
           elseif i == n2
            Lblock[i+n1,i+n1]   = L2[n2,n2]
            Lblock[i+n1,i-1+n1] = L2[n2,n2-1]        
           else
            Lblock[i+n1,i-1+n1] = L2[i,i-1];
            Lblock[i+n1,i+0+n1] = L2[i,i+0];
            Lblock[i+n1,i+1+n1] = L2[i,i+1];
        end
        Rblock[i+n1] = R2[i]
    end
    return Lblock, Rblock
end

"""
    calculate_dt(D, dx, CFL)

Calculate the time step `dt` for a diffusion process.

# Arguments
- `D::Float64`: The diffusion coefficient in [m^2/s].
- `dx::Float64`: The spatial step size in [m].
- `CFL::Float64`: The Courant-Friedrichs-Lewy (CFL) number.

# Returns
- `dt::Float64`: The calculated time step.
"""
function calculate_dt(D,dx,CFL)
    #Calculate dt----------------------------------------------------------
    dt = CFL * (dx .^ 2) *inv(D)
    return dt
end

"""
    calc_mass_err(Mass, Mass0)

Calculate the mass error between the final mass `Mass[end]` and the initial mass `Mass0`.

# Arguments
- `Mass::Vector`: A vector containing the mass values in [mol].
- `Mass0::Number`: The initial mass value in [mol].

# Output
- `ErrM::Number`: The calculated mass error.

"""
function calc_mass_err(Mass,Mass0)
    ErrM = (Mass[end] - Mass0) * inv(Mass0)
    println("The total mass difference is $(ErrM[end]*100)% (relevant for closed systems).")
end

"""
    calc_mass_vol(x_left, x_right, C_left, C_right, n, rho)

Calculate the total mass based on the volume.

# Arguments
- `x_left::Float64`: Left grid points.
- `x_right::Float64`: Right grid points.
- `C_left::Vector{Float64}`: Concentration values of the left phase in [mol].
- `C_right::Vector{Float64}`: Concentration values of the right phase in [mol].
- `n::Int`: Number which defines the geometry.
- `rho::Vector{Float64}`: Densities of the left and right phase [g/mol].

# Returns
- `Mtot::Float64`: The total mass.
"""
function calc_mass_vol(x_left,x_right,C_left,C_right,n,rho)
    #Calculate total mass based on the volume---------------------------------------------
    V_left_ini, V_right_ini, dVC_ini = calc_volume(x_left,x_right,n)
    M_left  = trapezoidal_integration(rho[1]*V_left_ini,C_left)
    M_right = trapezoidal_integration(rho[2]*V_right_ini,C_right)
    Mtot    = M_left + M_right
    return Mtot
end

function calc_mass_vol_simple_diff(x1,C1,ndim,rho)
    #Calculate total mass based on the volume---------------------------------------------
    #preallocations
    V = zeros(length(x1),1)
    dV = zeros(length(x1),1)
    #Calculate volumes
    for (i,_) in enumerate(1:length(x1))
        V[i]   = x1[i] ^ ndim			 #Volume left phase
    end
    for (i,_) in enumerate(1:length(V)-1)
        dV[i] = (V[i+1] - V[i])          #Volume change left phase
    end
    dV1 = [dV; 0] * 0.5
    dV2 = [0; dV] * 0.5
    dVC    = dV1 + dV2           #Total volume change
    Mtot  = trapezoidal_integration(rho[1]*V,C1)
    return Mtot
end

"""
    calc_volume(x1, x2, ndim)

Calculation of all volumes. The density in both phases is constant,
which is why the density does not need to be included in the calculations.

# Arguments
- `x1`: Grid points of the left phase.
- `x2`: Grid points of the right phase.
- `ndim`: Dimension of the system (representing geometry).

# Returns
- `V1`: Array of volumes for the left phase.
- `V2`: Array of volumes for the right phase.
- `dVC`: Array of total volume changes.

"""
function calc_volume(x1,x2,ndim)
    #=Calculation of all volumes. We assume a constant density in both phases,
    which is why the density does not need to be included in the calculations.
    ==========================================================================
    =#
    #preallocations
    V1 = zeros(length(x1),1)
    V2 = zeros(length(x2),1)
    dV = zeros(length(x1)+length(x2)-1,1)
    #Calculate volumes
    for (i,_) in enumerate(1:length(x1))
        V1[i]   = x1[i] ^ ndim			 #Volume left phase
    end
    for (i,_) in enumerate(1:length(x2))
        V2[i]   = x2[i] ^ ndim			 #Volume right phase
    end
    V = [V1; V2]                         #Total volume
    for (i,_) in enumerate(1:length(V)-1)
        dV[i] = (V[i+1] - V[i])          #Volume change left phase
    end
    dV1 = [dV; 0] * 0.5
    dV2 = [0; dV] * 0.5
    dVC    = dV1 + dV2           #Total volume change
    return V1, V2, dVC
end

"""
    construct_matrix_fem(x_left, x_right, C_left, C_right, D_l, D_r, dt, n, nels_l, nels_r, Mloc, Kloc, Lloc, res)

Constructs the global matrix for the FEM solver in a diffusion couple problem.

# Arguments
- `x_left::Vector{Float64}`: Left grid points.
- `x_right::Vector{Float64}`: Right grid points.
- `C_left::Vector{Float64}`: Concentration values of the left phase in [mol].
- `C_right::Vector{Float64}`: Concentration values of the right phase in [mol].
- `D_l::Float64`: Diffusion coefficient on the left side in [m^2/s].
- `D_r::Float64`: Diffusion coefficient on the right side in [m^2/s].
- `dt::Float64`: Time step size in [s].
- `n::Int`: Geometry definition.
- `nels_l::Int`: Number of elements on the left side.
- `nels_r::Int`: Number of elements on the right side.
- `Mloc::Matrix{Float64}`: Local mass matrix.
- `Kloc::Matrix{Float64}`: Local stiffness matrix.
- `Lloc::Matrix{Float64}`: Local load matrix.
- `res::Vector{Float64}`: Resolution.

# Returns
- `L_g::SparseMatrixCSC{Float64, Int}`: Global stiffness matrix.
- `R_g::Vector{Float64}`: Global RHS vector.
- `Co_l::Vector{Float64}`: Stores left side concentration values before the update.
- `Co_r::Vector{Float64}`: Stores right side concentration values before the update.
"""
function construct_matrix_fem(x_left,x_right,C_left,C_right,D_l,D_r,dt,n,nels_l,nels_r,Mloc,Kloc,Lloc,res)
    #store old values for RHS----------------------------------------------
    Co_l = copy(C_left)
    Co_r = copy(C_right)
    #Left side matrix---------------------------------------------------
    #Fill matrix---------------------------------------------------
    L1_g, R1_g = fill_matrix!(C_left,x_left,D_l,dt,n,nels_l)
    #Right side matrix---------------------------------------------------
    #Fill matrix---------------------------------------------------
    L2_g, R2_g = fill_matrix!(C_right,x_right,D_r,dt,n,nels_r)
    #Merge systems---------------------------------------------------
    L_g = blkdiag(L1_g,L2_g)
    R_g = [R1_g; R2_g]
    #L_g2,R_g2 = blocktest(L1_g,R1_g,L2_g,R2_g)
    #check_block = maximum(abs.(L_g2[:] .- L_g[:]));
    #println("Error on blk ",check_block)
    return L_g, R_g, Co_l, Co_r
end

"""
    create_grid!(Ri, nr, MRefin)

Create grid with or without variable spacing.

# Arguments
- `Ri::Vector{Float64}`: Initial radii [interface    total length] in [m].
- `nr::Vector{Int}`: Resolution vector.
- `MRefin::Int`: The refinement factor for the grid.

# Returns
- `x_left::Matrix{Float64}`: Left grid points.
- `x_right::Matrix{Float64}`: Right grid points.
- `dx1::Float64`: The grid spacing on the left side.
- `dx2::Float64`: The grid spacing on the right side.
- `x0::Matrix{Float64}`: Initial grid spacing for the whole domain.

"""
function create_grid!(Ri,nr,MRefin,verbose)
    Ri, nr, x_left, x_right, dx_left, dx_right = define_new_grid(Ri,nr,MRefin,verbose)
    dx1 = x_left[end] - x_left[end-1]
    dx2 = x_right[2]  - x_right[1]
    x0  = [copy(x_left); copy(x_right)]
    return x_left, x_right, dx1, dx2, x0
end

"""
    define_new_grid(Ri, nr, Rfact, verbose)

This function defines a new grid based on the given parameters.

## Arguments
- `Ri`: Radii [interface    total length] in [m]..
- `nr`: Resolution grid points on the left and right sides.
- `Rfact`: Grid refinement factor.
- `verbose`: A boolean indicating whether to print additional information.

## Returns
- `Ri`: Radii [interface    total length] in [m]..
- `nr`: Resolution grid points on the left and right sides.
- `x_left`:Left grid points.
- `x_right`: Right grid points.
- `dx_left`: Left grid spacing.
- `dx_right`: Right grid spacing.
- `Sc_left`: The scaling factor for the left side.
- `Sc_right`: The scaling factor for the right side.

"""
function define_new_grid(Ri,nr,Rfact,verbose)
    # Parameters
    if Rfact == 1.0         #Equally spaced grid
        x_left   = collect(LinRange(0.0, Ri[1], nr[1]))
        x_right  = collect(LinRange(Ri[1],Ri[2], nr[2]))
        dx_left  = diff(x_left)
        dx_right = diff(x_right)
    elseif Rfact > 1.0      #Refine both
        dx_left = collect(LinRange(1.0, 1.0 * inv(Rfact), nr[1]-1))     #spacing 
        x_left = [0; cumsum(dx_left)] * Ri[1] * inv(sum(dx_left))       #appling to left size ratio cumsum(dx_left)/sum(dx_left) is the spacing
        dx_left = diff(x_left)
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]
        d = dx_left[end]
        R = newton_solver(S, d, nr[2]-1, 1e-8, 200,verbose)
        dx_right = make_dx_right(R, d, nr[2]-1)
        x_right = [Ri[1]; Ri[1] .+ cumsum(dx_right)]
        dx_right = diff(x_right)
    else    #Refine only right grid
        Rfact   = abs(Rfact)
        x_left  = collect(LinRange(0.0, Ri[1], nr[1]-1))
        #x_right = collect(LinRange(Ri[1], Ri[2], nr[2]))
        dx_left = diff(x_left)
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]
        d = dx_left[end]
        R = newton_solver(S, d, nr[2]-1, 1e-8, 200,verbose)
        dx_right = make_dx_right(R, d, nr[2]-1)
        x_right  = [Ri[1]; Ri[1] .+ cumsum(dx_right)]
        dx_right = diff(x_right)
        if verbose == true
            println("MRefin: - $Rfact ; replaced by $Rfact on the right side, kept equal spacing on the left")
        end
    end
    Sc_left  = dx_left[1]*inv(dx_left[end])
    Sc_right = dx_right[end]*inv(dx_right[1])
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
advection and diffusion properties of the system. Usually, advection time scale `dtV` are more dominat than diffusion time scale `dtD`.
If the advection velocity `V_ip` is zero, `dtD` is used instead.


## Arguments
- `dx1`: Left spatial step size.
- `dx2`: Right spatial step size.
- `V_ip`: The advection velocity in [m/s].
- `D_l`: The diffusion coefficient on the left side in [m^2/s]. 
- `D_r`: The diffusion coefficient on the right side in [m^2/s].
- `CFL`: The Courant-Friedrichs-Lewy number.

## Returns
- `dt`: The calculated time step.

"""
function find_dt(dx1,dx2,V_ip,D_l,D_r,CFL)
    dt_drop = 0.02
    #Find the important dt-------------------------------------------------
    dtV   = minimum([dx1,dx2]) ^1 * inv(abs(V_ip))               #Advection time
    dtD   = minimum([dx1,dx2]) ^2 * inv(maximum([D_l,D_r]))      #Diffusion time
    dtV1  = dtV * dt_drop
    dtV2  = dtV * CFL * 5.0
    dt    = minimum([dtV1,dtV2])
    @show dt dtD
    if V_ip == 0.0
        dt   = dtD * CFL
    elseif dt > dtD
        dt   = dtD * CFL
    end
    return dt
end

"""
    fill_matrix!(C, Co, x, D, dt, ndim, nels, Mloc, Kloc, Lloc, res)

fill_matrix! function fills the global matrices L_g and R_g with the corresponding local matrices and vectors. 

# Arguments
- `C`: Concentration matrix in [mol].
- `Co`: Initial concentration matrix in [mol].
- `x`: Array of node positions.
- `D`: Diffusion coefficient in [m^2/s].
- `dt`: Time step in [s].
- `ndim`: Number of dimensions (geometry).
- `nels`: Number of elements.
- `Mloc`: Local mass matrix.
- `Kloc`: Local stiffness matrix.
- `Lloc`: Local L matrix.
- `res`: Residual vector.

# Returns
- `L_g`: Global LHS matrix.
- `R_g`: Global RHS vector.
"""
function fill_matrix!(C,x,D,dt,ndim,nels)
    #Check dimensions----------------------------------------------------------
    if ndim != 1 && ndim != 2 && ndim != 3
        error("Error in geometry (n must be 1,2 or 3).")
    end
    #Reset matrices------------------------------------------------------------
    L_g     = spzeros(length(x),length(x))                #size changes every iteration
    R_g     = zeros(length(x),1)                          #size changes every iteration
    Co      = copy(C)
    _dt     = inv(dt)
    #Make global matrices------------------------------------------------------
    for (iel,_) in enumerate(1:nels)
        x_1 = copy(x[iel])
        x_2 = copy(x[iel+1])
        #Define local matrices using linear shape functions---------------------
        if ndim  == 3
            Mloc = [-((x_1 - x_2)*(6*x_1^2 + 3*x_1*x_2 + 1*x_2^2))*inv(30)    -((x_1 - x_2)*(3*x_1^2 + 4*x_1*x_2 + 3*x_2^2))*inv(60);       #Local M (sub matrix)
                    -((x_1 - x_2)*(3*x_1^2 + 4*x_1*x_2 + 3*x_2^2))*inv(60)    -((x_1 - x_2)*(1*x_1^2 + 3*x_1*x_2 + 6*x_2^2))*inv(30)]
            Kloc = [-(D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))              (D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2));          #Local K (sub matrix)
                     (D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))             -(D*(x_1^2 + x_1*x_2 + x_2^2))*inv(3*(x_1 - x_2))]
        elseif ndim == 2
            Mloc  = [-((x_1 - x_2)*(3*x_1 + x_2))*inv(12)     x_2^2*inv(12) - x_1^2*inv(12);       #Local M (sub matrix)
                        x_2^2*inv(12) - x_1^2*inv(12)      -((x_1 - x_2)*(x_1 + 3*x_2))*inv(12)]
            Kloc  = [-(D*(x_1 + x_2))*inv(2*(x_1 - x_2))     (D*(x_1 + x_2))*inv(2*(x_1 - x_2));      #Local K (sub matrix)
                      (D*(x_1 + x_2))*inv(2*(x_1 - x_2))    -(D*(x_1 + x_2))*inv(2*(x_1 - x_2))]
        elseif ndim == 1
            Mloc  = [x_2*inv(3) - x_1*inv(3)   x_2*inv(6) - x_1*inv(6);       #Local M (sub matrix)
                     x_2*inv(6) - x_1*inv(6)   x_2*inv(3) - x_1*inv(3)]
            Kloc  = [-D*inv(x_1 - x_2)   D*inv(x_1 - x_2);      #Local K (sub matrix)
                      D*inv(x_1 - x_2)  -D*inv(x_1 - x_2)]
        end
        #Local matrices--------------------------------------------------------
        Lloc = Mloc .* _dt .+ Kloc
        Rloc = Mloc .* _dt  * Co[iel:iel+1]  
        #Global matrices-------------------------------------------------------
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

## Arguments
- `L1`: Length of the left side in [m]
- `L2`: Length of the right side in [m]
- `LIP`: Interface position in [m]
- `nx1`: Number of grid points on the left side
- `nx2`: Number of grid points on the right side
- `dX1_dXN`: Ratio of the first grid spacing to the last grid spacing

## Returns
- `x_left`: Array of grid points on the left side
- `x_right`: Array of grid points on the right side

"""
function linspace_interface(L1,L2,LIP,nx1,nx2,dX1_dXN)
    #This function calculates the adaptive grid depending on the position of
    #the interface.
    #preallocations----------------------------------------------------------
    mag_left    = 0.0
    c           = 0.0
    Fa          = zeros(nx2)
    Fb          = zeros(nx2)
    Fc          = zeros(nx2)
    Resa        = zeros(nx2)
    Resb        = zeros(nx2)
    Resc        = zeros(nx2)
    #=
    ==========================================================================
    Create Left Grid----------------------------------------------------------
    mag is the magnification factor between sequential dx. 
    Note that dx_i+1 = mag*dx_i -> Thus dx_n = mag^(nx-1)*dx_i
    The dxs are found by constructing a system of equations as shown below
    L1: length left side, L2: length right side, LIP: interface position
    |1     1      1|   |dx1|   |L2-L1|
    |-mag  1      0| = |dx2| = |  0  |
    |0     -mag   1| = |dx3| = |  0  | 
    =#
    #Solve the SoE from above
    mag_left  = (dX1_dXN) .^ (1*inv((1-nx1)))       #Calculate mag factor for left side
    ndx       = nx1 - 1                        #Number of dxs
    LHS       = Matrix{Float64}(I,ndx,ndx)          
    LHS[1,:] .= 1.0                       
    for (i,_) in enumerate(1:ndx-1)
        LHS[i+1,i] = -mag_left           
    end
    RHS       = zeros(ndx,1)
    RHS[1]    = LIP-L1
    dx_left   = (LHS\RHS)'       
    x1        = [0 cumsum(dx_left, dims = 2);]
    x_left    = (x1 .+ L1)'                      #Calculate left grid
    dx_left_last  = x_left[end] - x_left[end-1]  #Calculate last dx
    #=
    --------------------------------------------------------------------------
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
    #Calculate Right Grid------------------------------------------------------
    #Same approach as before with the left grid
    ndx      = nx2 - 1
    LHS      = Matrix{Float64}(I,ndx,ndx)
    LHS[1,:] .= 1.0
    for (i,_) in enumerate(1:ndx-1)
        LHS[i+1,i] = - mag_right    
    end
    RHS      = zeros(ndx,1)
    RHS[1]   = L2 - LIP
    dx_left  = (LHS\RHS)'
    x2       = [0 cumsum(dx_left, dims=2);]
    x_right  = (x2 .+ LIP)'
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
    #Newton Solver--------------------------------------
    x   = 1.1                                                       #Do not change this value (initial guess)
    Res = 1e23                                                      #Some large number
    for i in 1:max_iter
        Fx  = (1 - x^n) * inv(1 - x) - (S *inv(d))                         #Function to solve Eq.21 in the paper
        _dFx = inv(fma(-n * x^(n - 1), (1 - x), (1 - x^n)) *inv((1 - x)^2))    #Derivative
        Res = abs(Fx)                                              #Residual
        if Res < tol
            if verbose == true
                println("Newton converged in $i iterations - Res: $Res")
            end
            return x
        end
        x = x - 1.0* Fx * _dFx                                       #New value
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
    #Define parameters---------------------------------------------------
    msl = 3.0           #Mon. Slope Limiter
    n   = length(x)
    @views h   = x[2:end] .- x[1:end-1]
    @views y1  = y[1:end-1]
    @views y2  = y[2:end]
    m   = (y2 .- y1) .* inv.(h)
    #Calculate slopes-----------------------------------------------------
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
            d[i] = min(abs(d[i]),msl * min(abs(m[i-1]),abs(m[i]))) * sign(d[i])     #Monotonicity adjustment
        end
    end
    #Interpolation--------------------------------------------------------
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
    preallocations(x, C, C_left, C_right, res)

Pre-allocates memory for variables used in a diffusion couple simulation.

# Arguments
- `x`: Array of spatial coordinates.
- `C`: Concentration array in [mol].
- `C_left`: Concentration array for the left side of the diffusion couple in [mol].
- `C_right`: Concentration array for the right side of the diffusion couple in [mol].
- `res`: Resolution array.

# Returns
- `Co`: Pre-allocated array for the old concentration `C` in [mol].
- `Co_left`: Pre-allocated array for old `C_left` in [mol].
- `Co_right`: Pre-allocated array for old `C_right` in [mol].
- `dt`: Pre-allocated variable for time step in [s].
- `dx`: Pre-allocated array for spatial step sizes in [m].
- `Kloc`: Pre-allocated 2x2 matrix for local stiffness matrix.
- `Lloc`: Pre-allocated 2x2 matrix for local load matrix.
- `L_g`: Pre-allocated sparse matrix for global load matrix.
- `Mass`: Pre-allocated array for mass values.
- `Mloc`: Pre-allocated 2x2 matrix for local mass matrix.
- `nels`: Number of elements in the diffusion couple.
- `nels_l`: Number of elements on the left side of the diffusion couple.
- `nels_r`: Number of elements on the right side of the diffusion couple.
- `R_g`: Pre-allocated global RHS vector.
- `x_1`: Pre-allocated array for `x` (excluding the last element).
- `x_2`: Pre-allocated array for `x` (excluding the first element).
- `y_interp`: Pre-allocated array for interpolated values.

"""
function preallocations(x,C,C_left,C_right,res)
    #Pre-allocations----------------------------------------------------
    Co          = zeros(size(C))
    Co_left     = zeros(size(C_left))
    Co_right    = zeros(size(C_right))
    dt          = 0.0
    dx          = zeros(length(x) - 1,1)
    Kloc        = zeros(2, 2)
    Lloc        = zeros(2, 2)
    L_g         = spzeros(length(x),length(x))
    Mass        = Float64[]
    Mloc        = zeros(2, 2)
    nels        = length(x) - 1
    nels_l      = res[1] - 1
    nels_r      = res[2] - 1
    R_g         = zeros(length(x),1)
    x_1         = zeros(length(x) - 1,1)
    x_2         = zeros(length(x) - 1,1)
    y_interp    = zeros(1)
    return Co, Co_left, Co_right, dt, dx, Kloc, Lloc, L_g, Mass, Mloc, nels, nels_l, nels_r, R_g, x_1, x_2, y_interp
end

"""
    regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, nr, nmin, MRefin)

Regrid the grid and interpolate the concentration profiles.

# Arguments
- `Fl_regrid::Int`: Flag indicating whether to regrid or not.
- `x_left::Vector`: Vector of left grid points.
- `x_right::Vector`: Vector of right grid points.
- `C_left::Vector`: Vector of left concentration values in [mol].
- `C_right::Vector`: Vector of right concentration values in [mol].
- `Ri::Vector`: Radii [interface    total length] in [m].
- `V_ip::Float64`: Velocity of the interface in [m/s].
- `nr::Vector`: Resoluzion.
- `nmin::Int`: Minimum grid size.
- `MRefin::Int`: Refinement factor.

# Returns
- `x_left::Matrix`: Matrix of left grid points.
- `x_right::Matrix`: Matrix of right grid points.
- `C_left::Vector`: Vector of left concentration values.
- `C_right::Vector`: Vector of right concentration values.
- `dx1::Float64`: Grid spacing at the left interface.
- `dx2::Float64`: Grid spacing at the right interface.
- `nr::Vector`: Updated resolution.
"""
function regrid!(Fl_regrid, x_left, x_right, C_left, C_right, Ri, V_ip, nr, nmin, MRefin,verbose)
    x_left  = copy(vec(x_left))
    x_right = copy(vec(x_right))
    if Fl_regrid == 1
        if V_ip > 0.0            #Store new grid size for positive velocities
            nr[1] = round(Ri[1] * inv(Ri[2] - Ri[1]) * nr[2])
            if nr[1] < nmin[1]
                nr[1] = nmin[1]
            end
            if nr[2] < nmin[2]
                nr[2] = nmin[2]
            end
        elseif V_ip < 0.0        #Store new grid size for negative velocities
            nr[2] = round(Ri[2] * inv(Ri[2] - Ri[1]) * nr[1])
            if nr[1] < nmin[1]
                nr[1] = nmin[1]
            end
            if nr[2] < nmin[2]
                nr[2] = nmin[2]
            end
        end
        #Calculate new grid
        #X_left, X_right = linspace_interface(0, Ri[2], Ri[1], nr[1], nr[2], MRefin)
        Ri, nr, X_left, X_right = define_new_grid(Ri,nr,MRefin,verbose)
        dx1     = X_left[end] - X_left[end-1]
        dx2     = X_right[2] - X_right[1]
        C_left  = pchip(x_left, C_left, vec(X_left))
        C_right = pchip(x_right, C_right, vec(X_right))
        #Update grid
        x_left  = copy(collect(X_left))
        x_right = copy(collect(X_right))
    else
        dx1 = x_left[end] - x_left[end-1]
        dx2 = x_right[2] - x_right[1]
    end
    return x_left, x_right, C_left, C_right, dx1, dx2, nr
end

function rescale(Ri0, Ri_input, x_left_input, x_right_input, x0_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input, Lsc, Dsc, Vsc, tsc)
    #rescaling of numbers
    #dependent numbers
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

function scaling(Ri_input, Di_input, D0_input, V_input, t_tot_input, t_ar_input)
    #non-dimensionalization of input parameters
    #independent scales
    Lsc = 1e-3                              #[m]
    if Di_input == [-1 -1]
        Dsc = 1e-15
        #Dsc = (D0_input[1]+D0_input[2]) * inv(2.0)  #[m^2/s] 
    else
        #Dsc = maximum(Di_input)
        Dsc = (Di_input[1]+Di_input[2]) * inv(2.0)  #[m^2/s] 
    end                                       
    #dependent scales
    tsc   = Lsc^2 * inv(Dsc)
    Vsc   = Dsc   * inv(Lsc)
    #independent numbers
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

Save the current figure to a file if `save_file` is `true`. The file will be saved in the directory specified by `save_path` with a filename `save_name` that includes the current date and time.

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

Set the inner boundary conditions at the interface using fluxes.

# Arguments
- `L_g::Matrix`: The global matrix representing the system of equations (LHS).
- `R_g::Vector`: The global vector representing the right-hand side (RHS) of the system of equations.
- `KD::Float64`: The partition coefficient.
- `D_l::Float64`: The diffusion coefficient on the left side in [m^2/s].
- `D_r::Float64`: The diffusion coefficient on the right side in [m^2/s].
- `x_left::Vector`: Left grid points.
- `x_right::Vector`: Right grid points.
- `V_ip::Float64`: The interface velocity in [m/s].
- `rho::Vector`: The density value in [g/mol].
- `nr::Vector`: Resolution.

# Returns
- `L_g::Matrix`: Updated LHS matrix.
- `R_g::Vector`: Updated RHS vector.
- `ScF::Float64`: The scale factor used to reduce the condition number.

"""
function set_inner_bc_flux!(L_g,R_g,KD,D_l,D_r,x_left,x_right,V_ip,rho,nr)
    @show KD, D_l, D_r, V_ip, rho, nr
    #Reduce the condition Number---------------------------------------
    ScF = 1.0
    #ScF = mean(diag(L_g)) / length(diag(L_g))
    #ScF = maximum(abs.(diag(L_g)))
    #ScF =  norm(L_g,Inf)
    #ScF = norm(L_g,Frobenius)
    #inner BC1---------------------------------------------------------------
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],nr[1]]     = 1.0 * ScF
    L_g[nr[1],nr[1]+1]   = - KD * ScF
    R_g[nr[1]]           = 0.0
    @show L_g[nr[1],nr[1]] L_g[nr[1],nr[1]+1] R_g[nr[1]] 
    #inner BC2---------------------------------------------------------------
    L_g[nr[1]+1,:] .= 0.0
    L_g[nr[1]+1,nr[1]+1] = (-V_ip + rho[2] * D_r * inv(x_right[2] - x_right[1])) * ScF
    L_g[nr[1]+1,nr[1]+2] =       (- rho[2] * D_r * inv(x_right[2] - x_right[1])) * ScF
    L_g[nr[1]+1,nr[1]+0] = (+V_ip + rho[1] * D_l * inv(x_left[end] - x_left[end-1])) * ScF
    L_g[nr[1]+1,nr[1]-1] =       (- rho[1] * D_l * inv(x_left[end] - x_left[end-1])) * ScF
    R_g[nr[1]+1]         = 0.0
    @show L_g[nr[1]+1,nr[1]+1] L_g[nr[1]+1,nr[1]+2] L_g[nr[1]+1,nr[1]+0] L_g[nr[1]+1,nr[1]-1] R_g[nr[1]+1]
    return L_g, R_g, ScF
end

"""
    set_inner_bc_mb!(L_g, R_g, dVolC, Mtot, KD, nr)

Set the inner boundary conditions at the interface using mass balance (MB).

# Arguments
- `L_g::Matrix`: The global coefficient matrix.
- `R_g::Vector`: The global right-hand side vector.
- `dVolC::Vector`: The volume change vector in [m³].
- `Mtot::Float64`: The total mass in [mol].
- `KD::Float64`: The partition coefficient.
- `nr::Vector`: Resolution.

# Returns
- `L_g::Matrix`: Updated LHS matrix.
- `R_g::Vector`: Updated RHS vector.
- `ScF::Float64`: The scale factor used to reduce the condition number.

"""
function set_inner_bc_mb!(L_g,R_g,dVolC,Mtot,KD,nr)
    #Reduce the condition number-------------------------------------------
    #ScF      = sum(diag(L_g)) / length(diag(L_g))
    ScF = 1.0
    #Inner BC1 (MB)--------------------------------------------------------
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],:] = copy(dVolC)
    R_g[nr[1]]   = copy(Mtot)
    #Inner BC2 (KD)--------------------------------------------------------
    L_g[nr[1]+1,:] .= 0.0                      
    L_g[nr[1]+1,nr[1]+0] = - 1.0 * ScF                
    L_g[nr[1]+1,nr[1]+1] = + KD * ScF               
    R_g[nr[1]+1]         = 0.0
    return L_g, R_g, ScF
end

"""
    set_inner_bc_Lasaga!(Cl_i, beta, t, KD, D_r, D_l, D0, C_left, C_right, dx1, dx2, rho, L_g, R_g, nr)

Set inner boundary conditions for the special case of Lasaga (1983).

# Arguments
- `Cl_i::Float64`: Initial concentration of Cl in [mol].
- `beta::Float64`: Pre-defined variables in Lasaga (1983) -> energy. 
- `t::Float64`: Time in [s].
- `KD::Float64`: Partition coefficient.
- `D_r::Float64`: Diffusion coefficient on the right side in [m^2/s].
- `D_l::Float64`: Diffusion coefficient on the left side in [m^2/s].
- `D0::Array{Float64}`: Pre-exponential factor within the equation for the diffusion coefficient (D at T0) in [m^2/s].
- `C_left::Array{Float64}`: Array of concentrations on the left side in [mol].
- `C_right::Array{Float64}`: Array of concentrations on the right side in [mol].
- `dx1::Float64`: Grid spacing on the left side.
- `dx2::Float64`: Grid spacing on the right side.
- `rho::Array{Float64}`: Array of densities in [g/mol].
- `L_g::Array{Float64}`: Global left-hand side matrix.
- `R_g::Array{Float64}`: Global right-hand side vector.
- `nr::Array{Int64}`: Resolution.

# Returns
- `L_g::Array{Float64}`: Updated global left-hand side matrix.
- `R_g::Array{Float64}`: Updated global right-hand side vector.
- `ScF::Float64`: Scaling factor.
- `BC_left::Float64`: Left inner boundary condition (interface).
- `BC_right::Float64`: Right inner boundary condition (interface).
- `BC_left_Las::Float64`: Left inner boundary condition following Lasaga (1983).
- `BC_right_Las::Float64`: Right inner boundary condition following Lasaga (1983).
""" 
function set_inner_bc_Lasaga!(Cl_i,beta,t, KD,D_r,D_l,D0,C_left,C_right,dx1,dx2,rho,L_g,R_g,nr)
    #Set inner boundary conditions for major elements following Lasaga
    #Semi-analytical solution
    BC_left_Las  = Cl_i * exp(- beta * t)                   #Conentration at the left side of the interface
    C_right_Las  = BC_left_Las * inv(1 - BC_left_Las) * inv(KD)     
    BC_right_Las = C_right_Las * inv(1 + C_right_Las)          #Conentration at the right side of the interface
    #Numerical solution-----------------------------------------------
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
    if sol2>1 || sol2<0
        BC_right  = copy(sol1)
    else
        BC_right  = copy(sol2)
    end
    BC_left  = C_left[end-1] - D_r * dx1 *inv(D_l * dx2) * (BC_right - C_right[2])
    #Reduce the condition Number---------------------------------------
    ScF = 1.0
    #inner BC1---------------------------------------------------------------
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],nr[1]]     = 1.0 * ScF
    R_g[nr[1]]           = BC_left * ScF
    #inner BC2---------------------------------------------------------------
    L_g[nr[1]+1,:] .= 0.0
    L_g[nr[1]+1,nr[1]+1] = 1.0 * ScF
    R_g[nr[1]+1]         = BC_right * ScF
    return L_g, R_g, ScF, BC_left, BC_right, BC_left_Las, BC_right_Las
end

"""
    set_outer_bc!(BCout, L_g, R_g, C_left, C_right, ScF)

Set the outer boundary conditions for the diffusion problem.

# Arguments
- `BCout`: An array indicating the type of boundary condition at the outer boundaries.
- `L_g`: The global left-hand side matrix of the diffusion problem.
- `R_g`: The global right-hand side vector of the diffusion problem.
- `C_left`: Concentration vector of the left side in [mol].
- `C_right`: Concentration vector of the right side in [mol].
- `ScF`: A scaling factor.

# Details
- If `BCout[1]` is equal to 1, the left outer boundary condition is Dirichlet, otherwise it is Neumann.
- If `BCout[2]` is equal to 1, the right outer boundary condition is Dirichlet, otherwise it is Neumann.

# Returns
- `L_g`: The updated global left-hand side matrix.
- `R_g`: The updated global right-hand side vector.
"""
function set_outer_bc!(BCout,L_g,R_g,C_left,C_right,ScF)
    #Set boundary conditions

    if BCout[1] == 1            #Enter the loop for Dirichlet BC at the left outer BC, otherwise Neuman
        L_g[1,:]      .=  0.0
        L_g[1,1]       =  1.0 * ScF
        R_g[1]         =  C_left[1] * ScF
    end
    if BCout[2] == 1            #Enter the loop for Dirichlet BC at the right outer BC, otherwise Neuman
        L_g[end,:]      .=  0.0
        L_g[end,end]     =  1.0 * ScF
        R_g[end]         =  C_right[end] * ScF
    end
    return L_g, R_g
end

"""
    sinusoid_profile(C0, n, L, D, t, G)

Calculates a sinusoidal concentration profile.

This function takes the initial concentration `C0`, the number of sinusoidal modes `n`, the length of the system `L`, the diffusion coefficient `D`, the time `t`, and the amplitude `G` as input parameters. It calculates the concentration profile at a given time `t` using the sinusoidal equation.

# Arguments
- `C0`: Initial concentration at position x at time 0.
- `n`: Number of sinusoidal modes.
- `L`: Length of the modelling domain in [m].
- `D`: Diffusion coefficient in [m²/s].
- `t`: Time in [s].
- `G`: Amplitude.

# Returns
- `C`: Concentration profile at time `t`.
"""
function sinusoid_profile(C0,n,L,D,t,G,x)
    #Calculates a sinusoidal concentration profile
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
- `res`: A tuple specifying the number of variables on the left-hand side and right-hand side.

# Returns
- `C_left`: The updated concentration vector of the left side.
- `C_right`: The updated concentration vector of the right side.
"""
function solve_soe(L_g,R_g,res)
    #Solve the system of equations
    #prob = LinearProblem(L_g, R_g)          # Define the linear system
    #sol = solve(prob)                   # Solve the system
    #CN = sol.u                           # Extract the solution
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
- `Int2`: The integrated value.
"""
function trapezoidal_integration(x,fx)
    #Trapezoidal integration of f(x) over x
    nx = length(x) - 1     #number of intervals
    dx  = zeros(nx,1)       #preallocation to store dx
    Int = zeros(nx,1)       #preallocation to store integrals
    #Calculate dx and Int
    for (i,_) in enumerate(dx)
        dx[i]  = x[i+1] - x[i]
    end
    for (i,_) in enumerate(Int)
        Int[i] = 0.5 * dx[i] * (fx[i] + fx[i+1])
    end
    # Int2 = cumsum(Int, dims = 1)
    return cumsum(Int, dims = 1)[end]
end

"""
    update_t_dependent_param!(D0, Di, dt, Ea1, Ea2, KD_ar, R, T_ar, t_ar, t, t_tot)

Update the time dependent parameters `D_l`, `D_r`, `KD`, and `T` based on the given inputs.
If Di = [-1.0 -1.0], the diffusion coefficient will be calculated based on the Arrhenius relation.

# Arguments
- `D0::Vector{Float64}`: Pre-exponential factor within the calculation of the diffusion coefficient in [m^2/s].
- `Di::Vector{Float64}`: Diffusion coefficients in [m^2/s].
- `dt::Float64`: Time step in [s].
- `Ea1::Float64`: Activation energy for the left phase in [J/mol].
- `Ea2::Float64`: Activation energy for the right phase in [J/mol].
- `KD_ar::Vector{Float64}`: Array of partition coefficients.
- `R::Float64`: Gas constant in [J/(mol*K)].
- `T_ar::Vector{Float64}`: Array of temperatures in [K].
- `t_ar::Vector{Float64}`: Array of time in [s].
- `t::Float64`: Current time in [s].
- `t_tot::Float64`: Total time in [s].

# Returns
- `D_l::Float64`: Updated diffusion coefficient for the left side.
- `D_r::Float64`: Updated diffusion coefficient for the right side.
- `KD::Float64`: Updated partition coefficient.
- `T::Float64`: Updated temperature.
"""
function update_t_dependent_param!(D0,Di,dt,Ea1,Ea2,KD_ar,R,T_ar,t_ar,t,t_tot)
    #Interpolate KD and T-------------------------------------------------
    KD  = linear_interpolation_1D(t_ar,KD_ar,t)  #New partition coefficient  
    T   = linear_interpolation_1D(t_ar,T_ar,t)   #New temperature 
    #Check for boundaries-------------------------------------------------
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
    #Calculate D----------------------------------------------------------
    if Di == [-1 -1]
        D_l = D0[1] * exp(-Ea1 * inv(R * T))            #Diffusion coefficient left
        D_r = D0[2] * exp(-Ea2 * inv(R * T))            #Diffusion coefficient right
    elseif Di != [-1 -1]
        D_l = Di[1]                                 #Diffusion coefficient left
        D_r = Di[2]                                 #Diffusion coefficient right
    else
        error("Check initial diffusivities. If you want to use Di, the values should unequal to -1. If you want to calculate D based on the Arrhenius relation both entries of Di should be -1.")
    end
    return D_l, D_r, KD, T
end

"""
    update_t_dependent_param_simple!(D0, Di, dt, Ea1, R, T_ar, t_ar, t, t_tot)

Update the dependent parameters `D` and `T` based on the given inputs. If Di = [-1.0], the 
diffusion coefficient will be calculated based on the Arrhenius relation.

# Arguments
- `D0`: Pre-exponential factor within the calculation of the diffusion coefficient in [m^2/s].
- `Di`: Initial diffusion coefficient in [m^2/s].
- `dt`: Time step [J/(mol*K)].
- `Ea1`: Activation energy in [J/mol].
- `R`: Gas constant in [J/(mol*K)].
- `T_ar`: Array of temperatures in [K].
- `t_ar`: Array of time in [s].
- `t`: Current time in [s].
- `t_tot`: Total time in [s].

# Returns
- `D_l`: Updated diffusion coefficient.
- `T`: Updated temperature.

"""
function update_t_dependent_param_simple!(D0,Di,dt,Ea1,R,T_ar,t_ar,t,t_tot)
    #Interpolate KD and T-------------------------------------------------
    T   = linear_interpolation_1D(t_ar,T_ar,t)   #New temperature 
    #Check for boundaries-------------------------------------------------
    tol = 1e-12
    if t == t_tot
        T  = copy(T_ar[end])
        println("t = t_tot -> T is set to the last values of the arrays.")
    elseif minimum(T) < minimum(T_ar) - tol || maximum(T) > maximum(T_ar) + tol
        @show T
        error("Temperature exceeds the values of T_ar.")
    end
    #Calculate D----------------------------------------------------------
    if Di == -1
        D_l = D0 * exp(-Ea1 * inv(R * T))            #Diffusion coefficient 
    elseif Di != -1 
        D_l = Di                                 #Diffusion coefficient
    else
        error("Check initial diffusivity. If you want to use Di, the value should unequal to -1. If you want to calculate D based on the Arrhenius relation Di should be -1.")
    end
    return D_l, T
end

"""
    update_time!(t, dt, it, t_tot)

Update the time and time  related variables `t`, `dt`, and `it` for a given total time `t_tot`.

# Arguments
- `t::Number`: Current time in [s].
- `dt::Number`: Time step in [s].
- `it::Integer`: Time iterations.
- `t_tot::Number`: Total time in [s].

# Returns
- `t::Number`: Updated time.
- `dt::Number`: Updated time step.
- `it::Integer`: Updated time iterations.
"""
function update_time!(t,dt,it,t_tot)
    t   = t + dt                #New time
    it  = it + 1                #Time iterations
    if t > t_tot                #Correct for over shooting
        dt  = dt - (t - t_tot)
        t   = t_tot
    end
    return t, dt, it
end
