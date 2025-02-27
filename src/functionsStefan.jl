using BenchmarkTools, Revise, GLMakie, FileIO
export calculate_density, coeff_trans_line, composition, ndgrid, set_inner_bc_stefan!,values_between_known_indices!
#Functions---------------------------------------------------

"""
    calculate_density(X_A, Y_A, rho_left, rho_right, C_leftB, C_rightB, T)

Calculate the density of the phases at a given temperature `T` using interpolation in 2D.

## Arguments
- `X_A`: X-axis values for interpolation (composition X in [mol])
- `Y_A`: Y-axis values for interpolation (temperature T in [K])
- `rho_left`: Left density values for interpolation in [g/mol]
- `rho_right`: Right density values for interpolation in [g/mol]
- `C_leftB`: Left concentration values for interpolation in [mol]
- `C_rightB`: Right concentration values for interpolation in [mol]
- `T`: Temperature at which to calculate the density in [K]

## Returns
- `rho`: Array of normalized densities

"""
function calculate_density(X_A,Y_A,rho_left,rho_right,C_leftB,C_rightB,T)
    rho_LEFT    = interp2(X_A,Y_A,rho_left,C_leftB,T)                  
    rho_RIGHT   = interp2(X_A,Y_A,rho_right,C_rightB,T)
    rho_0       = [copy(rho_LEFT) copy(rho_RIGHT)] 
    rho_norm_L  = rho_LEFT .* inv(rho_0[1])                        #Normalized densities
    rho_norm_R  = rho_RIGHT .* inv(rho_0[1])                       #Normalized densities
    rho         = [copy(rho_norm_L) copy(rho_norm_R)]
    return rho
end


"""
    coeff_trans_line(eq_values)

Extracts coefficients for linear least squares from the input `eq_values`.

# Arguments
- `eq_values`: A 2x3 matrix containing the coefficients for the upper and lower transition lines.

# Returns
- `coeff_up`: A 1x3 vector containing the coefficients for the upper transition line.
- `coeff_do`: A 1x3 vector containing the coefficients for the lower transition line.
"""
function coeff_trans_line(eq_values)
    #Extract coefficients for linear least squares from input (eq_values) 
    coeff_up = eq_values[1,:]'      #Coefficients for X(T) upper transition line
    coeff_do = eq_values[2,:]'      #Coefficients for X(T) lower transition line
    return coeff_up,coeff_do
end

"""
    composition(coeff_up, coeff_do, T)

Compute the composition of two components A and B at a given temperature.

# Arguments
- `coeff_up::Vector{Float64}`: Coefficients for the composition of component B.
- `coeff_do::Vector{Float64}`: Coefficients for the composition of component A.
- `T::Float64`: Temperature at which to compute the composition in [K].

# Returns
- `C_left::Float64`: Composition of component A in [mol].
- `C_right::Float64`: Composition of component B in [mol].
"""
function composition(coeff_up,coeff_do,T)
    C_left  = coeff_do[1] .+ coeff_do[2] .* T .+ coeff_do[3] .* (T) .^ 2      #Composition of A
    C_right = coeff_up[1] .+ coeff_up[2] .* T .+ coeff_up[3] .* (T) .^ 2      #Composition of B
    return C_left, C_right
end

"""
    interp2(X1d, Y1d, Z2d, xi, yi)

Perform 2D bilinear interpolation within bounds.

# Arguments
- `X1d::Vector`: 1D monotonic (increasing) array representing the x-coordinates.
- `Y1d::Vector`: 1D monotonic (increasing) array representing the y-coordinates.
- `Z2d::Matrix`: 2D array representing the values to be interpolated.
- `xi::Vector`: 1D array of x-values to interpolate.
- `yi::Vector`: 1D array of y-values to interpolate.

# Returns
- `zi::Vector`: 1D array of interpolated values corresponding to `xi` and `yi`.

# Notes
- The lengths of `xi` and `yi` must be the same.
- The values in `xi` and `yi` should be within the bounds of `X1d` and `Y1d` respectively.
- The function uses 2D bilinear interpolation to compute the interpolated values.
"""
function interp2(X1d,Y1d,Z2d,xi,yi)
    @assert length(xi) == length(yi)           "lengths of xi, yi must be the same"
    zi = zeros(length(xi)) * NaN
    for (i,_) in enumerate(xi)
        @assert (xi[i] >= X1d[1]) && (xi[i] <= X1d[end]) "Value should be within x bounds to interpolate"
        @assert (yi[i] >= Y1d[1]) && (yi[i] <= Y1d[end]) "Value should be within y bounds to interpolate"
        ix_cart = findlast(X1d .< xi[i])
        iy_cart = findlast(Y1d .< yi[i])
        ix_lin  = LinearIndices(X1d)[ix_cart]
        iy_lin  = LinearIndices(Y1d)[iy_cart]
        dx      = X1d[ix_lin+1] - X1d[ix_lin]
        dy      = Y1d[iy_lin+1] - Y1d[iy_lin]
        w11     = (X1d[ix_lin+1]-xi[i]) * (Y1d[iy_lin+1] - yi[i]  ) * inv(dx*dy)
        w12     = (X1d[ix_lin+1]-xi[i]) * (yi[i] - Y1d[iy_lin]) * inv(dx*dy)
        w21     = (xi[i] - X1d[ix_lin]) * (Y1d[iy_lin+1] - yi[i]) * inv(dx*dy)
        w22     = (xi[i] - X1d[ix_lin]) * (yi[i] - Y1d[iy_lin]) * inv(dx*dy)
        zi[i]   = w11*Z2d[ix_lin,iy_lin] + w12*Z2d[ix_lin,iy_lin+1] + w21*Z2d[ix_lin+1,iy_lin] + w22*Z2d[ix_lin+1,iy_lin+1]
    end
    return zi
end

"""
    set_inner_bc_stefan!(L_g, R_g, C_left, C_right, nr)

Uses the Stefan boundary conditions to set the inner boundary conditions of the system of equations.

# Arguments
- `L_g::Matrix`: The matrix representing the left-hand side of the system of equations.
- `R_g::Vector`: The vector representing the right-hand side of the system of equations.
- `C_left::Vector`: Concentration of the left phase in [mol].
- `C_right::Vector`: Concentration of the right phase in [mol].
- `nr::Vector`: Resolution of the model.

# Returns
- `L_g::Matrix`: The modified matrix `L_g` with Stefan conditions applied.
- `R_g::Vector`: The modified vector `R_g` with Stefan conditions applied.
- `ScF::Float64`: Scaling factor for reducing the condition number of the matrix.
"""
function set_inner_bc_stefan!(L_g,R_g,C_left,C_right,nr)
    #Reduce the condition number-------------------------------------------
    ScF = Float64(maximum(abs.(C_left)))
    #ScF      = sum(diag(L_g)) / length(diag(L_g))
    #ScF = 1.0
    #Apply Dirichlet BC at interface---------------------------------------
    #inner BC1---------------------------------------------------------------
    L_g[nr[1],:] .= 0.0
    L_g[nr[1],nr[1]]     = 1.0 * ScF
    R_g[nr[1]]           = C_left[end] * ScF
    #Inner BC2 (KD)--------------------------------------------------------
    L_g[nr[1]+1,:] .= 0.0                               
    L_g[nr[1]+1,nr[1]+1] = 1.0 * ScF                
    R_g[nr[1]+1]         = C_right[1] * ScF
    return L_g, R_g, ScF
end

function values_between_known_indices!(vec1, vec2, val1, val2)
    idx1 = findfirst(x -> x <= val1, vec1)        #Find indices, where x is smaller than value 1
    idx2 = findlast(x -> x >= val2, vec1)        #Find indices, where x is larger than value 2

    if idx1 === nothing || idx2 === nothing
        error("The starting or end temperature is outside of the defined temperature range.")
    end
    
    start_idx = min(idx1, idx2)
    end_idx = max(idx1, idx2)

    first_val = vec2[start_idx]
    last_val  = vec2[end_idx]

    return first_val,last_val
end