using LinearAlgebra
export pchip_interpolation, pchip_evaluate
# Define the PCHIP function
function pchip_interpolation(x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real
    n = length(x)
    if n != length(y)
        throw(ArgumentError("x and y must have the same length"))
    elseif n < 2
        throw(ArgumentError("At least two points are needed for interpolation"))
    end
    
    # Calculate the differences and slopes
    h = diff(x)
    delta = diff(y) ./ h  # Slope of the secant line between points

    # Compute derivatives
    m = zeros(T, n)
    for i in 2:n-1
        if sign(delta[i-1]) * sign(delta[i]) > 0  # Same sign -> monotonic
            w1 = 2h[i] + h[i-1]
            w2 = h[i] + 2h[i-1]
            m[i] = (w1 + w2) > 0 ? (w1 * delta[i-1] + w2 * delta[i]) / (w1 + w2) : 0.0
        end
    end

    # End conditions (use one-sided differences)
    m[1] = delta[1]
    m[end] = delta[end]

    # Return the interpolating function
    return xi -> pchip_evaluate(xi, x, y, h, delta, m)
end

# Helper function to evaluate the PCHIP at a given point
function pchip_evaluate(xi, x, y, h, delta, m)
    n = length(x)
    # Find the interval containing xi
    idx = searchsortedfirst(x, xi) - 1
    idx = max(min(idx, n - 1), 1)  # Clamp idx to [1, n-1]

    # Hermite cubic interpolation formula
    t = (xi - x[idx]) / h[idx]
    h00 = (1 + 2t) * (1 - t)^2
    h10 = t * (1 - t)^2
    h01 = t^2 * (3 - 2t)
    h11 = t^2 * (t - 1)

    return h00 * y[idx] +
           h10 * h[idx] * m[idx] +
           h01 * y[idx + 1] +
           h11 * h[idx] * m[idx + 1]
end
#=------------------------------------------------------------------------------------------
# Use LinRange inputs
res = 100
x = LinRange(0, 10, res)      # x values: [0.0, 2.5, 5.0, 7.5, 10.0]
y = LinRange(1, 3, res) .^ 2  # y values: [1.0, 1.5625, 4.0, 7.5625, 9.0]

pch = pchip_interpolation(x, y)
xi = 4.0
println("Interpolated value at x = $xi: ", pch(xi))

# Evaluate at multiple points for plotting
xs = LinRange(0, 10, res-10)
ys = [pch(xi) for xi in xs]

Int1 = trapezoidal_integration(x,y)
Int2 = trapezoidal_integration(xs,ys)
println("Integral of the original data: $Int1 and integral of the interpolated data: $Int2")

# Plot
plot(xs, ys, label="PCHIP Interpolation", lw=2)
scatter!(x, y, label="Data Points", color=:red)=#
