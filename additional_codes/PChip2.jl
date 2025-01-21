using Plots
function pchip(x,y,X)
    #Shape-preserving piecewise Cubic Hermite Interpolating Polynomial
    #Set initial values
    msl = 3.0                           #Mon. Slope Limiter
    n   = length(x)                     #Number of data points
    h   = x[2:end] .- x[1:end-1]        #dx
    y1  = y[1:end-1]
    y2  = y[2:end]
    m   = (y2 .- y1) ./ h               #local slopes
    #Limit slopes
    d   = zeros(n,1)
    for i = 1:n
        if i == 1      
            d[i] = ((2.0 * h[1] + h[2]) * m[1] - h[1] * m[2]) / (h[1] + h[2])
            if m[i] * m[i+1] < 0.0 || m[i] == 0.0 || m[i+1] == 0.0
                d[i] = 0.0
            end
        elseif i == n
            d[i] = ((2.0 * h[n-1] + h[n-2]) * m[n-1] - h[n-1] * m[n-2]) / (h[n-1] + h[n-2])
            if m[i-1] * m[i-2] < 0.0 || m[i-1] == 0.0 || m[i-2] == 0.0
                d[i] = 0.0
            end
        else
            wim1    = 2.0 / (x[i]   - x[i-1])
            wi      = 2.0 / (x[i+1] - x[i])
            d[i]    = (wim1 * m[i-1] + wi * m[i]) / (wim1 + wi)
            if m[i-1] * m[i] < 0.0 || m[i-1] == 0.0 || m[i] == 0.0
                d[i] = 0.0
            end
            d[i] = min(abs(d[i]),msl * min(abs(m[i-1]),abs(m[i]))) * sign(d[i])     #Monotonicity adjustment
        end
    end
    #Calculate interpolated values
    P = zeros(length(X),1)
    for j = 1:length(X)
        if !isempty(findall(x -> x == X[j],x))
            k = findall(x -> x == X[j],x)
            P[j] = y[k][1]
            @show X[j] y[k]
        elseif X[j] > maximum(x) 
            k = maximum(findall(x -> x == maximum(x),x))
            @show k x[k]
            P[j] = y[k][1]
        elseif X[j] < minimum(x) 
            k = minimum(findall(x -> x == minimum(x),x))
            @show k x[k]
            P[j] = y[k][1]
        else
            k = maximum(findall(x -> x < X[j],x))
            if k > length(x) - 1 
                error("k out of bounds")
            end
            x1 = x[k]
            x2 = x[k+1]
            y1 = y[k]
            y2 = y[k+1]
            d1 = d[k]
            d2 = d[k+1]
            s   = X[j] - x1
            h   = x2   - x1
            P[j] = (        3.0 * h * s ^ 2.0   .- 2.0 * s ^ 3.0) / h ^ 3.0 * y2 .+
                   (h^3.0 - 3.0 * h * s ^ 2.0   .+ 2.0 * s ^ 3.0) / h ^ 3.0 * y1 .+
                   (s ^ 2.0  * (s - h) ^ 1.0)                    / h ^ 2.0 * d2 .+
                   (s ^ 1.0  * (s - h) ^ 2.0)                    / h ^ 2.0 * d1
        end
    end
    return P
end



xan = LinRange(1.5,pi*1.5,100)
yan = cos.(xan)
x   = LinRange(1.5,pi*1.5,20)
y   = cos.(x)

#Set interpolation points
X   = [3.5;2.0;pi*1.5;1.5+(x[2]-x[1])*0.4;pi*1.5-(x[end]-x[end-1])*0.2;1.5;8.4;1]#1.5;8.4;]
#X   = [3.5,2.0,pi*1.5,1.5+(x[2]-x[1])*0.4,pi*1.5-(x[end]-x[end-1])*0.2,1.5]
#X = [-0.5] 
P   = pchip(x,y,X)
P2  = vec(P)

plot(xan,yan,label="Analytical")
plot!(x,y,grid=:on,label="original data")
scatter!(X,P2,label="Interpolated")