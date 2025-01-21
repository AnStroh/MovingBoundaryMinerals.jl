using Plots, LinearAlgebra
# Newton Solver
function newton_solver(S, d, n, tol, max_iter,verbose)
    x   = 1.1    #Do not change this value
    Res = 1e23   #Some large number
    for i in 1:max_iter
        Fx  = (1 - x^n) / (1 - x) - (S / d)
        dFx = (-n * x^(n - 1) * (1 - x) + (1 - x^n)) / (1 - x)^2
        Res = abs(Fx);
        if Res < tol
            if verbose == true
                println("Newton converged in $i iterations - Res: $Res")
            end
            return x
        end
        x = x - 1.0* Fx / dFx
        if verbose == true
            println("Newton iteration $i: x = $x with Res: $Res")
        end
    end
    if Res>1e-3
        error("Cannot proceed (Newton Failure) - increase Ri[1] or increase grid resolution")
    else
        println("Warning Newton RES:",Res,"- increase Ri[1] or increase grid resolution")
    end
    return x    
end
# Function to Make dx_right
function make_dx_right(R, d1, n)
    dx    = zeros(n)
    dx[1] = d1
    for i in 2:n
        dx[i] = dx[i - 1] * R
    end
    return dx
end
# Main Function
function define_new_grid(Ri,nr,Rfact,verbose)
    # Parameters
    if Rfact == 1.0         #Equally spaced grid
        x_left   = collect(LinRange(0.0, Ri[1], nr[1]))
        x_right  = collect(LinRange(Ri[1],Ri[2], nr[2]))
        dx_left  = diff(x_left)
        dx_right = diff(x_right)
    elseif Rfact > 1.0      #Refine both
        dx_left = collect(LinRange(1.0, 1.0 / Rfact, nr[1]))
        x_left = [0; cumsum(dx_left)] * Ri[1] / sum(dx_left)
        dx_left = diff(x_left)
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]
        d = dx_left[end]
        R = newton_solver(S, d, nr[2], 1e-8, 100,verbose)
        dx_right = make_dx_right(R, d, nr[2])
        x_right = [Ri[1]; Ri[1] .+ cumsum(dx_right)]
        dx_right = diff(x_right)
    else    #Refine only right grid
        Rfact   = abs(Rfact)
        x_left  = collect(LinRange(0.0, Ri[1], nr[1]))
        x_right = collect(LinRange(Ri[1], Ri[2], nr[2]))
        dx_left = diff(x_left)
        #Set Non-linear Problem
        S = Ri[2] - Ri[1]
        d = dx_left[end]
        R = newton_solver(S, d, nr[2], 1e-8, 100,verbose)
        dx_right = make_dx_right(R, d, nr[2])
        x_right  = [Ri[1]; Ri[1] .+ cumsum(dx_right)]
        dx_right = diff(x_right)
        if verbose == true
            println("Rfact: - $Rfact ; replaced by $Rfact on the right side, kept equal spacing on the left")
        end
    end
    Sc_left  = dx_left[1]/dx_left[end];
    Sc_right = dx_right[end]/dx_right[1];
    if verbose == true
        #Print results to check
        println("dx left  : ", dx_left[end])
        println("dx right : ", dx_right[1])
        println("x_left   : ", x_left[end])
        println("x_right  : ", x_right[end])
        println("Refinement factor for left is ",Sc_left)
        println("Refinement factor for right is ",Sc_right) 
        println("dx1/dx2 (left)  = ",dx_left[1]/dx_left[2])
        println("dx2/dx1 (right) = ",dx_right[2]/dx_right[1])
    end
    return Ri, nr, x_left, x_right, dx_left, dx_right, Sc_left, Sc_right
end
function main()
    # Run Main
    verbose =  true
    Ri      =  [0.5,   1.0]
    nr      =  [100,   100]
    Rfact   = -1.5;               #If negative, it uses R = 1 on the left, and abs(R) on the right
    Ri, nr, x_left, x_right, dx_left, dx_right , Sc_left, Sc_right = define_new_grid(Ri,nr,Rfact,verbose);

    # Plot results to check (for debugging purposes - uses Plots)
    p1 = plot(x_left, seriestype=:scatter, title="x_left")
    p2 = plot(dx_left, seriestype=:scatter, title="dx_left")
    p3 = plot(x_right, seriestype=:scatter, title="x_right")
    p4 = plot(dx_right, seriestype=:scatter, title="dx_right")
    plot(p1, p2, p3, p4, layout=(2, 2))
    display(plot!())
end
main()