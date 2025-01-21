using Plots
#Routine to test a linear diffusion code with Neumann BC
#The solution is based on the eigenfuction expansion method
function calc_sinus(nx,L,D,tot,nterms)
    x   = collect(LinRange(0,L,nx));
    C   = ones(nx);
    Cin = ones(nx);
    C1   = 1.0; 
    C1in = 1.0;
    for n = 1:nterms
        C   .= C   .+ 2*L/pi./x .* ((-1).^n)./n .* sin.(n .* pi .*x./L) .* exp(-D*n^2*pi^2*tot/L^2);
        Cin .= Cin .+ 2*L/pi./x .* ((-1).^n)./n .* sin.(n .* pi .*x./L) .* exp(-D*n^2*pi^2*1e-5/L^2);
        C1   = C1   + 2*((-1).^n).* exp(-D*n^2*pi^2*tot/L^2)
    end
    C[1]   = C1;
    Cin[1] = 0.0;
    return x,C,Cin
end
function fd_num(x,Cin,D,tot)
    dx   = x[2]-x[1]
    dt   = dx^2/D*0.3
    nx   = length(Cin)
    Cnum = zeros(nx);   #Actual Initial
    Cnum[end] = 1.0;
    xc        = 0.5*x[1:end-1] .+ 0.5*x[2:end];
    t         = 0.0;
    while t<tot
        t = t + dt
        if t>tot
            t  = tot
            dt = dt -(t-tot)
        end

        Cnum[2:nx-1] = Cnum[2:nx-1] + dt/dx/dx*D*diff(xc.^2 .*diff(Cnum))./x[2:nx-1].^2
        Cnum[1]      = Cnum[2]
        Cnum[nx]     = 1.0;
    end
    Cnum[1] = Cnum[2];
    return Cnum 
end
function main()
    D     = 1.1
    L     = 1.0
    tot   = 1e-1
    nx    = 200
    x,C,Cin = calc_sinus(nx,L,D,tot,nx*5)
    Cnum    = fd_num(x,Cin,D,tot)
    return x, Cin, C, Cnum
end
x, Cin, C, Cnum = main()
plot(x, C)
plot!(x,Cnum,color=:black,linestyle=:dash)
plot!(x, Cin)