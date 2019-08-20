

using SpecialFunctions

include("spindown_besselzeros.jl")

function u_spindown(r, t, Re)
    # Assuming unit radius, speed and viscosity
    Omega, a = 1.0, 1.0
    nu = 1.0/Re
    S = 0.0
    ds = 0.0
    relerr = 0.0
    for n=1:length(lambda)
        ds = besselj(1, lambda[n]*r/a) / (lambda[n]*besselj(0, lambda[n]))*exp(-lambda[n]^2*nu*t/a)
        S += ds
        relerr = abs(ds/S) # If convergence is nice enough            
        if relerr < eps()
            break
        end
    end
    if relerr > eps()
        warn("Spindown solution not converged, relerr=", relerr)
    end
    u = -2*Omega*a*S
end


function Lu_spindown(r, t, Re)
    # Laplacian of the above
    Omega, a = 1.0, 1.0
    nu = 1.0/Re    
    S = 0.0
    ds = 0.0
    relerr = 0.0
    for n=1:length(lambda)
        b = lambda[n]/a
        ds = -b^2*besselj(1, b*r) / (lambda[n]*besselj(0, lambda[n]))*exp(-lambda[n]^2*nu*t/a)
        S += ds
        relerr = abs(ds/S) # If convergence is nice enough            
        if relerr < eps()
            break
        end
    end
    if relerr > eps()
        warn("Spindown solution not converged, relerr=", relerr)
    end
    u = -2*Omega*a*S
end

function Du_spindown(r, t, Re)
    # r derivative
    Omega, a = 1.0, 1.0
    nu = 1.0/Re    
    S = 0.0
    ds = 0.0
    relerr = 0.0
    for n=1:length(lambda)
        b = lambda[n]/a
        ds = b/2*(besselj(0, b*r)-besselj(2, b*r)) / (lambda[n]*besselj(0, lambda[n]))*exp(-lambda[n]^2*nu*t/a)
        S += ds
        relerr = abs(ds/S) # If convergence is nice enough            
        if relerr < eps()
            break
        end
    end
    if relerr > eps()
        warn("Spindown solution not converged, relerr=", relerr)
    end
    u = -2*Omega*a*S
end


function plot_spindown(t, Re)
    r = linspace(0,1,200)
    u = map(x -> u_spindown(x, t, Re), r)
    plot(r, u, label="Reference")
end
