import LinearMaps
import IterativeSolvers

function gmres_solve(dcurve, rhs, alpha;
                     maxiter=100, tol=1e-15,
                     verbose=true, flhs=nothing)
    if flhs==nothing
        flhs = system_matvec(dcurve, alpha)        
    end
    LHS = LinearMaps.LinearMap(flhs, 2*dcurve.numpoints)    
    # GMRES solve
    sol, gmlog = IterativeSolvers.gmres(LHS, rhs; tol=tol, restart=maxiter, maxiter=maxiter, verbose=false, log=true)
    msg = string("in ", gmlog.iters, " iterations, residual=", gmlog.data[:resnorm][end])
    if gmlog.isconverged
        if verbose
            println("GMRES converged $msg")
        end
    else
        warn("GMRES did not converge $msg")
    end
    density = reshape(sol, 2, dcurve.numpoints)
end
