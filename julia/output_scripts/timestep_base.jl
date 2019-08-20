# Base code for timestepping, precomputes a solver and tests it on a simple problem
push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using PyPlot
using ModStokesFastDirect
using Humanize
import ModifiedStokesSolver
import AnalyticDomains
import PUX

mss = ModifiedStokesSolver
include("../src/plottools.jl")
include("../src/timestepping.jl")

Solution = mss.Solution
PrecomputedSolver = mss.PrecomputedSolver
unitcircle = AnalyticDomains.starfish(n_arms=0, amplitude=0.0)


function estimate_fourier_resolution(U)
    Uhat = fftshift(fft(U))
    Ures = max( norm(Uhat[1:2,:], Inf),
                norm(Uhat[end-1:end,:], Inf),
                norm(Uhat[:,1:2], Inf),
                norm(Uhat[:,end-1:end], Inf) ) / maximum(abs.(Uhat))
    return Ures
end

function print_footprint(var)
    println("Footprint: ", datasize(Base.summarysize(var)))
end

# Precompute a solver
function precompute_solver(;curve=unitcircle, alpha=2.0,
                           numpanels=100, panelorder=16,
                           method=:heuristic, fastdirect_TOL=1e-10,
                           Ngrid=150, pux_ep=2, pux_R=0.4,
                           dotest=true)
    println("======= SETUP")
    @show alpha
    # Discretize
    println("* Discretizing")
    dcurve = mss.discretize(curve, numpanels, panelorder, equal_arclength=true)
    clf()
    plot_boundary(dcurve)
    # PUX prep
    println("* PUX precompute")
    @time PUXStor, X, Y, Lgrid, Ngrid, interior = PUX.pux_precompute(curve, dcurve, Ngrid, pux_ep, R=pux_R)
    print_footprint(PUXStor)
    figure()
    plot(dcurve.points[1,:], dcurve.points[2,:], ".k")
    for s=PUXStor
        PUX.plot_partitions(s)
    end
    axis("image")
    if method==:heuristic
        # Heuristic for small alpha and smooth problem
        # Other problems may change fastdirect build speed and GMRES convergence
        if sum(numpanels) <= 500
            method = :dense_gmres
        else
            method = :fastdirect
        end
        println("* Heuristic solver choice: $method")
    end

    if method==:fastdirect
        newmethod = :dense_gmres
        warn("Fast direct solver not available, switching to $newmethod")
        method = newmethod
    end
    
    if method==:fastdirect
        # Fast direct prep
        println("* Fast direct build")
        @time Ainv = compressed_system_matrix_inverse(dcurve, alpha, 64, fastdirect_TOL)
        print_footprint(Ainv)
    elseif method==:fmm
        # GMRES prep
        println("* GMRES FMM prep")
        @time flhs = mss.system_matvec(dcurve, alpha)
        print_footprint(flhs)
    elseif method==:dense_lu
        println("* Dense matrix assembly")
        @time MAT = mss.system_matrix(dcurve, alpha)
        print_footprint(MAT)       
        println("* Dense matrix LU")
        @time MATL, MATU, MATp = lu(MAT)
    elseif method==:dense_gmres
        println("* Dense matrix assembly")
        @time MAT = mss.system_matrix(dcurve, alpha)
        print_footprint(MAT)
    elseif method==:none
        println("* Dummy solver")
        dotest = false
    else
        error("unknown method")
    end    
    # Useful vars
    xbdry = dcurve.points[1,:]
    ybdry = dcurve.points[2,:]        
    interior_points = [X[interior]' ; Y[interior]']
    numeval = size(interior_points, 2)
    if method == :none
        solve = nothing
    else
        println("* FMM tree precompute")
        @time fmm_precomp = mss.doublelayer_precomp(dcurve, interior_points, alpha)
        print_footprint(fmm_precomp)       
        println("* Near eval precompute")
        @time near_weights = mss.doublelayer_near_weights(dcurve, interior_points, alpha; ifgrad=true)
        print_footprint(near_weights)
        # Form function that does the solve required in a timestep
        function solve(F1, F2, ubc1, ubc2)
            # Extend RHS
            #println("* PUX eval")
            tic()
            FEXT1, FEXT2 = PUX.pux_eval_vec(F1, F2, PUXStor)
            time_pux = toq()
            # Compute particular solution
            #println("* Periodic FFT Solve")
            tic()
            U1, U2, upbdry1, upbdry2, dU1dx, dU1dy, dU2dx, dU2dy =
                mss.per_modstokes(FEXT1, FEXT2, Lgrid, alpha, xbdry, ybdry, ifgrad=true)
            time_fft = toq()
            # Estimate particular solution resolution
            UP_resolution =  max(estimate_fourier_resolution(U1), 
                                 estimate_fourier_resolution(U2))
            UPgrad_resolution = max( estimate_fourier_resolution(dU1dx),
                                     estimate_fourier_resolution(dU1dy),
                                     estimate_fourier_resolution(dU2dx),
                                     estimate_fourier_resolution(dU2dy) )
            # Integral equation RHS
            rhs1 = ubc1 - upbdry1
            rhs2 = ubc2 - upbdry2
            rhs = vec([rhs1'; rhs2'])
            # Solve integral equation
            tic()
            if method==:fastdirect
                #println("* Fast direct solve")
                dsol = apply_compressed_matrix(Ainv, rhs)
                density = reshape(dsol, 2, :)
            elseif method==:fmm
                #println("* GMRES FMM solve")
                density = mss.gmres_solve(dcurve, rhs, alpha, flhs=flhs)
            elseif method==:dense_lu
                #println("* LU apply")
                dsol = MATU\(MATL\rhs[MATp])
                density = reshape(dsol, 2, :)            
            elseif method==:dense_gmres
                #println("* GMRES dense solve")
                dsol, gmlog = IterativeSolvers.gmres(MAT, rhs; tol=1e-15, restart=100, maxiter=100,
                verbose=false, log=true)
                if !gmlog.isconverged
                    warn("GMRES did not converge, residual = $(gmlog.data[:resnorm][end])")
                end
                density = reshape(dsol, 2, :)                        
            end
            time_solve = toq()
            # Estimate density resolution
            L = mss.legendre_matrix(panelorder)
            density_resolution = 0.0
            for i=1:dcurve.numpanels
                xcoeff = L*density[1, (1:panelorder) + panelorder*(i-1)]
                ycoeff = L*density[2, (1:panelorder) + panelorder*(i-1)]
                xres = max(abs(xcoeff[end]), abs(xcoeff[end-1]))
                yres = max(abs(ycoeff[end]), abs(ycoeff[end-1]))
                res = max(xres, yres)
                density_resolution = max(density_resolution, res)
            end        
            # Evaluate layer potential at interior points
            #println("* FMM eval")
            tic()
            u, u1grad, u2grad = mss.doublelayer_fast(dcurve, density, interior_points, alpha; specquad=true, weights=near_weights, ifgrad=true, precomp=fmm_precomp)
            time_fmm = toq()
            @. begin
                # Add hom sol to part sol on grid            
                U1[interior]   += u[1,:]     
                U2[interior]   += u[2,:]     
                dU1dx[interior] += u1grad[1,:]
                dU1dy[interior] += u1grad[2,:]
                dU2dx[interior] += u2grad[1,:]
                dU2dy[interior] += u2grad[2,:]
                # Clear exterior
                U1[!interior] = 0.0
                U2[!interior] = 0.0
                dU1dx[!interior] = 0.0
                dU1dy[!interior] = 0.0
                dU2dx[!interior] = 0.0
                dU2dy[!interior] = 0.0            
            end
            println(@sprintf("Timings: pux=%.3f, fft=%.3f, solve=%.3f, fmm=%.3f", time_pux, time_fft, time_solve, time_fmm))
            info(@sprintf("Resolution est.: den=%.1e, uP=%.1e, uP'=%.1e",
                          density_resolution, UP_resolution, UPgrad_resolution))        
            return Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy,
            copy(F1), copy(F2), copy(ubc1), copy(ubc2))
        end
    end
    solver = PrecomputedSolver(alpha, dcurve, X, Y, Lgrid, PUXStor, interior, solve)
    if dotest
        test_solver(solver)
    end
    return solver
end

# Useful functions
function pux_gradient(U1, U2, PUXStor, Lgrid)
    U1EXT, U2EXT = PUX.pux_eval_vec(U1, U2, PUXStor)    
    dU1dx, dU1dy = mss.per_gradient(U1EXT, Lgrid)
    dU2dx, dU2dy = mss.per_gradient(U2EXT, Lgrid)
    return dU1dx, dU1dy, dU2dx, dU2dy
end
function ugradu(sol::Solution)
    ugradu(sol.U1, sol.U2, sol.dU1dx, sol.dU1dy, sol.dU2dx, sol.dU2dy)
end
function ugradu(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy)
    F1 = U1.*dU1dx + U2.*dU1dy
    F2 = U1.*dU2dx + U2.*dU2dy
    return F1, F2
end


# Simple tester for precomputed solver
function test_solver(precomp)    
    # Unpack
    dcurve = precomp.dcurve
    alpha = precomp.alpha
    X = precomp.X
    Y = precomp.Y
    interior = precomp.interior
    xbdry = dcurve.points[1,:]
    ybdry = dcurve.points[2,:]        
    interior_points = [X[interior]' ; Y[interior]']
    numeval = size(interior_points, 2)
    # Define ref problem (linear oscillatory)
    k1 = 1.0
    k2 = 2.0
    ksq = k1^2 + k2^2
    phi_k(x, y) = exp(1im*(k1*x + k2*y))
    fhat1 = 3.0
    fhat2 = 5.0
    kdotfhat = k1*fhat1 + k2*fhat2
    # rhs
    ffunc1(x, y) = real(fhat1*phi_k(x,y))
    ffunc2(x, y) = real(fhat2*phi_k(x,y))
    # solution
    ufunc1(x, y) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y) )
    ufunc2(x, y) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y) )
    # derivatives of solution
    ufunc1x(x, y) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k1 )
    ufunc1y(x, y) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k2 )
    ufunc2x(x, y) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k1 )
    ufunc2y(x, y) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k2 )    
    
    # Setup test problem: BCs and RHS
    F1 = ffunc1.(X,Y)
    F2 = ffunc2.(X,Y)    
    ubc1 = ufunc1.(xbdry, ybdry)
    ubc2 = ufunc2.(xbdry, ybdry)
    
    ### BEGIN STUFF TO PUT IN TIMESTEP LOOP
    println("======= TEST SOLVE")
    @time begin        
        # Solve inhom mod Stokes
        sol = precomp.solve(F1, F2, ubc1, ubc2)
        newF1 = sol.U1.*sol.dU1dx + sol.U2.*sol.dU1dy
        newF2 = sol.U1.*sol.dU2dx + sol.U2.*sol.dU2dy
        println("======= TEST SOLVE DONE")
    end
    ### END STUFF TO PUT IN TIMESTEP LOOP    
    # Check errors in test problem
    println("* Computing error")
    # Reference
    xint = interior_points[1,:]
    yint = interior_points[2,:]
    uref1 = ufunc1.(xint, yint)
    uref2 = ufunc2.(xint, yint)

    uref = [uref1 uref2]'
    unorm = norm(vec(uref), Inf)
    # Error
    u = [sol.U1[interior] sol.U2[interior]]'
    err = u .- uref
    maxnorm(v) = maximum(abs.(v))
    l2norm(v) = sqrt(sum(vec(v.^2)))
    println("Max rel err on grid: ", maxnorm(err)/maxnorm(uref))
    println("L2 rel err on grid: ", l2norm(err)/l2norm(uref) )    

    println("* Derivative errors (max rel)")
    gradnorm = max( norm(sol.dU1dx[interior], Inf),
                    norm(sol.dU1dy[interior], Inf),
                    norm(sol.dU2dx[interior], Inf),
                    norm(sol.dU2dy[interior], Inf) )
    E_du1dx = norm(ufunc1x.(xint, yint) - sol.dU1dx[interior], Inf) / gradnorm
    E_du1dy = norm(ufunc1y.(xint, yint) - sol.dU1dy[interior], Inf) / gradnorm
    E_du2dx = norm(ufunc2x.(xint, yint) - sol.dU2dx[interior], Inf) / gradnorm
    E_du2dy = norm(ufunc2y.(xint, yint) - sol.dU2dy[interior], Inf) / gradnorm
    Egrad_max = max(E_du1dx, E_du1dy, E_du2dx, E_du2dy)
    @show Egrad_max

    println("* Displaying")
    clf()
    # Put on grid
    E = zeros(X)
    # Max error
    maxerr = maximum(abs.(err), 1)
    E[interior] = maxerr / unorm
    # Or vorticity
    vort_comp = sol.dU2dx[interior] - sol.dU1dy[interior]
    vort_ref = ufunc2x.(xint, yint) - ufunc1y.(xint, yint)      
    E[interior] = abs.(vort_comp - vort_ref) / maxnorm(vort_ref)
    title("Vorticity error")

    logE = log10.(E+1e-100)
    #interior_pcolor(X, Y, logE, interior, vmin=-16, vmax=0)
    interior_pcolor(X, Y, logE, interior, vmax=0, cmap=PyPlot.cm[:coolwarm])
    cbar = colorbar()
    cbar[:set_label]("Rel. error")
    plot(dcurve.points[1,:], dcurve.points[2,:], ".k")
    for s=precomp.PUXStor
        PUX.plot_partitions(s)
    end
    
    axis("image")
    xmax = maximum(dcurve.points[1,:])*1.1
    xmin = minimum(dcurve.points[1,:])*1.1
    ymax = maximum(dcurve.points[2,:])*1.1
    ymin = minimum(dcurve.points[2,:])*1.1
    axis([xmin, xmax, ymin, ymax])
end

function timestep(Nstep, dt, sdc_order, hook, solver, t0, sol, ubc1, ubc2; istart=0)
    dt_substep = sdc_substep(dt, sdc_order)
    Re = solver.alpha^2*dt_substep
    ti = t0 + istart*dt
    for i=istart+1:Nstep
        println("=== STEP ", i, "/", Nstep)
        @time begin
            ti = t0 + i*dt
            newsol = step_sdc(dt_substep, sdc_order, sol, solver, ubc1, ubc2)    
            prevsol, sol = sol, newsol
            hook(solver, sol, ti, i)
        end
    end
    return ti, sol
end

hook_null(a1, a2, a3, a4) = nothing


;
