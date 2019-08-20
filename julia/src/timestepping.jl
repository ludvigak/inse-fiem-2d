# Routines for taking steps using either Euler or SDC

using FastGaussQuadrature

function step_euler(dt, sol, solver, ubc1, ubc2)
    alpha = solver.alpha
    F1, F2 = ugradu(sol)
    RHS1 = alpha^2*(sol.U1 - dt*F1)
    RHS2 = alpha^2*(sol.U2 - dt*F2)
    newsol = solver.solve(RHS1, RHS2, ubc1, ubc2)
    return newsol
end


function step_bdf(dt, prevsol, sol, solver, ubc1, ubc2)
    alpha = solver.alpha
    F1, F2 = ugradu(sol)
    prevF1, prevF2 = ugradu(prevsol)
    RHS1 = 1/3*alpha^2*(4*sol.U1 - prevsol.U1 - 2*dt*(2*F1 - prevF1))
    RHS2 = 1/3*alpha^2*(4*sol.U2 - prevsol.U2 - 2*dt*(2*F2 - prevF2))
    newsol = solver.solve(RHS1, RHS2, ubc1, ubc2)
    return newsol
end

function sdc_substep(dt, sdc_order)
    return sdc_order > 1 ? dt / (sdc_order-1) : dt
end

function sdc_integration(M)
    # Setup substeps
    x = collect(linspace(-1, 1, M))
    # Vandermonde matrix (transposed)
    At = ones(M,M)
    for m=2:M
        @. At[m, :] = At[m-1, :]*x
    end
    # Use Gauss-Legendre for quadrature
    xleg, wleg = gausslegendre(M)
    # Compute substep integration weights
    I = zeros(M,M-1)
    for m=1:M-1
        a, b = x[m], x[m+1]
        wint = wleg*(b-a)/2
        xint = (xleg+1)*(b-a)/2 + a
        xint_i = ones(M)
        p = zeros(M)
        p[1] = sum(wint)
        for i=2:M
            xint_i .*= xint            
            p[i] = sum(wint .* xint_i)
        end
        I[:,m] = At\p
    end
    return I/2
end

function step_sdc(dt, order, sol, solver, ubc1, ubc2)
    # TODO: Need to implement BC time dependence
    alpha = solver.alpha
    Re = alpha^2*dt
    P = order # order-1 can also work
    K = order
    if order>1
        W = sdc_integration(P)
    else
        # Fallback to Euler
        P = 2
        K = 1
    end
    # Set initial conditions for all passes
    solutions = Array{Solution}(P, K)    
    for k=1:K
        solutions[1, k] = sol
    end    
    # First pass (k=0)
    # Eulers steps
    for m=1:P-1
        init = solutions[m, 1]
        F1init, F2init = ugradu(init)    
        RHS1 = alpha^2*init.U1 - Re*F1init
        RHS2 = alpha^2*init.U2 - Re*F2init
        solutions[m+1, 1] = solver.solve(RHS1, RHS2, ubc1, ubc2)        
    end

    # Subsequent passes
    for k=1:K-1        
        for m=1:P-1
            # Load historic stages that we need
            init = solutions[m,   k+1] # Substep initial condition
            pred = solutions[m+1, k  ] # Prediction
            hist = solutions[m  , k  ] # Prediction history
            # Post-process stages
            F1init, F2init = ugradu(init)
            #
            L1pred = alpha^2*pred.U1 - pred.RHS1
            L2pred = alpha^2*pred.U2 - pred.RHS2
            #           
            F1hist, F2hist = ugradu(hist)    
            # Integrate
            s = size(init.U1)
            I1, I2 = zeros(s), zeros(s)
            L = dt*(order-1)
            w = W[:, m]*L
            for n=1:P
                tmp = solutions[n, k]
                F1tmp, F2tmp = ugradu(tmp)
                L1tmp = alpha^2*tmp.U1 - tmp.RHS1
                L2tmp = alpha^2*tmp.U2 - tmp.RHS2
                I1 += w[n]*( L1tmp - Re*F1tmp )
                I2 += w[n]*( L2tmp - Re*F2tmp )
            end

            # f(t) = 1 + t
            # @show t1 = quadgk(f, dt*(m-1), dt*m)[1]
            # @show t2 = sum(w .* f(0:dt:dt*(P-1)))
            # @show t1-t2
            # @assert abs(t1-t2)/abs(t1) < 1e-14
            
            # Create new RHS
            RHS1 = alpha^2*init.U1 - Re*(F1init-F1hist) - L1pred + I1/dt
            RHS2 = alpha^2*init.U2 - Re*(F2init-F2hist) - L2pred + I2/dt
            #            
            solutions[m+1, k+1] = solver.solve(RHS1, RHS2, ubc1, ubc2) 
        end
    end

    if K>1
        # Output error est (UNTESTED)
        maxnorm(v) = maximum(abs.(v))
        l2norm(v) = sqrt(sum(vec(v.^2)))
        corr1 = solutions[P, K].U1 - solutions[P, K-1].U1
        corr2 = solutions[P, K].U2 - solutions[P, K-1].U2
        sol = solutions[P,K]        
        info("U1 SDC err est: l2=", l2norm(corr1)/l2norm(sol.U1),", max=", maxnorm(corr1)/maxnorm(sol.U1))
        info("U2 SDC err est: l2=", l2norm(corr2)/l2norm(sol.U2),", max=", maxnorm(corr2)/maxnorm(sol.U2))
    end    
    return solutions[P,K]

    # Original order 2 code
    # # First take Euler step
    # F1, F2 = ugradu(sol)    
    # RHS1 = alpha^2*U1 - Re*F1
    # RHS2 = alpha^2*U2 - Re*F2
    # prov = solver.solve(RHS1, RHS2, ubc1, ubc2)
    # # Post process Euler step
    # F1prov, F2prov = ugradu(prov)
    # L1prov = alpha^2*prov.U1 - prov.RHS1
    # L2prov = alpha^2*prov.U2 - prov.RHS2
    # # Post process initial conds
    # L1 = alpha^2*sol.U1 - sol.RHS1
    # L2 = alpha^2*sol.U2 - sol.RHS2
    # # Integrate
    # I1 = (L1-Re*F1 + L1prov-Re*F1prov)/2
    # I2 = (L2-Re*F2 + L2prov-Re*F2prov)/2
    # # Take correction step
    # RHS1 = Re/dt*U1 - L1prov + Re*I1
    # RHS2 = Re/dt*U2 - L2prov + Re*I2
    # newsol = solver.solve(RHS1, RHS2, ubc1, ubc2)        
    # return newsol    
end
