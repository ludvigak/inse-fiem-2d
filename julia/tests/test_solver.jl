# Tests solving a full problem and evaluating solution near boundary

push!(LOAD_PATH, string(pwd(),"/src"))

import IterativeSolvers
using ModifiedStokesSolver
import AnalyticDomains
using Base.Test

@testset "IntegralEquation" begin
    # Run over large alpha span
    for alpha=[1.0, 10.0, 100.0, 1000.0]
        @testset "alpha=$alpha" begin
            # Discretize
            numpanels = 22
            panelorder = 16
            curve = AnalyticDomains.starfish(n_arms=3, amplitude=0.1)
            dcurve = discretize(curve, numpanels, panelorder)
            N = dcurve.numpoints
            
            # Setup a problem
            xsrc = [2.0, 2.0]
            srand(1)
            qsrc = 1-2*rand(2)
            nsrc = 1-2*rand(2)
            fsrc = 1-2*rand(2)
            ubc(x) = stresslet(x-xsrc, qsrc, nsrc, alpha) +
                stokeslet(x-xsrc, fsrc, alpha)

            src = [2.0 2.0]'
            str = [-3.0 5.0]'
            ureffunc(targ) = fmm_stokeslet_targ(src, targ, str, alpha)
            urefgradfunc(targ) = fmm_stokeslet_targ(src, targ, str, alpha; ifgrad=true)    
            ubc(x) = ureffunc(hcat(x,x))[:,1]

            
            # Right hand side
            rhs = zeros(2*N)
            for i=1:dcurve.numpoints
                rhs[2(i-1)+(1:2)] = ubc(dcurve.points[:,i])
            end

            # Dense matrix lhs
            LHS = system_matrix(dcurve, alpha)

            # GMRES solve
            maxiter = 50
            tol = 1e-15
            sol, gmlog = IterativeSolvers.gmres(LHS, rhs; tol=tol, restart=maxiter, maxiter=maxiter, verbose=false, log=true)
            density = reshape(sol, 2, N)
            @test gmlog.isconverged
            
            ## POST PROCESS
            # Interior point
            zt = [0.1, 0.2]
            ref = ubc(zt)
            dlp = doublelayer(dcurve, density, zt, alpha)
            relerr = maximum(abs.(ref.-dlp)) / maximum(abs.(ref))
            @test relerr < 1e-14
            # Near boundary
            @testset "Near evaluation" begin    
                # Single point
                t1 = dcurve.t_edges[:,1]        
                h1 = t1[2]-t1[1]
                #ztmp = curve.tau(t1[1]+h1/2 + 1im*h1*0.1)
                # Test points very close to panel edge and panel nodes
                ztmp = curve.tau(t1[2]+1e-5 + 1im*h1*1e-4)
                zclose = [real(ztmp); imag(ztmp)]
                zclose = hcat(zclose, dcurve.points[:,3] + 1e-4*dcurve.normals[:,3])
                ref, grad1ref, grad2ref = urefgradfunc(zclose)
                dlp, g1, g2 = doublelayer_fast(dcurve, density, zclose, alpha, ifgrad=true)
                err = ref .- dlp
                relerr_near = maximum(abs.(err)) / maximum(abs.(ref))
                @test relerr_near < 2e-12
                g = [g1 g2]
                gradref = [grad1ref grad2ref]
                grad_relerr_near = maximum(abs.(g-gradref)) / maximum(abs.(gradref))
                @test grad_relerr_near < 2e-10
                # Grid
                Ngrid = 20
                xmax = maximum(dcurve.points[1,:])
                xmin = minimum(dcurve.points[1,:])
                ymax = maximum(dcurve.points[2,:])
                ymin = minimum(dcurve.points[2,:])
                x = linspace(xmin, xmax, Ngrid)
                y = linspace(ymin, ymax, Ngrid)
                X, Y = ndgrid(x, y)
                zt = [vec(X) vec(Y)]'

                interior = interior_points(dcurve, zt)
                numeval = sum(interior)
                zt_interior = zt[:, interior]
                # Compute on grid
                u, u1grad, u2grad = doublelayer_fast(dcurve, density, zt_interior, alpha, ifgrad=true)
                # Reference
                uref, grad1ref, grad2ref = urefgradfunc(zt_interior)
                unorm = norm(vec(uref), Inf)
                # Error
                err = u .- uref
                maxerr = Array{Float64}(numeval)
                for i=1:numeval
                    maxerr[i] = norm(err[:, i], Inf)
                end
                max_relerr_grid = norm(vec(maxerr), Inf) / unorm
                @test max_relerr_grid < 4e-14

                refgrad = hcat(grad1ref, grad2ref)
                grad = hcat(u1grad, u2grad)
                max_grad_relerr_grid = maximum(abs.(refgrad-grad)) / maximum(abs.(refgrad)) 
                @test max_grad_relerr_grid < 7e-12
            end
        end
    end
end
