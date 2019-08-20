push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using PyPlot
using ModStokesFastDirect
import AnalyticDomains
import PUX
include("../src/plottools.jl")

mss = ModifiedStokesSolver
unitcircle = AnalyticDomains.starfish(n_arms=3, amplitude=0.00)

function estimate_fourier_resolution(U)
    Uhat = fftshift(fft(U))
    Ures = max( norm(Uhat[1:2,:], Inf),
                norm(Uhat[end-1:end,:], Inf),
                norm(Uhat[:,1:2], Inf),
                norm(Uhat[:,end-1:end], Inf) ) / maximum(abs.(Uhat))
    return Ures
end


curve = AnalyticDomains.starfish()
Ngrid = 800
numpanels = 400
pux_R = 0.15

# Dev
# curve = unitcircle
# Ngrid = 100
# numpanels = 100

panelorder = 16
method = :dense_gmres
pux_ep = 2.0


plotext = false
ploterror=true


# Define ref problem (linear oscillatory)
k1 = 1.0
k2 = 1.0
ksq = k1^2 + k2^2
phi_k(x, y) = exp(1im*(k1*x + k2*y))
fhat1 = 1.0
fhat2 = -2.0
kdotfhat = k1*fhat1 + k2*fhat2
# rhs
ffunc1(x, y) = real(fhat1*phi_k(x,y))
ffunc2(x, y) = real(fhat2*phi_k(x,y))
# solution
ufunc1(x, y, alpha) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y) )
ufunc2(x, y, alpha) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y) )
# derivatives of solution
ufunc1x(x, y, alpha) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k1 )
ufunc1y(x, y, alpha) = real( (ksq*fhat1 - k1*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k2 )
ufunc2x(x, y, alpha) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k1 )
ufunc2y(x, y, alpha) = real( (ksq*fhat2 - k2*kdotfhat)/(ksq*(ksq+alpha^2))*phi_k(x,y)*1im*k2 )    

function solve(pux_R, alpha)
    info(">>>>>>>>>>>>>>>>>>>>>> alpha=$alpha")
    @show Ngrid
    
    # Discretize
    println("* Discretizing")
    dcurve = mss.discretize(curve, numpanels, panelorder, equal_arclength=true)

    if method==:fastdirect
        # Fast direct prep
        println("* Fast direct build")
        @time Ainv = compressed_system_matrix_inverse(dcurve, alpha, 64, fastdirect_TOL)
    elseif method==:fmm
        # GMRES prep
        println("* GMRES FMM prep")
        @time flhs = mss.system_matvec(dcurve, alpha)
    elseif method==:dense_lu
        println("* Dense matrix assembly")
        @time MAT = mss.system_matrix(dcurve, alpha)
        println("* Dense matrix LU")
        @time MATL, MATU, MATp = lu(MAT)
    elseif method==:dense_gmres
        println("* Dense matrix assembly")
        @time MAT = mss.system_matrix(dcurve, alpha)
    else
        error("unknown method")
    end
    
    println("* PUX precompute")
    @time PUXStor, X, Y, Lgrid, Ngrid_new, interior = PUX.pux_precompute(curve, dcurve, Ngrid, pux_ep, R=pux_R)
    xbdry = dcurve.points[1,:]
    ybdry = dcurve.points[2,:]
    interior_points = [X[interior]' ; Y[interior]']
    numeval = size(interior_points, 2)
    println("* FMM tree precompute")
    @time fmm_precomp = mss.doublelayer_precomp(dcurve, interior_points, alpha)
    println("* Near eval precompute")
    @time near_weights = mss.doublelayer_near_weights(dcurve, interior_points, alpha; ifgrad=true)
    hgrid = Lgrid/Ngrid_new
    
    # Setup test problem: BCs and RHS
    F1 = ffunc1.(X,Y)
    F2 = ffunc2.(X,Y)    
    ubc1 = ufunc1.(xbdry, ybdry, alpha)
    ubc2 = ufunc2.(xbdry, ybdry, alpha)
    
    # Extend RHS
    #println("* PUX eval")
    tic()
    FEXT1, FEXT2 = PUX.pux_eval_vec(F1, F2, PUXStor)
    time_pux = toq()
    # Compute particular solution
    #println("* Periodic FFT Solve")
    tic()
    U1, U2, upbdry1, upbdry2, dU1dx, dU1dy, dU2dx, dU2dy =
        mss.per_modstokes(FEXT1, FEXT2, Lgrid, alpha, xbdry, ybdry, ifgrad=true, zerok0=true)
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
        println("* GMRES dense solve")
        @time dsol, gmlog = IterativeSolvers.gmres(MAT, rhs; tol=1e-15, restart=100, maxiter=100,
                                             verbose=false, log=true)
        if !gmlog.isconverged
            warn("GMRES did not converge, residual = $(gmlog.data[:resnorm][end])")
        end
        density = reshape(dsol, 2, :)                        
    end
    time_solve = toq()
    # Estimate density resolution
    L = mss.legendre_matrix(panelorder)
    density_resolution = 0.0;
    for i=1:numpanels
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
    if plotext
        figure(4)
        clf()
        pcolor(X, Y, FEXT1)
        plot(xbdry, ybdry,"r")
        axis("equal")
    end    
    sol = mss.Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy,
                       copy(F1), copy(F2), copy(ubc1), copy(ubc2))
    return hgrid, errcheck(sol, X, Y, interior_points, interior, alpha)
end

function errcheck(sol, X, Y, interior_points, interior, alpha)
    ### END STUFF TO PUT IN TIMESTEP LOOP    
    # Check errors in test problem
    println("* Computing error")
    # Reference
    xint = interior_points[1,:]
    yint = interior_points[2,:]
    uref1 = ufunc1.(xint, yint, alpha)
    uref2 = ufunc2.(xint, yint, alpha)
    numeval = size(interior_points, 2)

    uref = [uref1 uref2]'

    maxnorm(v) = maximum(abs.(v))
    l2norm(v) = sqrt(sum(vec(v.^2)))
    
    # Error
    @show size(sol.U1)
    @show size(interior)
    u = [sol.U1[interior] sol.U2[interior]]'
    err = u .- uref
    maxerr = Array{Float64}(numeval)
    for i=1:numeval
        maxerr[i] = norm(err[:, i], Inf)
    end
    max_rel_err = maxnorm(err) / maxnorm(uref)
    l2_rel_err = l2norm(err) / l2norm(uref)
    println("Max rel err on grid: ", max_rel_err)
    println("L2 rel err on grid: ",  l2_rel_err)

    println("* Derivative errors (max rel)")
    grad = max(     norm(vec(sol.dU1dx[interior]), Inf),
                    norm(vec(sol.dU1dy[interior]), Inf),
                    norm(vec(sol.dU2dx[interior]), Inf),
                    norm(vec(sol.dU2dy[interior]), Inf) )
    

    gradnorm = max( norm(vec(sol.dU1dx[interior]), Inf),
                    norm(vec(sol.dU1dy[interior]), Inf),
                    norm(vec(sol.dU2dx[interior]), Inf),
    norm(vec(sol.dU2dy[interior]), Inf) )


    grad = [sol.dU1dx[interior]; sol.dU1dy[interior]; sol.dU2dx[interior]; sol.dU2dy[interior]]
    gradref = [ufunc1x.(xint, yint, alpha)
               ufunc1y.(xint, yint, alpha)
               ufunc2x.(xint, yint, alpha)
               ufunc2y.(xint, yint, alpha)]
    
    Egrad_max = maxnorm(grad.-gradref) / maxnorm(gradref)
    Egrad_l2 = l2norm(grad.-gradref) / l2norm(gradref)

    @show Egrad_max
    @show Egrad_l2    

    if ploterror
        println("* Displaying")
        figure(1, figsize=(4,3))
        clf()
        # Put on grid
        E = zeros(sol.U1)
        E[interior] = maxerr / maxnorm(uref)

        logE = log10.(E+1e-16)
        #interior_pcolor(X, Y, logE, interior, vmin=-16, vmax=0)
        interior_pcolor(X, Y, logE, interior, cmap=PyPlot.cm[:coolwarm])
        clim(log10(l2_rel_err)-1, log10(l2_rel_err))        
        #plot(dcurve.points[1,:], dcurve.points[2,:], ".k")
        axis("off")
        axis("equal")
        xlim(0.6, 1.3)
        ylim(-0.35, 0.35)
        cbar = colorbar(pad=0.0, shrink=0.8)
        cbar[:set_label]("Relative error")
        tight_layout(0.2)
        savefig(@sprintf("figs/starfish_errfield_alpha%.3f.png", alpha))
    end
    
    return max_rel_err, l2_rel_err, Egrad_max, Egrad_l2
end    


alpha_list = 10 .^ (-2.0:0.25:4.0)

Ntests = length(alpha_list)
max_rel_err = zeros(Ntests)
l2_rel_err = zeros(Ntests)
grad_max_rel_err = zeros(Ntests)
grad_l2_rel_err = zeros(Ntests)
hgrid = zeros(Ntests)


include("alpha_variation_plot.jl")

tic()
for i=1:Ntests
    hgrid[i], errors = solve(pux_R, alpha_list[i])
    max_rel_err[i], l2_rel_err[i], grad_max_rel_err[i], grad_l2_rel_err[i] = errors
    plot_results()
end
toq()


include("alpha_variation_plot.jl")
plot_results(paperplots=true)
savefig("../docs/paper/figs/alpha_variation.pdf")
