# Show near evaluation inside starfish

push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using PyPlot
import IterativeSolvers
import LinearMaps
import AnalyticDomains
using ModifiedStokesSolver

include("../src/plottools.jl")

# Params
alpha = 5.0
numpanels = 100
panelorder = 16
solve_inteq = true
# Discretize
curve = AnalyticDomains.starfish(n_arms=5, amplitude=0.3)
println("* Discretizing")
@time dcurve = discretize(curve, numpanels, panelorder,
                          equal_arclength=true)
N = dcurve.numpoints

if solve_inteq
    # Setup a problem
    xsrc = [5.0, 5.0]
    srand(1)
    qsrc = 1-2*rand(2)
    nsrc = 1-2*rand(2)
    fsrc = 1-2*rand(2)
    ubc(x) = stresslet(x-xsrc, qsrc, nsrc, alpha) +
        stokeslet(x-xsrc, fsrc, alpha)
    # Right hand side
    rhs = zeros(2*N)
    for i=1:dcurve.numpoints
        rhs[2(i-1)+(1:2)] = ubc(dcurve.points[:,i])
    end

    # FMM lhs
    @time flhs = system_matvec(dcurve, alpha)
    LHS = LinearMaps.LinearMap(flhs, 2*N)
    
    # GMRES solve
    println("* GMRES solve")
    maxiter = 50
    tol = 1e-15
    @time sol, gmlog = IterativeSolvers.gmres(LHS, rhs; tol=tol, restart=maxiter, maxiter=maxiter, verbose=false, log=true)
    println( (gmlog.isconverged ? "  Converged" : "Did NOT converge"),
             " in $(gmlog.iters) iterations, residual=",
             gmlog.data[:resnorm][end])
    density = reshape(sol, 2, N)
    gmlog.isconverged
else
    # Just compare to a unit density on a refined curve
    println("* Discretizing ref curve")    
    refmul = 5
    @time refgrid = discretize(curve, numpanels, refmul*panelorder, equal_arclength = true)    
    density = ones(size(dcurve.points))
    refden = ones(size(refgrid.points))
    ubc(zt) = doublelayer(refgrid, refden, zt, alpha; specquad=false)  
end

## POST PROCESS
xmax = maximum(dcurve.points[1,:])
xmin = minimum(dcurve.points[1,:])
ymax = maximum(dcurve.points[2,:])
ymin = minimum(dcurve.points[2,:])

xmin = 0.0
ymin = 0.0
xmin = 0.4
ymax = 0.6

N = 200
# Points sent to pcolor
xplot = linspace(xmin, xmax, N)
yplot = linspace(ymin, ymax, N)
# Create ndgrid using centroids
xcen = (xplot[1:end-1] + xplot[2:end])/2
ycen = (yplot[1:end-1] + yplot[2:end])/2
X, Y = ndgrid(xcen, ycen)
# Vector of grid points
zt = [vec(X) vec(Y)]'
# Get interior point mask
interior = interior_points(dcurve, zt)
numeval = sum(interior)
zt_interior = zt[:, interior]

println("* Computing near weights")
@time weights = ModifiedStokesSolver.doublelayer_near_weights(dcurve, zt_interior, alpha)
@time weights = ModifiedStokesSolver.doublelayer_near_weights(dcurve, zt_interior, alpha)

println("* Computing field on grid, using precomputed weights")
@time u = doublelayer_fast(dcurve, density, zt_interior, alpha; weights=weights)
@time u = doublelayer_fast(dcurve, density, zt_interior, alpha; weights=weights)

println("* Computing field on grid, on the fly")
@time uotf = doublelayer_fast(dcurve, density, zt_interior, alpha)
@time uotf = doublelayer_fast(dcurve, density, zt_interior, alpha)

println("* Computing field on grid, without specquad")
@time udummy = doublelayer_fast(dcurve, density, zt_interior, alpha; specquad=false)
@time udummy = doublelayer_fast(dcurve, density, zt_interior, alpha; specquad=false)

println("* Difference between on the fly and precomputed weights:")
println("  ", norm(vec(u)-vec(uotf),Inf))

println("* Computing error")
# Reference
uref = Array{Float64}(2, numeval)    
for i=1:numeval
    uref[:,i] = ubc(zt_interior[:, i])
end
unorm = norm(vec(uref), Inf)
# Error
err = u .- uref
maxerr = Array{Float64}(numeval)
for i=1:numeval
    maxerr[i] = norm(err[:, i], Inf)
end
println("Maxerr on grid: ", norm(vec(maxerr), Inf))
# Put on grid
E = zeros(N-1,N-1)
E[interior] = maxerr
Erel = E / unorm
logE = log10.(Erel+1e-100)



println("* Displaying")
clf()

interior_pcolor(X, Y, logE, interior, vmin=-16, vmax=0)

cbar = colorbar()
cbar[:set_label]("Rel. error")
plot(dcurve.points[1,:], dcurve.points[2,:], ".-k")
axis("image")
axis([xmin, xmax, ymin, ymax])
show()
