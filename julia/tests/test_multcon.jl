push!(LOAD_PATH, string(pwd(),"/src"))

using ModifiedStokesSolver
import AnalyticDomains
using Base.Test
using PyPlot
import IterativeSolvers


@testset "MultiplyConnected" begin

    # Setup discretization
    R1 = 1.0
    R2 = 0.25
    alpha = 3.0
    panelorder = 16    
    numpanels1 = 40
    numpanels2 = 8
    tol = 1e-8
                                              
    # Setup problem
    xsrc1 = [0.01, 0.02]
    xsrc2 = [1.0, 1.0]
    srand(1)
    qsrc = 1-2*rand(2)
    nsrc = 1-2*rand(2)
    fsrc = 1-2*rand(2)
    ubc(x) = stresslet(x-xsrc1, qsrc, nsrc, alpha) +
        stokeslet(x-xsrc1, fsrc, alpha) +
        stresslet(x-xsrc2, qsrc, nsrc, alpha) +
        stokeslet(x-xsrc2, fsrc, alpha)

    # Discretize
    curve1 = AnalyticDomains.starfish(;amplitude = 0.1, radius = R1)
    curve2 = AnalyticDomains.starfish(;amplitude = 0.0, radius = R2, interior=true)
    dcurve = discretize([curve1, curve2], [numpanels1, numpanels2], panelorder)
    N = dcurve.numpoints

    # Right hand side
    rhs = zeros(2*N)
    for i=1:dcurve.numpoints
        rhs[2(i-1)+(1:2)] = ubc(dcurve.points[:,i])
    end

    # Dense matrix 
    println("* Matrix assembly and solve")
    @time LHS = system_matrix(dcurve, alpha)
    @time sol = LHS\rhs
    density = reshape(sol, 2, N)

    # Interior test point
    zt = [(R1+R2)/2, 0.01]
    ref = ubc(zt)
    dlp = doublelayer(dcurve, density, zt, alpha)
    relerr = maximum(abs.(ref.-dlp)) / maximum(abs.(ref))
    println("* Interior eval")
    @show relerr
    @test relerr < 1e-15
    
    # Field test
    Ngrid = 50
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
    @time u = doublelayer_fast(dcurve, density, zt_interior, alpha, specquad=true)
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

    max_relerr_grid = norm(vec(maxerr), Inf) / unorm
    @show max_relerr_grid
    @test max_relerr_grid < 1e-13

    # Plot
    E = zeros(X)
    E[interior] = maxerr
    clf()
    pcolor(X, Y, log10.(E))
    colorbar()
end
