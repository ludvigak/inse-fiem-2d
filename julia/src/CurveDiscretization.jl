__precompile__()
module CurveDiscretization

using Compat
using Compat.Printf

using FastGaussQuadrature
using FMMLIB2D

export DiscreteCurve
export discretize
export multiply_connected
export traparclen
export minmax_panel
export interior_points
export nearest_nb_R

# Curve discretization
struct DiscreteCurve
    panelorder::Int64
    numpanels::Int64
    numpoints::Int64
    edges::Array{Float64, 2}
    n_edges::Array{Float64, 2}        
    t_edges::Array{Float64, 2}
    points::Array{Float64, 2}
    normals::Array{Float64, 2}
    weights::Array{Float64}
    dS::Array{Float64}
    curvature::Array{Float64}
    prevpanel::Array{Int64}
    nextpanel::Array{Int64}
    curvenum::Array{Int64}    
end

include("legendre.jl")
include("mappings.jl")

"""
    Discretize curve using composite Gauss-Legendre quadrature
    """
function discretize(curve, numpanels, panelorder;
                    equal_arclength=true)
    if isa(curve, Array)
        # TODO: Would be nicer to use multiple dispatch
        curves = curve
        Ncurves = length(curves)
        @assert length(numpanels)==Ncurves
        dcurves = Array{DiscreteCurve}(Ncurves)
        for i=1:Ncurves
            dcurves[i] = discretize(curves[i], numpanels[i], panelorder, equal_arclength=equal_arclength)
        end
        return multiply_connected(dcurves)
    end    
    start = 0.0
    stop = 2*pi
    # Setup division in parameter
    if equal_arclength
        tsplit, _, _ = traparclen(curve, numpanels)
        push!(tsplit, 2*pi)
    else
        # Equal parameter length
        tsplit = linspace(start, stop, numpanels+1);
    end
    # Prepare structures
    numpoints = numpanels*panelorder
    glpoints, glweights = gausslegendre(panelorder)
    points = zeros(2, numpoints)
    normals = zeros(2, numpoints)
    derivatives = zeros(2, numpoints)
    derivatives2 = zeros(2, numpoints)
    weights = zeros(numpoints)
    dS = zeros(numpoints)
    curvature = zeros(numpoints)
    t_edges = zeros(2, numpanels)
    edges = zeros(4, numpanels)
    n_edges = zeros(4, numpanels)    
    splitcomplex(z) = [real(z), imag(z)]
    L = legendre_matrix(panelorder)
    resolution = 0.0
    # Populate
    for i in 1:numpanels
        # setup quadrature for t in [a,b]
        a, b = tsplit[i:i+1]
        t_edges[:, i] = [a, b]
        edges[1:2, i] = splitcomplex(curve.tau(a))
        edges[3:4, i] = splitcomplex(curve.tau(b))
        n_edges[1:2, i] = splitcomplex(curve.normal(a))
        n_edges[3:4, i] = splitcomplex(curve.normal(b))                        
        dt = b-a
        for j=1:panelorder
            idx=j + panelorder*(i-1)
            t = a + (1+glpoints[j])*dt/2
            w = glweights[j]*dt/2
            z = curve.tau(t)
            zp = curve.dtau(t)
            zpp = curve.d2tau(t)
            n = curve.normal(t)
            points[:, idx] = splitcomplex(z)
            derivatives[:, idx] = splitcomplex(zp)
            derivatives2[:, idx] = splitcomplex(zpp)
            normals[:, idx] = splitcomplex(n)
            weights[idx] = w
            dS[idx] = w*abs(zp)
            curvature[idx] = imag(conj(zp)*zpp/abs(zp)^3)
        end
        # Get resolution estimate
        idx = (1:panelorder) + panelorder*(i-1)
        coeff = L*(dS[idx]./weights[idx])
        thisresolution = abs(coeff[end]) + abs(coeff[end-1])
        resolution = max(resolution, thisresolution)
    end

    prevpanel = circshift(1:numpanels, 1)
    nextpanel = circshift(1:numpanels, -1)
    curvenum = ones(Int64, numpoints)
    
    info(@sprintf("  Grid resolution: %.2e\n", resolution))

    DiscreteCurve(
        panelorder,
        numpanels,
        numpoints,
        edges,
        n_edges,
        t_edges,
        points,
        normals,
        weights,
        dS,
        curvature,
        prevpanel,
        nextpanel,
        curvenum
    )
end

"""
Created a multiply connected curve discretization based on several 
simply connected ones.
"""
function multiply_connected(grids)
    # Init
    panelorder = grids[1].panelorder
    numpanels = 0
    numpoints = 0
    points = zeros(2, 0)
    normals = zeros(2, 0)
    derivatives = zeros(2, 0)
    derivatives2 = zeros(2, 0)
    weights = zeros(0)
    dS = zeros(0)
    curvature = zeros(0)
    edges = zeros(4, 0)
    n_edges = zeros(4, 0)    
    t_edges = zeros(2, 0)    
    prevpanel = Int64[]
    nextpanel = Int64[]
    curvenum = Int64[]
    for i=1:length(grids)
        assert(grids[i].panelorder == panelorder)
        # Concatenate
        edges = hcat(edges, grids[i].edges)
        n_edges = hcat(n_edges, grids[i].n_edges)        
        t_edges = hcat(t_edges, grids[i].t_edges)
        points = hcat(points, grids[i].points)
        normals = hcat(normals, grids[i].normals)
        weights = vcat(weights, grids[i].weights)
        dS = vcat(dS, grids[i].dS)
        curvature = vcat(curvature, grids[i].curvature)
        prevpanel = vcat(prevpanel, grids[i].prevpanel+numpanels)
        nextpanel = vcat(nextpanel, grids[i].nextpanel+numpanels)
        append!(curvenum, i*ones(Int64, grids[i].numpoints))
        # Update count
        numpanels += grids[i].numpanels
        numpoints += grids[i].numpoints        
    end
    # Construct output
    return DiscreteCurve(
        panelorder,
        numpanels,
        numpoints,
        edges,
        n_edges,
        t_edges,
        points,
        normals,
        weights,
        dS,
        curvature,
        prevpanel,
        nextpanel,
        curvenum
    )
end

function glpanels(a, b, numpanels, pedges, order)
    # Discretize interval using composite Gauss-Legendre quadrature
    # 
    # Original code by Rikard Ojala
    T, W = gausslegendre(order)
    t = zeros(order*numpanels)
    w = zeros(order*numpanels)
    ptr = 1
    for j = 1:numpanels
        t[ptr:ptr+order-1] = (pedges[j+1]-pedges[j])*T/2+(pedges[j]+pedges[j+1])/2;
        wj = W*(pedges[j+1]-pedges[j])/2
        w[ptr:ptr+order-1] = wj
        ptr = ptr + order
    end
    return t, w
end

function traparclen(curve, N, NL=1000)
    # t,z,W,L = traparclen(f, N, NL=1000)
    #
    # Discretizes the curve given by f(t) t = [0,2*pi] using the trapezoidal rule
    # equidistant in arc-length .
    # N is the number of quadrature points requested and if specified,
    # NL is the number of 16-point Gauss-Legendre panels used to compute the 
    # length of the curve. Default is NL = 1000.
    # Non-adaptive. Make sure that the curve is well resolved by N points.
    #
    # Original code by Rikard Ojala
    if N*2 > NL
        NL = 2*N
        info("Using $NL Gauss-Legendre points in traparclen")
    end
    order = 16
    pedges = linspace(0, 2*pi, NL+1)
    T, W = glpanels(0, 2*pi, NL, pedges, order)
    zp = curve.dtau.(T)
    L = W'*abs.(zp)    
    L2 = linspace(1/N,1-1/N,N-1)*L
    #Initial guess
    t = linspace(0,2*pi,N+1)
    t = t[2:end-1]
    dt = 1
    iter = 0
    maxiter = 30
    while norm(dt)/norm(t) > 1e-13 && iter < maxiter
        # Set up n-point GL quadrature on each segment.
        T, W = glpanels(0,2*pi,N-1,[0;t],order)
        zpt = curve.dtau.(t)
        zpT = curve.dtau.(T)        
        #Compute the cumulative sum of all the segment lengths
        F = cumsum(sum(reshape(W.*abs.(zpT),order,N-1), 1)')
        dt = (F-L2)./abs.(zpt)
        # Sort the parameters just in case. 
        t = vec(t - dt)
        sort!(t)        
        iter = iter + 1
    end
    if iter == maxiter
        warn("Newton for equal arclength did not converge, bad discretization")
    end
    t = [0;t]
    z = curve.tau.(t)
    W = 2*pi/N*ones(N,1)
    return t, z, W, L
end

function minmax_panel(grid::DiscreteCurve)
    hmin, hmax = Inf, 0.0
    for i=1:grid.numpanels
        idx = (i-1)*grid.panelorder + (1:grid.panelorder)
        h = sum(grid.dS[idx])
        hmin = min(hmin, h)
        hmax = max(hmax, h)
    end
    return hmin, hmax
end

function interior_points(grid::DiscreteCurve, zt::Array{Float64})
    # First pass: Stokes DLP identity works to O(1) very close to bdry
    ndS = grid.normals .* grid.dS'
    # # Hack: Approximate stokes DLP by modstokes DLP with small alpha    
    # density = ones(2, grid.numpoints)
    # u = fmm_stresslet_targ(grid.points, zt, density, ndS, 1e-8)
    # marker = u[1,:]
    # Better, use Laplace DLP
    density = ones(grid.numpoints)
    U = rfmm2d(source=grid.points, target=zt, dipstr=density, dipvec=ndS,
                     ifpot=false, ifpottarg=true)
    marker = -1/(2*pi)*U.pottarg    
    interior = marker .> 0.5
    # Second pass: find points close to boundary and distrust them
    _, hmax = minmax_panel(grid)
    dx = hmax*3/grid.panelorder # Estimate of largest point gap on bdry
    R = dx
    nnb, _ = nearest_nb_R(grid, zt, R)
    # listed points are within R of a boundary point
    # Third pass: check sign of projection against near normal
    veryclose = zeros(interior)
    N = size(zt, 2)
    for i=1:N
        if nnb[i] != 0
            rvec = zt[:,i] - grid.points[:,nnb[i]]
            n = grid.normals[:,nnb[i]]
            rdotn = rvec[1]*n[1] + rvec[2]*n[2]
            interior[i] = rdotn > 0
            veryclose[i] = abs(rdotn) < R/10
        end
    end
    # Fourth pass: Find projection using parametrization
    coeffs = map_panels(grid)
    for i=1:N
        if veryclose[i] != 0
            # Find corresponding panel
            near_panel = Int64(ceil(nnb[i]/grid.panelorder))
            z = zt[1,i] + 1im*zt[2,i]
            tloc, _ = invert_map(grid, coeffs, near_panel, z)
            interior[i] = imag(tloc) > 0
        end
    end
    return interior
end

function nearest_nb_R(grid::DiscreteCurve, zt::Array{Float64}, R::Float64)
    # Bin sort of boundary points
    delta = 1e-8
    xmax = max(maximum(zt[1,:]), maximum(grid.points[1,:])) + delta
    xmin = min(minimum(zt[1,:]), minimum(grid.points[1,:])) - delta
    ymax = max(maximum(zt[2,:]), maximum(grid.points[2,:])) + delta
    ymin = min(minimum(zt[2,:]), minimum(grid.points[2,:])) - delta
    Lx = xmax-xmin
    Ly = ymax-ymin
    Nx = Int64(floor(Lx/R))
    Ny = Int64(floor(Ly/R))
    hx = Lx/Nx
    hy = Ly/Ny
    Nx += 1
    Ny += 1
    bin_count, bin_lists = binsort_points(grid.points, xmin, ymin, hx, hy, Nx, Ny)
    # For each target points, find neareast neighbor
    R2 = R^2
    N = size(zt, 2)        
    nearest_nb = zeros(Int64, N)
    dist = zeros(N)    
    for trg_idx=1:N
        z = zt[:,trg_idx]        
        ih, jh = bin_idx(z, xmin, ymin, hx, hy) # home bin
        r2min = Inf
        imin = 0
        # Iterate over nb bins
        for inb = ih+(-1:1)
            for jnb = jh+(-1:1)
                if inb>0 && jnb>0 && inb<Nx+1 && jnb<Ny+1
                    # Compare to target points in bin
                    list = bin_lists[inb, jnb]
                    for k=1:length(list)
                        i = list[k]
                        r2 = (grid.points[1,i]-z[1])^2 + (grid.points[2,i]-z[2])^2
                        if r2<r2min && r2<R2
                            imin = i
                            r2min = r2
                        end
                    end
                end
            end
        end
        nearest_nb[trg_idx] = imin
        dist[trg_idx] = sqrt(r2min)
    end # nearest nb list built
    return nearest_nb, dist
end

function binsort_points(zt, xmin, ymin, hx, hy, Nx, Ny)    
    N = size(zt, 2)
    # 1. Count points in each bin
    bin_count = zeros(Int64, Nx, Ny)
    for i=1:N
        ix, iy = bin_idx(zt[:,i], xmin, ymin, hx, hy)
        bin_count[ix, iy] += 1
    end
    # 2. Setup bin lists
    fill_ptr = ones(Int64, Nx, Ny)
    bin_lists = Array{Array{Int64}}(Nx, Ny)
    for i=1:Nx
        for j=1:Ny
            bin_lists[i,j] = zeros(Int64, bin_count[i,j])
        end
    end
    # 3. Fill bins
    for i=1:N
        ix, iy = bin_idx(zt[:,i], xmin, ymin, hx, hy)
        ptr = fill_ptr[ix,iy]
        bin_lists[ix,iy][ptr] = i
        fill_ptr[ix,iy] += 1
    end
    return bin_count, bin_lists
end

@inline function bin_idx(z, xmin, ymin, hx, hy)
    ix = Int64(floor( (z[1]-xmin)/hx )) + 1
    iy = Int64(floor( (z[2]-ymin)/hy )) + 1
    return ix, iy
end

end
