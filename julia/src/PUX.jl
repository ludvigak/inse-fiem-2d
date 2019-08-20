# Based on Matlab source by Fredrik Fryklund and Erik Lehto

module PUX

mydir = Base.Filesystem.dirname(Base.source_path())
const PUX_MATLAB_SRC = "$(mydir)/../../external/pux"
const RBFQR_MATLAB_SRC = "$(mydir)/../../external/rbf-qr"

using MATLAB
using PyPlot
using CurveDiscretization

export pux_params, pux_setup!, pux_eval

function __init__()
    mat"clear()"
    mat"addpath(genpath($PUX_MATLAB_SRC))"    
    mat"addpath(genpath($RBFQR_MATLAB_SRC))"    
end

mutable struct PUXStorage
    numcenters::Int64
    PUXParams::Dict
    PUXData::Dict
    PUXDataInit::Bool
end

function pux_precompute(curves, dcurve, Ngrid, ep; P=0, R=0.0)
    Ngrid_int = Ngrid
    if P==0
        # Determine P from radius
        hgrid = 2*maximum(abs.(vec(dcurve.points))) / Ngrid
        P = ceil(R/hgrid)
    end
    # Setup volume grid
    Lgrid, Ngrid, X, Y = setup_volgrid(dcurve, Ngrid_int, P)
    xe = [vec(X) vec(Y)]
    xet = xe'
    interior = interior_points(dcurve, xet)        
    # Setup PUX params, one setup per boundary
    if !isa(curves, Array)
        curves = [curves]
    end
    Ncurves = length(curves)
    PUXStor = Array{PUXStorage}(Ncurves)
    for i=1:Ncurves
        curve_i = curves[i]
        idx = findin(dcurve.curvenum, i)
        arclength = sum(dcurve.dS[idx])    
        PUXStor[i] = pux_params(Lgrid/2, Ngrid, ep, P, arclength)
        # Get equispaced points
        _, z_equi = traparclen(curve_i, PUXStor[i].numcenters)
        z_grid = dcurve.points[1,idx] + 1im*dcurve.points[2,idx]
        # Precompute
        pux_setup!(xe, interior, z_grid, z_equi, PUXStor[i])
    end
    # Return everything that the solver needs
    return PUXStor, X, Y, Lgrid, Ngrid, interior
end

function pux_params(Lgrid, Ngrid, shape, support, arclength)
    M = Float64(Ngrid)
    LBox = Float64(Lgrid)    
    ep = Float64(shape)
    P = Float64(support)
    arcL = Float64(arclength)    
    # Percentage of parition radius when a zero-paritions will be removed
    R_ZeroPart_threshold = 0.8
    # radius of partition ratio (1 = coverage, >1 overlap)
    R_ratio = 2.5
    hgrid = 2.0*LBox/M

    R_part = P*hgrid
    R = R_part/R_ratio

    info("PUX PARAMS:")
    info("P = ", P)
    info("R = ", R)
    info("R_part = ", R_part)
    
    regularity = 0
    if regularity == 0
        regularity = min(floor(sqrt(P)-0.9), 5) # test (?)
    end
    # Number of Vogel points, i.e. RBF-centres. 
    n_RBFCentres = round(min(0.8*pi/4*P^2, 4*P))
    n_RBFCentres = min(n_RBFCentres, 200)
    info("Number of Vogel nodes: ", n_RBFCentres)

    # Number of partitions
    n_Part = ceil(arcL/(2*R)) + 1    
    
    # Pack into Dict
    PUXParams = Dict{String,Any}(
        "R_ZeroPart_threshold" => R_ZeroPart_threshold, 
        "R_part" =>   R_part, 
        "regularity" => regularity,  
        "n_RBFCentres" =>         n_RBFCentres,        
        "n_Part" =>                n_Part, 
        "M" => M, 
        "LBox" => LBox, 
        "ep" => ep, 
        "arcL" => arcL
    )
    # Return PUXStorgae
    stor = PUXStorage(PUXParams["n_Part"]*2, PUXParams, similar(PUXParams), false)
    return stor
end

function pux_setup!(xe::Array{Float64,2},
                    interior::BitArray{1},
                    z_grid::Vector{Complex{Float64}},
                    z_equi::Vector{Complex{Float64}},
                    stor::PUXStorage)
    params = stor.PUXParams
    mat"$PUXData = pux_setup($xe, $interior, $z_grid, $z_equi, $params)"
    stor.PUXData = PUXData
    stor.PUXDataInit = true
end

function pux_eval(f1, stor)
    # Dummy wrapper for vector PUX
    f2 = zeros(f1)
    u1, u2 = pux_eval_vec(f1, f2, stor)
    return u1
end

function pux_eval_vec(f1::Array{Float64},
                      f2::Array{Float64},
                      stor::PUXStorage);
                      limit_od::Float64=3.5
    pux_eval_vec(f1, f2, [stor])
end

function pux_eval_vec(f1::Array{Float64},
                      f2::Array{Float64},
                      stor_list::Array{PUXStorage};
                      limit_od::Float64=3.5)
    fe1 = zeros(size(f1))
    fe2 = zeros(size(f2))
    for curvenum = 1:length(stor_list)
        stor = stor_list[curvenum]
        @assert stor.PUXDataInit "Must run pux setup first!"
        data = stor.PUXData
        # Extraxt data from MATLAB struct
        # Everything is Float64 when we get it, convert explicitly
        A = data["A"]
        W = data["W"]
        If = map(Int64, data["If"])
        Jf = map(Int64, data["Jf"])
        n_Boxpnts = Int64(data["n_Boxpnts"])
        n_Part = Int64(data["n_Part"])
        n_InterpPart = Int64( data["n_InterpPart"] )
        n_extensionPnts = Int64(data["n_extensionPnts"])
        idx_xeOutOmegaInInterpPart = map(Int64, data["idx_xeOutOmegaInInterpPart"])
        idx_xeInOmega_boolean = map(Bool, data["idx_xeInOmega_boolean"])
        idx_InterpPartCentres = map(Int64, data["idx_InterpPartCentres"])
        idx_xeInInterpPart_stencil = map(Int64, data["idx_xeInInterpPart_stencil"])
        # Coarsen sets how to sample the evaluation grid to obtain the
        # interpolation grid. Set to 1 if they should be equal.
        coarsen = 1
        # Find values in collocation points
        f1_collocation = zeros(n_Boxpnts)
        f2_collocation = zeros(n_Boxpnts)    
        f1_collocation[idx_xeInOmega_boolean] = f1[idx_xeInOmega_boolean]
        f2_collocation[idx_xeInOmega_boolean] = f2[idx_xeInOmega_boolean]    
        # Create coarser grid for interpolation
        if (coarsen > 1)
            # NEEDS TO BE TESTED
            warn("Coarsening is untested in Julia implementation")
            M = Int64( stor.PUXParams["M"] )
            I, J = ndgrid(1:coarsen:M, 1:coarsen:M)
            I = I'
            J = J'
            i_crs = I+M*(J-1)
            i_crs = vec(i_crs)
        end
        # Vector for data used for building extension. The names S is
        # based on the input for MATLABS sparse:
        # S = sparse(i,j,s,m,n) uses vectors i, j, and s to generate an
        #     m-by-n sparse matrix such that S(i(k),j(k)) = s(k)
        # Sf stores information about the local extensions.
        Sf1 = zeros(n_extensionPnts)
        Sf2 = zeros(n_extensionPnts)    
        # Counters used to fill in Sf.
        idx_endf = 0
        max_resnorm = 0.0
        ext_rel_magn = 0.0
        for i = 1:n_InterpPart
            ## Global indices
            # Indices of point in xe within R_Ip of parition center Ip(i).
            idx_xeInInterpPart_i = idx_InterpPartCentres[i] + idx_xeInInterpPart_stencil
            idx_xeInOmegaInInterpPart_i = idx_xeInInterpPart_i[idx_xeInOmega_boolean[idx_xeInInterpPart_i]]
            idx_xeInInterpPartOutOmega = idx_xeInInterpPart_i[.~idx_xeInOmega_boolean[idx_xeInInterpPart_i]]
            n_local_evals = length(idx_xeInInterpPartOutOmega)
            # Local indices
            idx_local_xeInOmega = find(idx_xeInOmega_boolean[idx_xeInInterpPart_i])
            idx_local_xeOutOmega = find(.~idx_xeInOmega_boolean[idx_xeInInterpPart_i])
            
            if (coarsen > 1)
                # Coarser grid. Samples points set by coarsen.
                ic, ij = intersect(idx_xeInInterpPart_i[idx_xeInOmega_boolean[idx_xeInInterpPart_i]],i_crs) # global
                idx_xeInOmegaInInterpPart_i = ic[idx_xeInOmega_boolean[ic]]  # global
                idx_local_xeInOmega = idx_local_xeInOmega[ij]
            end
            # Solve overdetermined system to obtain f in RBF centres.
            MAT = A[idx_local_xeInOmega,:]
            rhs = hcat(f1_collocation[idx_xeInOmegaInInterpPart_i], f2_collocation[idx_xeInOmegaInInterpPart_i])

            # Limit overdeterminedness
            rows, cols = size(MAT)
            if rows/cols > limit_od
                newrows = limit_od*cols
                idx = 1:floor(Integer, rows/newrows):rows
                f_RBFCentres = MAT[idx, :]\rhs[idx, :]
            else
                f_RBFCentres = MAT\rhs
            end
            # Compute residual
            resnorm = rel = norm(MAT*f_RBFCentres - rhs) / norm(rhs)
            max_resnorm = max(max_resnorm, resnorm)
            # Extrapolate
            Af = A[idx_local_xeOutOmega,:] * f_RBFCentres
            Sf1[idx_endf+1:idx_endf+n_local_evals] = Af[:,1]
            Sf2[idx_endf+1:idx_endf+n_local_evals] = Af[:,2]        
            idx_endf = idx_endf + n_local_evals
            # Heuristic sanity check of extrapolation
            ext_rel_magn = max(ext_rel_magn, norm(Af, Inf) / norm(rhs, Inf))
        end
        #info("PUX: maximum residual: $max_resnorm")
        #info("PUX: extension rel. magnitude: $ext_rel_magn")
        if ext_rel_magn > 10
            #warn("PUX: Partitions too large?")
        end
        
        # cell_idx_xeOutOmegaInInterpPart is obtained by range search from 
        #partition centers on uniform grid. Thus they are not the same amount as 
        #idx_local_xeOutOmega_boolean. idx_local_xeOutOmega_boolean
        # comes from inpolygon. 
        Sf1 = Sf1[1:idx_endf]
        Sf2 = Sf2[1:idx_endf]    
        Imat1 = sparse(If, Jf, Sf1, n_Boxpnts, n_Part)
        Imat1 = Imat1[idx_xeOutOmegaInInterpPart,:]
        Imat2 = sparse(If, Jf, Sf2, n_Boxpnts, n_Part)
        Imat2 = Imat2[idx_xeOutOmegaInInterpPart,:]    
        ## Combine local extensions into global extension
        fe1[idx_xeOutOmegaInInterpPart] = full(sum(Imat1.*W, 2))
        fe1[idx_xeInOmega_boolean] = f1_collocation[idx_xeInOmega_boolean]
        fe2[idx_xeOutOmegaInInterpPart] = full(sum(Imat2.*W, 2))
        fe2[idx_xeInOmega_boolean] = f2_collocation[idx_xeInOmega_boolean]
    end
    return fe1, fe2
end

function plot_partitions(stor::PUXStorage)
    function partplot(clist, Rlist, style)
        th = linspace(0, 2*pi)
        for i=1:size(clist,1)
            c = clist[i,:]
            plot(c[1] + Rlist[i]*cos.(th), c[2] + Rlist[i]*sin.(th), style)
        end
    end
    PUXData = stor.PUXData
    partplot(PUXData["Ip"], PUXData["R_Ip"]*ones(PUXData["n_InterpPart"]), "r")
    partplot(PUXData["Op"], PUXData["R_Op"], "g-")
end    


function plot_partitions(arr::Array{PUXStorage})
    for s in arr
        plot_partitions(s)
    end
end

# Wrappers to original Matlab code
function pux_params_matlab(Lgrid, Ngrid, shape, support, arclength)    
    # Convert to Float64 (MATLAB default type)
    pux_ep = Float64(shape)
    pux_P = Float64(support)
    Lgrid = Float64(Lgrid)
    Ngrid = Float64(Ngrid)
    arclength = Float64(arclength)
    mat"$PUXParams = pux_params($Ngrid, $Lgrid, $pux_ep, $pux_P, $arclength)"
    stor = PUXStorage(PUXParams["n_Part"]*2, PUXParams, similar(PUXParams), false)
    return stor
end

function pux_eval_matlab(fin::Array{Float64},
                         stor::PUXStorage)
    data = stor.PUXData
    mat"$fext = pux_eval($fin, $data)"
    return reshape(fext, size(fin))
end

function setup_volgrid(grid::DiscreteCurve, Ngrid, pux_P)
    # Setup minimal grid size, ***assuming curve is centered at origin***
    Dmax = maximum(abs.(vec(grid.points)))
    Lgrid = 2*Dmax
    hgrid = Lgrid/Ngrid
    padding = 2*pux_P+2 # This should always accomodate PUs (which is right, 2*P or 4*P ???))
    Ngrid += padding
    Ngrid += Ngrid%2 # Even FFT grid is better
    info("After padding Ngrid=$Ngrid")
    Lgrid += padding*hgrid
    Lgrid *= exp(0.001)  # fudge factor to reduce number of grid points exactly on bdry
    xgrid = linspace(-Lgrid/2, Lgrid/2, Ngrid+1)
    ygrid = linspace(-Lgrid/2, Lgrid/2, Ngrid+1)
    xgrid = xgrid[1:end-1]
    ygrid = ygrid[1:end-1]
    X, Y = ndgrid(xgrid, ygrid)    
    Lgrid, Ngrid, X, Y
end

# Matlab compability 
function ndgrid(x, y)
    nx = length(x)
    ny = length(y)
    vartype = typeof(x[1])
    @assert typeof(y[1])==vartype
    X = Array{vartype}(nx, ny)
    Y = Array{vartype}(nx, ny)
    for i=1:nx
        for j=1:ny
            X[i,j] = x[i]
            Y[i,j] = y[j]
        end
    end
    return X,Y
end

end
