# Add const here when dev is over
NEARLIM_H_16 = 1.2
UPSAMPLIM_H_16 = 0.4
 
NEARLIM_RADIUS_16 = 3.5
UPSAMP_RADIUS_16 = 3.0

EXRAD_H = 0.05 # FMM exclusion radius / panel length
EDGELIM_H = 0.05 # Panel edge near criterion

# Bernstein radius computed in first quadrant
bernstein_radius_q1(z) = abs(z + sqrt(z^2 - 1))
q1(z) = abs(real(z)) + 1im*abs(imag(z))
bernstein_radius(z) = bernstein_radius_q1(q1(z))

rotate_and_scale(za, zb, z) = (za+zb)/(za-zb) + z*2/(zb-za)
complex2real(z) = [real(z), imag(z)]    
tovec(x) = vec(x)
fromvec(x) = reshape(x, 2, :)

# Evaluate DLP with correction
function doublelayer(grid::DiscreteCurve,
                     density::Array{Float64},
                     target::Array{Float64},
                     alpha::Float64;
                     specquad=true)
    u = [0.0, 0.0]
    n = grid.panelorder
    coeffs = CurveDiscretization.map_panels(grid)
    for panel_idx = 1:grid.numpanels
        (up, rmin2, h) = dlp_panel_direct(grid, density, target,
                                          alpha, panel_idx)
        rmin = sqrt(rmin2)
        if specquad && rmin < 2*NEARLIM_H_16 *h*16/n
            # Nearby, investigate further:
            # Compute target point Bernstein radius
            z = target[1] + 1im*target[2]
            troot, _, converged = CurveDiscretization.invert_map(grid, coeffs, panel_idx, z,
                                                                 fall_back_to_initial=false)
            
            rho = bernstein_radius(troot)
            # Empirical thresholds for n=16
            
            #if specquad && rmin < NEARLIM_H_16 *h*16/n
            #upsample = rmin < UPSAMPLIM_H_16 *h*16/n        
            if rho < NEARLIM_RADIUS_16^(16/n) # TODO: implement actual estimates here
                #@show troot, converged
                # Compute correction for this panel
                ucorr = dlp_panel_specquad(grid, density, target, alpha,
                                           panel_idx, troot)
                up = ucorr
            end # End near eval
        end
        u += up
    end # End panel loop
    return u
end

@inline function dlp_panel_direct(grid::DiscreteCurve,
                                  density::Array{Float64},
                                  target::Array{Float64},
                                  alpha::Float64,
                                  panel_idx::Int64)
    up = [0.0, 0.0]    
    rmin2 = Inf
    h = 0.0
    for idx=grid.panelorder*(panel_idx-1) + (1:grid.panelorder)
        rvec1 = grid.points[1,idx]-target[1]
        rvec2 = grid.points[2,idx]-target[2]
        rnorm2 = rvec1*rvec1 + rvec2*rvec2
        rmin2 = min(rmin2, rnorm2)
        h += grid.dS[idx]
        ui = stresslet([rvec1, rvec2],
                       density[:,idx],
                       grid.normals[:,idx],
                       alpha)
        up += ui*grid.dS[idx]
    end
    return (up, rmin2, h)
end

function dlp_panel_specquad_weights(grid,za, zb, zt, zj,
                                    nuj, kappa,
                                    zpwj, zpj,
                                    troot,
                                    alpha::Float64,
                                    panel_idx::Int64,
                                    ifgrad::Bool,                                    
                                    omega,
                                    DD,
                                    omega1, omega2,
                                    glpoints1, glweights1,
                                    glpoints2, glweights2,
                                    bcweights
                                    )
    n = length(glpoints1)
    n_sub = length(glpoints2)        
    panel_length = sum(abs.(zpwj))
    dt = real(sum(zpwj ./ zpj))
    # Subdivide panel
    dt_max = LARGEALPHA_LIMIT_32*2/(alpha*panel_length)    
    intervals = subdivide_interval_with_bisection(dt_max, troot, n_sub)
    # Compute weights on intervals
    S = zeros(2, 2*n)
    if ifgrad
        DUDXspec, DUDYspec = complex(zeros(1, 2*n)), complex(zeros(1, 2*n))
    end
    for idx = 1:length(intervals)-1
        ta = intervals[idx]
        tb = intervals[idx+1]
        troot_sub = (troot-ta)*2/(tb-ta) - 1.0
        # New source points, in original frame
        tsrc = ta + (1+glpoints2)*(tb-ta)/2
        # Interpolation to new source points
        Psub = bclag_interp_matrix(glpoints1, tsrc, bcweights)
        dt_sub = tb-ta
        if ta==-1 && tb==1
            za_sub, zb_sub = za, zb
        else
            za_sub = bclag_interp_eval(glpoints1, zj, ta, bcweights)
            zb_sub = bclag_interp_eval(glpoints1, zj, tb, bcweights)
        end
        nuj_sub = Psub*nuj
        zj_sub = Psub*zj
        zpj_sub = Psub*zpj
        kappa_sub = Psub*kappa
        zpwj_sub = zpj_sub .* glweights2 * dt/2 * (dt_sub/2)
        # nua and nub are not used, do this to make sure we notice if they are
        nua = nothing
        nub = nothing        
        # Everything in place, do computation
        near_sub = bernstein_radius(troot_sub) < NEARLIM_RADIUS_16^(16/n_sub)
        if near_sub
            S_sub, DUDXspec_sub, DUDYspec_sub =
                complex_specquad_weights(zj_sub, zpwj_sub, nuj_sub, kappa_sub, za_sub, zb_sub, nua, nub,
                                         omega, omega1, omega2, zt, alpha,
                                         DD, n_sub, ifgrad)
        else
            S_sub = complex_direct_weights(zj_sub, zpwj_sub, nuj_sub, zt, alpha)
            if ifgrad
                DUDXspec_sub, DUDYspec_sub =
                    gradient_kernel_direct_weights(zj_sub, zpwj_sub, nuj_sub, zt, alpha)
            end
        end
        S += S_sub*kron(Psub, eye(2))
        if ifgrad
            DUDXspec += DUDXspec_sub*kron(Psub, eye(2))
            DUDYspec += DUDYspec_sub*kron(Psub, eye(2))
        end
    end        
    # Join weights into one matrix and return
    if ifgrad
        W = vcat(S,
                 real(DUDXspec),
                 real(DUDYspec),
                 imag(DUDXspec),
                 imag(DUDYspec))
        return W
    else
        return S
    end
end

function complex_direct_weights(zj, zpwj, nuj, zt, alpha)
    # Build weights from panel to point
    n = length(zj)
    S = zeros(2, 2*n)
    for i=1:n
        r1 = real(zj[i]-zt)
        r2 = imag(zj[i]-zt)        
        n1 = real(nuj[i])
        n2 = imag(nuj[i])        
        dS = abs(zpwj[i])
        K11, K12, K21, K22 = dblkernel_noarrays(r1, r2, n1*dS, n2*dS, alpha)
        S[1, 1 + 2*(i-1)] = K11
        S[1, 2 + 2*(i-1)] = K12
        S[2, 1 + 2*(i-1)] = K21
        S[2, 2 + 2*(i-1)] = K22
    end
    return S
end

# Compute weights, using specquad if needed
function complex_specquad_weights(zj, zpwj, nuj, kappa, za, zb, nua, nub,
                                  omega, omega1, omega2, zt, alpha,
                                  DD, n, ifgrad)
    # Compute near eval coefficients
    wL, wC, wQ, wT = modifiedweights(za, zb, zt, zj)
    (wcorrL, wcorrC, wQ1, wQ2) =
        near_eval_weights(za, zb, zt, zj, nuj, zpwj, wL, wC, wQ)
    # Build weights from panel to point
    ucorr = [0.0, 0.0]    
    S = zeros(2, 2*n)
    Kcorr = zeros(2,2)
    for i=1:n
        r1 = real(zj[i]-zt)
        r2 = imag(zj[i]-zt)        
        n1 = real(nuj[i])
        n2 = imag(nuj[i])
        rnorm = sqrt(r1*r1 + r2*r2)
        rdotn = r1*n1 + r2*n2
        dS = abs(zpwj[i])
        # Direct quadrature:
        K11, K12, K21, K22 = dblkernel_noarrays(r1, r2, n1*dS, n2*dS, alpha)
        Kcorr[1,1] = K11
        Kcorr[1,2] = K12
        Kcorr[2,1] = K21
        Kcorr[2,2] = K22
        # Special quadrature correction:        
        # Components of split stresslet
        uS11, uS21, uL11, uL21, uC11, uC21, uQ = stresslet_split_noarrays(r1, r2, 1.0, 0.0, n1, n2, alpha)
        uS12, uS22, uL12, uL22, uC12, uC22, uQ = stresslet_split_noarrays(r1, r2, 0.0, 1.0, n1, n2, alpha)
        
        # # Compute value as direct sum + corrections
        # # This could in theory avoid some cancellation,
        # # since uS is not explicitly used, and uL only
        # # multiplies a quadrature correction.

        # Actually, this avoids cancellation for targets close to nodes
        Kcorr[1,1] = uS11*dS
        Kcorr[1,2] = uS12*dS
        Kcorr[2,1] = uS21*dS
        Kcorr[2,2] = uS22*dS

        # Log
        wLi = wcorrL[i] #- log(rnorm)*dS
        Kcorr[1,1] += wLi*uL11
        Kcorr[2,1] += wLi*uL21
        Kcorr[1,2] += wLi*uL12
        Kcorr[2,2] += wLi*uL22
        # Cauchy
        wCi = wcorrC[i] #- rdotn/rnorm^2*dS
        Kcorr[1,1] += wCi*uC11
        Kcorr[2,1] += wCi*uC21
        Kcorr[1,2] += wCi*uC12        
        Kcorr[2,2] += wCi*uC22
        # "Stresslet"
        A = wQ1[i]+wQ2[i]
        B = wQ1[i]-wQ2[i]
        Kcorr[1,1] += uQ*( real(A) )#- rdotn/rnorm^4*dS*r1*r1)
        Kcorr[1,2] += uQ*(-imag(B) )#- rdotn/rnorm^4*dS*r1*r2)
        Kcorr[2,1] += uQ*( imag(A) )#- rdotn/rnorm^4*dS*r1*r2)
        Kcorr[2,2] += uQ*( real(B) )#- rdotn/rnorm^4*dS*r2*r2)
        block_2x2_copy!(S, 1, i, Kcorr)            
    end
    ## Gradient
    if ifgrad
        # Gradient spec weights
        partial_integration = true            
        DUDX, DUDY =
            gradient_kernel_weights(zj, zpwj, nuj, kappa, za, zb, nua, nub,
                                    omega, omega1, omega2, zt, alpha,
                                    wcorrL, wC, wQ, wT, DD, partial_integration)
        return S, DUDX, DUDY
    else
        return S, [], []
    end        
end    

# On the fly computation of special quadrature
function dlp_panel_specquad(grid::DiscreteCurve,
                            density::Array{Float64},
                            target::Array{Float64},
                            alpha::Float64,
                            panel_idx::Int64,
                            troot)
    # For the correction weights we need complex arrays of
    # points on this panel
    n = grid.panelorder    
    idx = n*(panel_idx-1) + (1:n)            
    edges = grid.edges[:, panel_idx]
    za = edges[1] + 1im*edges[2]
    zb = edges[3] + 1im*edges[4]
    zt = target[1] + 1im*target[2]
    zj = grid.points[1,idx] + 1im*grid.points[2,idx]
    nuj = grid.normals[1,idx] + 1im*grid.normals[2,idx]
    weights = grid.weights[idx]    
    zpwj = grid.dS[idx] .*nuj./1im
    omega = density[1,idx] + 1im*density[2,idx]
    # Resampling stuff
    Nsec = 32 # Fixed order  
    glpoints1, glweights1 = gausslegendre(n)
    glpoints2, glweights2 = gausslegendre(Nsec)
    bcweights = bclag_interp_weights(glpoints1)    
    # Panel subdivision
    # Subdivide panel if too long to resolve split kernel,
    # and target point in near quad region
    panel_length = sum(abs.(zpwj))
    dt_max = LARGEALPHA_LIMIT_32*2/(alpha*panel_length)
    intervals = subdivide_interval_with_bisection(dt_max, troot, Nsec)
    #intervals = subdivide_interval(dt_max, troot)
    # Go through subpanels and add contribution
    u = [0.0, 0.0]
    for idx = 1:length(intervals)-1
        ta = intervals[idx]
        tb = intervals[idx+1]
        troot_sec = (troot-ta)*2/(tb-ta) - 1.0
        # New source points, in original frame
        tsrc = ta + (1+glpoints2)*(tb-ta)/2
        # Interpolation to new source points
        Psec = bclag_interp_matrix(glpoints1, tsrc, bcweights)
        zj_sec = Psec*zj
        nuj_sec = Psec*nuj
        omega_sec = Psec*omega
        za_sec = bclag_interp_eval(glpoints1, zj, ta, bcweights)
        zb_sec = bclag_interp_eval(glpoints1, zj, tb, bcweights)
        zpwj_sec = glweights2*(tb-ta)/2 .* (Psec*(zpwj./glweights1))
        u += dlp_specquad_complex(za_sec, zb_sec, zt, zj_sec, nuj_sec, zpwj_sec, omega_sec, alpha, troot_sec)
    end   
    return u
end

function dlp_specquad_complex(za, zb, zt, zj, nuj, zpwj, omega, alpha, troot)
    N = length(zj)
    rho = bernstein_radius(troot)
    near = rho < NEARLIM_RADIUS_16^(16/N)
    if near
        # Compute near eval coefficients
        (wcorrL, wcorrC, wQ1, wQ2) =
            near_eval_weights(za, zb, zt, zj, nuj, zpwj)
    end
    ucorr = [0.0, 0.0]
    for i=1:N
        r = complex2real(zj[i]-zt)
        f = complex2real(omega[i])
        n = complex2real(nuj[i])
        dS = abs(zpwj[i])
        u = stresslet(r, f, n, alpha)
        ucorr += u*dS                
        (uS, uL, uC, uQ) = stresslet_split(r, f, n, alpha)
        rnorm = norm(r)
        rdotf = dot(r,f)
        rdotn = dot(r,n)
        if near
            # Compute value as direct sum + corrections
            # This could in theory avoid som cancellation,
            # since uS is not explicitly used, and uL only
            # multiplies a quadrature correction.
            # Log
            ucorr += (wcorrL[i] - log(rnorm)*dS)*uL
            # Cauchy
            ucorr += (wcorrC[i] - rdotn/rnorm^2*dS)*uC
            # "Stresslet"        
            tmp = omega[i]*wQ1[i] + conj(omega[i])*wQ2[i]
            ucorr += (complex2real(tmp) - rdotf*r*rdotn/rnorm^4*dS)*uQ
        end
    end
    return ucorr
end


# FMM evaluation of DLP, with near corrections
function doublelayer_fast(grid::DiscreteCurve,
                          density::Array{Float64},
                          targets::Array{Float64},
                          alpha::Float64;
                          ifgrad::Bool=false,
                          specquad::Bool=true,
                          weights=[],
                          precomp=())
    # First compute double layer at all points using FMM
    ndS = grid.normals .* grid.dS'    
    if isempty(precomp)
        fmmpars, tree, sorted_pts = doublelayer_precomp(grid, targets, alpha)
    else
        fmmpars, tree, sorted_pts = precomp
    end
    FMM = fmm_stresslet_targ(fmmpars, tree, sorted_pts, density, ndS, alpha, ifgrad=ifgrad)
    if ifgrad
        u, u1grad, u2grad = FMM
    else
        u, u1grad, u2grad  = FMM, Float64[], Float64[]
    end
    
    if specquad
        if length(weights)==0
            weights = doublelayer_near_weights(grid, targets, alpha, ifgrad=ifgrad)
        end
        doublelayer_near_apply!(grid, density, u, u1grad, u2grad, weights, alpha, ifgrad)
    end
    if ifgrad
        return u, u1grad, u2grad
    else
        return u
    end
end

function doublelayer_precomp(grid, targets, alpha)
    exrad = Array{Float64}(undef, grid.numpoints)
    panel_lengths = Array{Float64}(undef, grid.numpanels)
    for i=1:grid.numpanels
        h = 0.0
        for j=1:grid.panelorder
            h += grid.dS[(i-1)*grid.panelorder + j]
        end
        for j=1:grid.panelorder
            exrad[(i-1)*grid.panelorder + j] = EXRAD_H*h
        end
    end
    return fmm_stresslet_prep(grid.points, targets, alpha, exrad=exrad)
end

function doublelayer_near_weights(grid::DiscreteCurve,
                                  targets::Array{Float64},
                                  alpha::Float64;
                                  ifgrad::Bool=false)
    # For each point on curve, keep a list of tuples: (target_idx, weights)
    WeightElement = Tuple{Int64,Array{Float64,2}}
    weights = Array{Array{WeightElement}}(grid.numpanels)
    for i=1:grid.numpanels
        weights[i] = Tuple[]
    end
    # Find which points are "near"
    # Limits experimentally determined for n=16
    nearlim_h = NEARLIM_H_16^(16/grid.panelorder)
    _, hmax = minmax_panel(grid)
    nnb, dist = nearest_nb_R(grid, targets, nearlim_h*hmax)
    numeval = size(targets, 2)
    # Compute Legendre coeffs of panels
    coeffs = CurveDiscretization.map_panels(grid)    
    # Precompute upsampling quantities
    n = grid.panelorder
    n2 = 2n
    glpoints1, glweights1 = gausslegendre(n)
    glpoints2, glweights2 = gausslegendre(n2)
    # Barycentric interpolation weights
    bcweights1 = bclag_interp_weights(glpoints1)
    bcweights2 = bclag_interp_weights(glpoints2)
    # Interpolation matrices for endpoints values
    B1n2 = vec(bclag_interp_matrix(glpoints2, -1)).'
    B2n2 = vec(bclag_interp_matrix(glpoints2,  1)).'
    # Real to complex matrix
    OMEGAn2 = kron(eye(n2), [1.0 1.0im])            
    # On-panel differentiation matrix
    DDn2 = legendre_diff_matrix(n2)
    # Matrices to send
    OMEGA1n2 = B1n2*OMEGAn2
    OMEGA2n2 = B2n2*OMEGAn2
    # Main loop over points
    for i=1:numeval
        if nnb[i] != 0
            near_point_idx = nnb[i];
            near_point = grid.points[:, near_point_idx]
            rmin = dist[i]
            target = targets[:,i]
            zt = target[1] + 1im*target[2]           
            near_panel_idx = Int64(ceil(near_point_idx/grid.panelorder))
            # Check three closest panels
            neighbors = (near_panel_idx,
                         grid.prevpanel[near_panel_idx],
                         grid.nextpanel[near_panel_idx])
            for nb_idx = 1:3
                # Load nodes on panel
                panel_idx = neighbors[nb_idx]
                idx = n*(panel_idx-1) + (1:n)            
                zj = grid.points[1,idx] + 1im*grid.points[2,idx]
                panel_length = sum(grid.dS[idx])
                # Check root (preimage) of target for this panel
                troot, _, converged = CurveDiscretization.invert_map(grid, coeffs, panel_idx, zt,
                                                                     fall_back_to_initial=false)                    
                rho = bernstein_radius(troot)
                rmin2 = minimum(abs2.(zt - zj))
                # Do nothing if ouside of near quadrature limit AND not excluded by FMM
                if rho > NEARLIM_RADIUS_16^(16/grid.panelorder) && rmin >= (EXRAD_H*panel_length)^2
                    continue
                end
                # Load panel data
                za = grid.edges[1, panel_idx] + 1im*grid.edges[2, panel_idx]
                zb = grid.edges[3, panel_idx] + 1im*grid.edges[4, panel_idx]
                idx = n*(panel_idx-1) + (1:n)                            
                nuj = grid.normals[1,idx] + 1im*grid.normals[2,idx]
                kappa = grid.curvature[idx]    
                zpwj = grid.dS[idx] .* nuj./1im
                zpj = zpwj ./ grid.weights[idx]
                # Get direct weights, for subtraction
                if ifgrad
                    DUDXdir, DUDYdir = gradient_kernel_direct_weights(zj, zpwj, nuj, zt, alpha)
                end
                D = complex_direct_weights(zj, zpwj, nuj, zt, alpha)               
                # Corner case: Points too close to discretization node can cause cancellation errors
                # FMM avoids points within EXRAD_H*panel_length, so no subtraction needed
                for j=1:n
                    if abs2(zj[j]-zt) < ( EXRAD_H*panel_length )^2
                        D[1:2, (1:2) + 2*(j-1)] = 0.0
                        if ifgrad
                            DUDXdir[1, (1:2) + 2*(j-1)] = 0.0
                            DUDYdir[1, (1:2) + 2*(j-1)] = 0.0
                        end
                    end
                end
                if ifgrad
                    Wsub = vcat(-D,
                                real(-DUDXdir),
                                real(-DUDYdir),
                                imag(-DUDXdir),
                                imag(-DUDYdir))
                else
                    Wsub = -D
                end
                # Check if too close to an edge
                prev_panel_length = sum(grid.dS[n*(grid.prevpanel[panel_idx]-1) + (1:n)])
                next_panel_length = sum(grid.dS[n*(grid.nextpanel[panel_idx]-1) + (1:n)])
                left_edge_rad  = EDGELIM_H * min(panel_length, prev_panel_length)
                right_edge_rad = EDGELIM_H * min(panel_length, next_panel_length)
                if abs(zt-za) > left_edge_rad && abs(zt-zb) > right_edge_rad
                    # Not close to an edge, compute special quadrature as usual
                    Wspec = dlp_panel_specquad_weights(grid,za, zb, zt, zj,
                                                       nuj, kappa,
                                                       zpwj, zpj,
                                                       troot, alpha, panel_idx, ifgrad,
                                                       OMEGAn2, DDn2,
                                                       OMEGA1n2, OMEGA2n2,
                                                       glpoints1, glweights1,
                                                       glpoints2, glweights2,
                                                       bcweights1
                                                       )
                    push!(weights[panel_idx], (i, Wspec+Wsub))
                else
                    if abs(zt-za) <= left_edge_rad
                        # Close to left edge, just subtract direct interaction
                        # Neighbor will deal with special quadrature
                        push!(weights[panel_idx], (i, Wsub))
                        continue
                    end
                    # Close to right edge, merge with right panel
                    next_panel_idx = grid.nextpanel[panel_idx]
                    zb = grid.edges[3, next_panel_idx] + 1im*grid.edges[4, next_panel_idx]                    
                    # index list spanning both panel
                    idx_wide = vcat(n*(panel_idx-1) + (1:n), n*(next_panel_idx-1) + (1:n))
                    # panels may be of different length in parameter
                    dt1 = grid.t_edges[2, panel_idx] - grid.t_edges[1, panel_idx]
                    dt2 = grid.t_edges[2, next_panel_idx] - grid.t_edges[1, next_panel_idx]
                    dt = dt1+dt2
                    ta, tb, tc = -1, dt1*2/dt-1, 1
                    g1 = ta + (1+glpoints1)*(tb-ta)/2
                    g2 = tb + (1+glpoints1)*(tc-tb)/2
                    # GL nodes of both panels, scaled into [-1,1]
                    glpoints_wide = vcat(g1, g2)
                    # Now we'll interpolate to 2*n nodes,
                    # find which fall into the first panel
                    pts1 = count(glpoints2 .< tb)
                    # Create interpolation from the two panels to nodes on new
                    P1 = bclag_interp_matrix(g1, glpoints2[1:pts1])
                    P2 = bclag_interp_matrix(g2, glpoints2[pts1+1:n2])
                    Pwide = zeros(n2,2n)
                    Pwide[1:pts1, 1:n] = P1
                    Pwide[pts1+1:n2, n+(1:n)] = P2
                    # Load boundary data for both panels
                    zj_wide = grid.points[1,idx_wide] + 1im*grid.points[2,idx_wide]
                    nuj_wide = grid.normals[1,idx_wide] + 1im*grid.normals[2,idx_wide]
                    kappa_wide = grid.curvature[idx_wide]    
                    zpwj_wide = grid.dS[idx_wide] .* nuj_wide./1im
                    zpj_wide = zpwj_wide ./ grid.weights[idx_wide]
                    # Interpolate to merged panel
                    zj = Pwide*zj_wide
                    nuj = Pwide*nuj_wide
                    kappa = Pwide*kappa_wide
                    zpj = Pwide*zpj_wide
                    zpwj = zpj*dt/2.*glweights2
                    # Find root of zt for merged panel (safer than scaling)
                    coeff = legendre_matrix(n2)*rotate_and_scale(za, zb, zj)
                    troot_guess = 1im*imag(troot)/2
                    troot_new, _, _ = newton_legendre(coeff, rotate_and_scale(za, zb, zt),
                                                      troot_guess, 30, 1e-13)
                    # Get specquad weights for merged panel
                    Wwide = dlp_panel_specquad_weights(grid,za, zb, zt, zj,
                                                       nuj, kappa,
                                                       zpwj, zpj,
                                                       troot_new, alpha, panel_idx, ifgrad,
                                                       OMEGAn2, DDn2,
                                                       OMEGA1n2, OMEGA2n2,
                                                       glpoints2, glweights2,
                                                       glpoints2, glweights2,
                                                       bcweights2
                                                       )
                    # Multiply weights with interpolation and split to the right source panels
                    # Also subtract direct interaction from current panel
                    Wwide = Wwide * kron(Pwide, eye(2))
                    Wwide1 = Wwide[:,1:2*n]
                    Wwide2 = Wwide[:,2*n+1:4*n]
                    push!(weights[panel_idx], (i, Wwide1+Wsub))
                    push!(weights[next_panel_idx], (i, Wwide2))                    
                end
            end
        end
    end
    return weights
end


# Apply precomputed near eval weights
function doublelayer_near_apply!(grid::DiscreteCurve,
                                 density::Array{Float64},
                                 u, u1grad, u2grad,
                                 weights,
                                 alpha::Float64,
                                 ifgrad)
    n = grid.panelorder
    for panel_idx = 1:grid.numpanels
        for i=1:length(weights[panel_idx])
            point_idx = weights[panel_idx][i][1]
            Wspec = weights[panel_idx][i][2]
            uspec1, uspec2 = 0.0, 0.0
            idx0 = 2*n*(panel_idx-1)
            if ifgrad
                uspec1, uspec2, du1dx, du1dy, du2dx, du2dy = matvec_6xN(Wspec, density, idx0, 2*n)
                u1grad[1,point_idx] += du1dx
                u1grad[2,point_idx] += du1dy
                u2grad[1,point_idx] += du2dx
                u2grad[2,point_idx] += du2dy                
            else
                uspec1, uspec2 = matvec_2xN(Wspec, density, idx0, 2*n)
            end
            u[1,point_idx] += uspec1
            u[2,point_idx] += uspec2                            
        end
    end
end
# Helper functions for the above
function matvec_2xN(A, x, idx, n)
    f1,f2 = 0.0, 0.0
    @simd for i=1:n
        @inbounds f1 += A[1,i]*x[idx+i]
        @inbounds f2 += A[2,i]*x[idx+i]
    end
    return f1,f2
end
function matvec_6xN(A, x, idx, n)
    f1,f2,f3,f4,f5,f6 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    @simd for i=1:n
        @inbounds f1 += A[1,i]*x[idx+i]
        @inbounds f2 += A[2,i]*x[idx+i]
        @inbounds f3 += A[3,i]*x[idx+i]
        @inbounds f4 += A[4,i]*x[idx+i]
        @inbounds f5 += A[5,i]*x[idx+i]
        @inbounds f6 += A[6,i]*x[idx+i]        
    end
    return f1,f2,f3,f4,f5,f6
end


# DLP on curve
# !This does not include correction weights for log kernel
function doublelayer_self(grid::DiscreteCurve,
                          density::Array{Float64},
                          alpha::Float64;
                          zerolimit::Bool=true)
    npoints = grid.numpoints
    u = zeros(2,npoints)
    for i=1:npoints
        ui = 0
        zi = grid.points[:,i]
        for j = 1:i-1
            rij = grid.points[:, j] - zi
            ui += stresslet(rij, density[:,j],
                            grid.normals[:,j], alpha) * grid.dS[j]
        end
        for j = i+1:npoints
            rij = grid.points[:, j] - zi
            ui += stresslet(rij, density[:,j],
                            grid.normals[:,j], alpha) * grid.dS[j]
        end
        u[:,i] = ui
    end
    if zerolimit
        for i=1:npoints        
            # Self
            ni = grid.normals[:,i]
            ti = [ni[2],-ni[1]]
            T2C0 = -1/pi
            u[:,i] += -T2C0*grid.curvature[i]*0.5*grid.dS[i] *
                dot(ti, density[:,i])*ti
        end
    end
    return u
end

"Return dense system matrix"
function system_matrix(grid::DiscreteCurve,
                       alpha::Float64)
    nvecdS = vec(grid.dS' .* grid.normals)
    nvec = vec(grid.normals)
    L = sum(grid.dS)
    Dmat = doublelayer_matrix(grid, alpha)
    ## Add add rank correction
    # NNT = nvec*nvecdS' / L
    # Dmat += NNT
    @. Dmat += nvec*nvecdS' / L
    # Add identity
    #Dmat += eye(2*grid.numpoints)/2
    for i=1:2*grid.numpoints
        Dmat[i,i] += 1/2
    end
    return  Dmat
end

"Return system matvec function handle"
function system_matvec(grid::DiscreteCurve,
                       alpha::Float64)
    ndS = grid.normals .* grid.dS'
    L = sum(grid.dS)    
    fmmpars, tree, sorted_pts = fmm_stresslet_prep(grid.points, grid.points, alpha)
    corrmat = correction_block_matrix(grid, alpha)
    return x -> x/2 + corrmat*x +
        tovec( fmm_stresslet_self(fmmpars, tree, sorted_pts, fromvec(x), ndS, alpha) ) +
        tovec(grid.normals)*dot(x, tovec(ndS))/L  
end

"Return block matrix that corrects a naive direct sum"
function correction_block_matrix(grid::DiscreteCurve,
                                 alpha::Float64)
    return specquad_block_matrix(grid, alpha; correction=true)
end

"Return block matrix containing elements of full system matrix for
neighboring panels"
function system_block_matrix(grid::DiscreteCurve,
                                 alpha::Float64)
    return specquad_block_matrix(grid, alpha;
                                 correction=false,
                                 system_matrix = true)
end

function specquad_block_matrix(grid::DiscreteCurve,
                               alpha::Float64;
                               correction::Bool=false,
                               system_matrix::Bool=false)
    # Compute sparse block matrix that contains all interactions requiring
    # special quadrature
    #
    # correction=true:
    # Subtract naive point-to-point sum (excluding self interactions) from
    # blocks, such that results can be used to correct an FMM
    #
    # system_matrix=true:
    # Include dbl layer jump and null space correction
    
    npoints = grid.numpoints
    panelorder = grid.panelorder
    L = sum(grid.dS)
    # Allocate arrays for constructing sparse matrix
    # only make room for nearest neighbor interactions,
    # this will have to grow for close geometries
    blocksize = 2*panelorder
    nzblocks = 3*grid.numpanels 
    indptr = Array{Int64}(grid.numpanels+1)
    indices = Array{Int64}(nzblocks)
    data = zeros(nzblocks, blocksize, blocksize)
    B1 = Array{Float64}(blocksize, blocksize) # Panel-to-panel block
    B2 = Array{Float64}(blocksize, blocksize) # Panel-to-panel block    
    idx = 1;
    for trg_panel = 1:grid.numpanels
        indptr[trg_panel] = idx -1 # -1 for Python indexing
        neighbors = (grid.prevpanel[trg_panel],
                     trg_panel,
                     grid.nextpanel[trg_panel])
        for src_panel = neighbors
            if correction
                # Compute block that represents naive sum
                p2p_block!(B1, grid, src_panel, trg_panel, alpha, specquad=false, upsampling=false)
                if src_panel==trg_panel
                    for i=1:panelorder
                        B1[2*(i-1)+(1:2), 2*(i-1)+(1:2)] = 0.0
                    end
                end
            end
            # Compute correct block
            p2p_block!(B2, grid, src_panel, trg_panel, alpha)
            if system_matrix
                # Add null space correction
                for i=1:panelorder
                    trg_idx = i + panelorder*(trg_panel-1)
                    ni1, ni2 = grid.normals[1:2, trg_idx]
                    for j=1:panelorder
                        src_idx = j + panelorder*(src_panel-1)
                        nj1, nj2 = grid.normals[1:2, src_idx]
                        dSoverL = grid.dS[src_idx] / L
                        B2[2*(i-1)+1, 2*(j-1)+1] += ni1*nj1*dSoverL
                        B2[2*(i-1)+2, 2*(j-1)+1] += ni2*nj1*dSoverL
                        B2[2*(i-1)+1, 2*(j-1)+2] += ni1*nj2*dSoverL
                        B2[2*(i-1)+2, 2*(j-1)+2] += ni2*nj2*dSoverL
                    end
                end
                if src_panel==trg_panel
                    # Add jump
                    for i=1:2*panelorder
                        B2[i,i] += 1.0/2                        
                    end
                end
            end
            # Insert into structure
            indices[idx] = src_panel -1
            data[idx,:,:] = B2
            if correction
                @. data[idx,:,:] -= B1
            end
            idx += 1
        end
        indptr[trg_panel+1] = idx -1       
    end
    # Assemble BSR matrix
    s = pyimport("scipy.sparse")
    A = BSR_Matrix(s[:bsr_matrix]((data, indices, indptr),
                                  blocksize=(blocksize, blocksize)))
    return A
end


function doublelayer_matrix(grid::DiscreteCurve,
                            alpha::Float64)
    npoints = grid.numpoints
    panelorder = grid.panelorder   
    A = Array{Float64}(2*npoints, 2*npoints) # Output matrix
    B = Array{Float64}(2*panelorder, 2*panelorder) # Panel-to-panel block
    for src_panel = 1:grid.numpanels
        for trg_panel = 1:grid.numpanels
            p2p_block!(B, grid, src_panel, trg_panel, alpha)
            block_nxn_copy!(A, trg_panel, src_panel, B, 2*panelorder)
        end
    end
    return A
end

@inline function p2p_block!(Bin, grid, src_panel, trg_panel, alpha;
                            upsampling=false,
                            specquad=true)
    panelorder = numpoints = grid.panelorder
    self_panel = (trg_panel==src_panel)
    prev_panel = (trg_panel==grid.prevpanel[src_panel])
    next_panel = (trg_panel==grid.nextpanel[src_panel])    
    if self_panel || prev_panel || next_panel
        # The upsample flag corresponds to Scheme B in Helsing & Holst,
        # which computes panel interactions using upsampled grids.
        # This improves accuracy for large alpha when subdivision are not used
        # Subdivisions make upsampling obsolete, so turned off by default
        upsample = upsampling
        near = specquad
        if self_panel
            d = 0.0
        elseif prev_panel
            d = 1.0
        else
            d = -1.0
        end        
    else
        upsample = false
        near = false
    end
    # Setup buffer block
    if upsample
        B = Array{Float64}(2*2*panelorder, 2*2*panelorder)
    else
        B = Bin;
    end    
    # Load sources
    src_idx = (1:panelorder) + panelorder*(src_panel-1)
    zsrc = grid.points[1:2, src_idx]
    nsrc = grid.normals[1:2, src_idx]
    wsrc = grid.dS[src_idx]
    # Load targets
    trg_idx = (1:panelorder) + panelorder*(trg_panel-1)
    ztrg = grid.points[1:2, trg_idx]
    if upsample
        P = upsampling_matrix(panelorder, 2*panelorder)
        Q = upsampling_matrix(2*panelorder, panelorder)
        zsrc = (P*zsrc')'
        nsrc = (P*nsrc')'
        ztrg = (P*ztrg')'
        numpoints *= 2
        # Take GL weights out of interpolation
        _, glweights1 = gausslegendre(panelorder)
        _, glweights2 = gausslegendre(2*panelorder)
        wsrc = glweights2 .* (P*(wsrc ./ glweights1))
    end
    # Compute interaction block between src and trg panels
    if near
        ## Near, need some kind of special quadrature
        # Transformations
        trg_dt = grid.t_edges[2, trg_panel] - grid.t_edges[1, trg_panel]
        src_dt = grid.t_edges[2, src_panel] - grid.t_edges[1, src_panel]
        scale = trg_dt / src_dt
        trans = -d*(1+scale)
        # If panel is short enough w.r.t. alpha, use direct computation
        # If panel is too long, compute quadrature using subdivisions of src panel        
        panel_max = LARGEALPHA_LIMIT_16/alpha
        panel_length = sum(wsrc)        
        if panel_length <= panel_max
            ## Direct
            dblkernel_block!(B, ztrg, zsrc, nsrc, wsrc, alpha)            
            add_logcorr_block!(B, ztrg, zsrc, nsrc, wsrc, scale, trans, alpha)            
            if self_panel
                # Correct 2x2 blocks on diagonal, is now NaN
                ksrc = grid.curvature[src_idx] # Curvature
                if upsample
                    ksrc = P*ksrc
                end
                for i=1:numpoints
                    ni = nsrc[:,i]
                    ti = [ni[2],-ni[1]]
                    T2C0 = -1/pi
                    K = (-T2C0*ksrc[i]*0.5*wsrc[i])*(ti*ti')
                    block_2x2_copy!(B, i, i, K)
                end                
            end                
        else
            ## Using subdivisions
            fill!(B, 0.0) # Need to initialize, since we add to it in loop
            # Seems subdivisions are fine with using 16 points
            numsec = 16
            gcpoints, gcweights = gausslegendre(numsec)            
            glpoints, glweights = gausslegendre(numpoints)
            bcweights = bclag_interp_weights(glpoints)
            for i=1:numpoints # Target point loop
                # Location of target in parametrization of source panel 
                ttrg = trans + scale*glpoints[i]
                # Allowed panel size, in parametrization                
                hsub = LARGEALPHA_LIMIT_16*2/(alpha*panel_length)
                # Don't let min subpanel size be too close to 2.0
                #hsub = min(1.9, hsub)
                # Get subdivision of source interval
                #intervals = subdivide_interval(hsub, ttrg)
                intervals = subdivide_interval_with_bisection(hsub, ttrg, numsec)
                
                R = zeros(2, 2*numsec) # Output row
                for idx = 1:length(intervals)-1
                    fill!(R, 0.0)                    
                    ta = intervals[idx]
                    tb = intervals[idx+1]
                    # New source points, in original frame
                    tsrc = ta + (1+gcpoints)*(tb-ta)/2
                    # Interpolation to new source points
                    Psec = bclag_interp_matrix(glpoints, tsrc, bcweights)
                    zsrc_sec = (Psec*zsrc')'
                    nsrc_sec = (Psec*nsrc')'
                    wsrc_sec = gcweights*(tb-ta)/2 .* (Psec*(wsrc ./ glweights))
                    # Direct double layer interaction
                    dblkernel_block!(R, ztrg[1:2,i], zsrc_sec, nsrc_sec, wsrc_sec, alpha)                    
                    # Target point, in new frame                   
                    ttrg_sec = (ttrg-ta)*2/(tb-ta) - 1.0                    
                    # tdist = Distance in parametrization of src panel                                        
                    tdist = abs.(ttrg_sec - gcpoints)
                    # Very empirical limit for when neighbors need special quad due to log kernel:
                    # 0.2 for 32, 0.4 for 16
                    # Too small gives nearly sing error
                    # Too large makes special quad inaccurate, though not very sensitive                    
                    if minimum(tdist)<0.2*32/numsec
                        # Compute Helsing log-correction and apply to row
                        WfrakL = WfrakLcomp(ttrg_sec, gcpoints)
                        KL = Array{Float64}(2,2)                        
                        for j=1:numsec
                            nb_wcorr = WfrakL[1, j]/gcweights[j] - log(tdist[j])
                            zi = ztrg[1:2, i]
                            zj = zsrc_sec[1:2, j]
                            nj = nsrc_sec[1:2, j]
                            wj = wsrc_sec[j]
                            rij = zj.-zi
                            dblkernel_log!(rij, nj, alpha, KL)                
                            block_2x2_add!(R, 1, j, KL.*(nb_wcorr*wj))
                        end
                    end
                    # Add row*interpolation to output block
                    B[2*(i-1)+(1:2), :] += R*kron(Psec, eye(2))
                end
            end
            

        end        
    else
        # Not near, direct interaction
        dblkernel_block!(B, ztrg, zsrc, nsrc, wsrc, alpha)        
    end
    if upsample
        Bin .= kron(Q, eye(2)) * B * kron(P, eye(2))
    end
end

@inline function dblkernel_block!(B, ztrg, zsrc, nsrc, wsrc, alpha)
    K = Array{Float64}(2,2)    
    rji = Array{Float64}(2)
    nj = Array{Float64}(2)
    numtrg = size(ztrg, 2)
    numsrc = size(zsrc, 2)
    for j=1:numsrc
        @. nj = nsrc[:,j]*wsrc[j]
        for i=1:numtrg
            @. rji = zsrc[:,j]-ztrg[:,i]
            dblkernel!(rji, nj, alpha, K)
            block_2x2_copy!(B, i, j, K)
        end
    end    
end

@inline function block_nxn_copy!(A, i, j, K, n)
    for k=1:n
        for l=1:n
            A[n*(i-1)+k, n*(j-1)+l] = K[k,l]
        end
    end
end

@inline function block_2x2_copy!(A, i, j, K)
    for k=1:2
        for l=1:2
            A[2*(i-1)+k, 2*(j-1)+l] = K[k,l]
        end
    end
end

@inline function block_2x2_add!(A, i, j, K)
    for k=1:2
        for l=1:2
            A[2*(i-1)+k, 2*(j-1)+l] += K[k,l]
        end
    end
end

