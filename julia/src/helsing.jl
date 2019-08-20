#
# Routines for modified quadrature weights originally written by Johan Helsing
#

using FastGaussQuadrature
using Vandermonde

@inline function add_logcorr_block!(B, ztrg, zsrc, nsrc, wsrc, scale, trans, alpha)
    numpoints = length(wsrc)
    glpoints, glweights = gausslegendre(numpoints)        
    if trans==0
        glpoints_scaled = glpoints            
    else
        glpoints_scaled = @. trans + scale*glpoints
    end
    #WfrakL = WfrakLcomp(trans, scale, glpoints)
    WfrakL = WfrakLcomp(glpoints_scaled, glpoints)
    # Assemble block (including diagonal)
    KL = Array{Float64}(2,2)
    tdist = Array{Float64}(numpoints)
    for i=1:numpoints
        # tdist = Distance in parametrization of src panel
        @. tdist = abs(glpoints_scaled[i] - glpoints[:])
        # Very empirical limit for when neighbors need special quad due to log kernel:
        # 0.2 for 32, 0.4 for 16
        # Too gives nearly sing error
        # Too large makes special quad inaccurate, though not very sensitive
        if (trans == 0) || (minimum(tdist) < 0.2*32/numpoints) 
            for j=1:numpoints            
                nb_wcorr = WfrakL[i, j]/glweights[j] - log(tdist[j])
                # i==j term on self panel omitted since GL(0)=0 for modified stresslet
                zi = ztrg[1:2, i]
                zj = zsrc[1:2, j]
                nj = nsrc[1:2, j]
                wj = wsrc[j]
                rij = zj.-zi
                dblkernel_log!(rij, nj, alpha, KL)                
                block_2x2_add!(B, i, j, KL.*(nb_wcorr*wj))
            end
        end
    end
end

"""
The choice of input arguments trans=0 and scale=1 gives WL
for ri and rj on the same quadrature panel γp. The choice
trans=±2 and scale=1 gives WL for ri on a neighboring panel
γp±1, assuming it is equal in parameter length. The input
argument tfrak is a column vector whose entries contain the
canonical nodes ti.
"""
function WfrakLcomp(trans::Float64,
                    scale::Float64,
                    tfrak::Array{Float64,1})
    npt = length(tfrak)
    tt=trans+scale*tfrak
    Q=zeros(npt, npt)
    p=zeros(1,npt+1)
    for j=1:npt
        p[1]=log(abs((1-tt[j])/(1+tt[j])))
        for k=1:npt
            c=(1-(-1)^k)/k
            p[k+1]=p[k]*tt[j]+c
        end
        for i=1:2:npt-1
            Q[i,j]=log(abs(1-tt[j]^2))-p[i+1] # 1658880
            Q[i+1,j]=p[1]-p[i+2]
        end
        for i=1:npt
            Q[i,j]=Q[i,j]/i
        end
    end
    
    # Direct solution:
    # A=ones(npt, npt)
    # for k=2:npt
    #     A[k,:]=tfrak.*A[k-1,:]
    # end    
    # WfrakL =(A\Q).'
    
    # Bjorck-Pereyra:
    WfrakL = zeros(npt, npt);
    for i=1:npt
        WfrakL[i,:] = pvand(tfrak, Q[:,i]);
    end
    
    return WfrakL
end

function WfrakLcomp(targets, sources)
    tt = targets
    tfrak = sources
    nsrc = length(sources)
    ntrg = length(targets)
    Q=zeros(nsrc, ntrg)
    p=zeros(1,nsrc+1)
    for j=1:ntrg
        p[1]=log(abs((1-tt[j])/(1+tt[j])))
        for k=1:nsrc
            c=(1-(-1)^k)/k
            p[k+1]=p[k]*tt[j]+c
        end
        for i=1:2:nsrc-1
            Q[i,j]=log(abs(1-tt[j]^2))-p[i+1] # 1658880
            Q[i+1,j]=p[1]-p[i+2]
        end
        for i=1:nsrc
            Q[i,j]=Q[i,j]/i
        end
    end
    # Direct solution:
    # A=ones(npt, npt)
    # for k=2:npt
    #     A[k,:]=tfrak.*A[k-1,:]
    # end    
    # WfrakL =(A\Q).'
    
    # Bjorck-Pereyra:
    WfrakL = zeros(ntrg, nsrc);
    for i=1:ntrg
        WfrakL[i,:] = pvand(tfrak, Q[:,i]);
    end
    
    return WfrakL
end

"""
    This function relies on complex arithmetic and takes, as
    input, points and vec- tors in R2 represented as points in
    C. Otherwise the notation follows Section 6: input
    parameters za and zb correspond to tau(ta) and tau(tb); z is the
    target point z ∈ E;and zj, nuj,and zpwj are column vectors
    whose entries contain the points zj , the exterior unit
    normals nu at rj , and the weighted velocity function zpwj
    , j = 1,...,npt.
    """
function modifiedweights(za::Complex{Float64},
                         zb::Complex{Float64},
                         z::Complex{Float64},
                         zj::Array{Complex{Float64}}
                         )
    npt = length(zj)
    # Transform endpoints to [-1,1]
    dz = (zb-za)/2
    ztr = (z-(zb+za)/2)/dz
    zjtr = (zj-(zb+za)/2)/dz
    p = Array{Complex{Float64}}(npt+1)
    q = Array{Complex{Float64}}(npt)
    r = Array{Complex{Float64}}(npt)
    s = Array{Complex{Float64}}(npt)
    
    c = (1-(-1).^(1:npt))./(1:npt)
    p[1] = log(1-ztr) - log(-1-ztr)
    p1 = log(1-ztr) + log(-1-ztr)
    # Signs here depend on how we have defined
    # direction of curve and normals. I think.
    if imag(ztr)<0 && abs(real(ztr))<1
        p[1] += 2im*pi
        p1 -= 2im*pi
    end
    for k=1:npt
	p[k+1] = ztr*p[k] + c[k]
    end

    r[1] = -1/(1+ztr)-1/(1-ztr)
    s[1] = -1/2*(1/(1-ztr)^2 - 1/(-1-ztr)^2)
    for k=1:npt-1
	r[k+1] = ztr*r[k] + p[k]
        s[k+1] = ztr*s[k] + r[k]
    end
    
    q[1:2:npt-1] = p1 - p[2:2:npt]
    q[2:2:npt] = p[1] - p[3:2:npt+1]
    q = q./(1:npt)
    # Bjorck-Pereyra solve
    wL = pvand(zjtr, q)*dz
    wC = pvand(zjtr, p[1:npt])
    wQ = pvand(zjtr, r) / dz
    wT = pvand(zjtr, s) / dz^2
    
    # Direct solve:
    # A = Array{Complex{Float64}}(npt, npt)
    # A[1,:] = 1.0
    # for k=2:npt
    #     for j=1:npt
    #         A[k,j]=zjtr[j]*A[k-1,j] # vandermonde transposed
    #     end
    # end    
    # wL = (A\q)*dz
    # wC = A\p[1:npt]
    # wQ = (A\r) / dz        
    # Only do one solve
    # sol = A\[q p[1:npt] r];    
    # wL = sol[:,1]*dz
    # wC = sol[:,2]
    # wQ = sol[:,3] / dz
    
    return (wL, wC, wQ, wT)
end

# Compute near eval weights for log, Cauchy and stresslet singularities
function near_eval_weights(za, zb, z, zj, nuj, zpwj)
    # Compute modified weights
    wL, wC, wQ = modifiedweights(za, zb, z, zj)
    # Create corrections
    return near_eval_weights(za, zb, z, zj, nuj, zpwj, wL, wC, wQ)
end

function near_eval_weights(za, zb, z, zj, nuj, zpwj, wL, wC, wQ)
    # Create corrections
    dz = (zb-za)/2
    nuj = -nuj # Helsing has n = -i*t
    ## Log
    wcorrL = @. imag(wL*conj(nuj)) + log(abs(dz))*abs(zpwj)
    ## Cauchy
    wcorrC = @. -imag(wC)
    ## Unmodified stresslet
    wcorrQ1 = -1/2*imag(wC)
    wcorrQ2 = @. 1/4im*conj(conj(nuj)^2*wC + conj(zj-z)*wQ)
    return wcorrL, wcorrC, wcorrQ1, wcorrQ2
end

