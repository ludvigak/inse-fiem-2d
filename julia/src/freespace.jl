include("../../external/mbh2dfmm/julia/modbhrouts.jl")

using FINUFFT

## FREESPACE SOLVER WITH NUFFT
"
fs_modstokes(F1, F2, L, lambda, xnu, ynu; precomp=GR)

Solves modified Stokes equation on NxN grid of length L with rhs with
freespace boundary conditions and compactly supported rhs (F1, F2).

GR is precomputed truncated Green's function to modified biharmonic:
GR = fs_modbh_precomp(N, L, lambda)
"
function fs_modstokes(F1, F2, L, lambda, xnu=[], ynu=[]; precomp=[])
    # Check input
    M = size(F1)
    @assert length(M)>1 && M[1]==M[2] "Grid must be square"
    @assert size(F1)==size(F2)
    N = M[1]
    # Get truncated Green's function
    if isempty(precomp)
        GRpre = fs_modbh_precomp(N, L, lambda)
    else
        GRpre = precomp
    end
    # Pad input to double size
    sf = 2
    N2 = sf*N
    info("FFT size is $(N2)x$(N2)")    
    Fpad1 = zeros(Complex128, N2, N2)
    Fpad2 = zeros(Complex128, N2, N2)
    Fpad1[1:N, 1:N] = F1
    Fpad2[1:N, 1:N] = F2
    # Do FFTs in-place
    fft!(Fpad1)
    fft!(Fpad2)
    Fhat1, Fhat2 = Fpad1, Fpad2
    # Convolve with truncated greens function    
    Uhat1, Uhat2 = Fhat1, Fhat2 # Operate in-place
    k1, k2 = k_vectors([N2, N2], [L*sf, L*sf])
    k1 = ifftshift(k1)
    k2 = ifftshift(k2)   
    k1sq, k2sq = k1.^2, k2.^2
    for j=1:N2
        for i=1:N2
            K1 = k1[i]
            K2 = k2[j]
            Ksq = k1sq[i] + k2sq[j]
            KdotFhat = K1*Fhat1[i,j] + K2*Fhat2[i,j]            
            Uhat1[i,j] = GRpre[i,j]*(Ksq*Fhat1[i,j] - K1*KdotFhat)
            Uhat2[i,j] = GRpre[i,j]*(Ksq*Fhat2[i,j] - K2*KdotFhat)            
        end
    end
    if length(xnu) > 0
        # NUFFT for off-grid points
        @assert length(xnu)==length(ynu)
        # TODO: Get this scaling right so we can use modeord=0
        origin = pi/sf
        scale = 2*pi/(sf*L)
        xn = origin + xnu*scale
        yn = origin + ynu*scale
        Nnu = length(xnu)
        unu1 = Array{Complex128}(Nnu)
        unu2 = Array{Complex128}(Nnu)
        opts = finufft_default_opts()
        opts.modeord = 1        
        nufft2d2!(xn, yn, unu1, 1, 1e-15, Uhat1, opts)
        nufft2d2!(xn, yn, unu2, 1, 1e-15, Uhat2, opts)
        uscale = 1/(N*sf)^2        
        unu1 = real(unu1)*uscale
        unu2 = real(unu2)*uscale
    else
        unu1, unu2 = [], []
    end
    # Transform back
    ifft!(Uhat1)
    ifft!(Uhat2)
    # Truncate results to original grid size
    U1 = zeros(N,N)
    U2 = zeros(N,N)
    for j=1:N
        for i=1:N
            U1[i,j] = real(Uhat1[i,j])
            U2[i,j] = real(Uhat2[i,j])
        end
    end  
    return U1, U2, unu1, unu2
end

"Precompute truncated Green's function to modified biharmonic equation"
function fs_modbh_precomp(N, L, lambda)
    sf = 1+sqrt(2) # Minimum oversampling in 2D    
    Nf = ceil(Int, sf*N)
    if mod(Nf,2) != 0
        Nf+=1 # Even grids are always faster (right?)
    end
    info("Precompute FFT size is $(Nf)x$(Nf)")
    sf = Nf/N # Effective oversampling
    R = sqrt(2)*L
    logR = log(R)
    rscale = 1.0
    n = 1
    qn = modbh_qn([R],lambda,rscale,n)
    Q0 = qn[1,1]
    Q1 = qn[2,1]
    lam2 = lambda*lambda    
    # Compute
    GR = Array{Complex128}(Nf, Nf)
    k1, k2 = k_vectors([Nf, Nf], [L*sf, L*sf])
    k1 = ifftshift(k1)
    k2 = ifftshift(k2)
    k1sq, k2sq = k1.^2, k2.^2    
    for j=1:Nf
        for i=1:Nf
            K1 = k1[i]
            K2 = k2[j]
            Ksq = k1sq[i] + k2sq[j]
            K = sqrt(Ksq)
            J1 = besselj(1, K*R)
            J0M1 = besselj0m1(K*R)
            J0 = J0M1+1
            GR[i,j] = lam2*(J0M1 + K*R*logR*J1)/(Ksq*(lam2+Ksq))
            GR[i,j] += (Q0*K*R*J1-lambda*R*Q1*J0)/(lam2+Ksq)
        end
    end
    # Correct at k=0
    mid = 1 # grid is even    
    GR[mid, mid] = R*R*(-0.25+0.5*logR)- R*Q1/lambda
    GR .*= 1/lam2
    ifft!(GR)    
    # GR now in real space and has rubbish in the center,
    # Truncate by picking out corner values
    GRtrunc = Array{Complex128}(2*N, 2*N)    
    idx = 1:N
    offset = Nf - N
    for i=idx
        for j=idx
            GRtrunc[i,  j  ] = GR[i,       j       ]
            GRtrunc[i+N,j  ] = GR[i+offset,j       ]
            GRtrunc[i,  j+N] = GR[i,       j+offset]
            GRtrunc[i+N,j+N] = GR[i+offset,j+offset]
        end
    end
    # Return Fourier transform    
    fft!(GRtrunc)
    return GRtrunc
end

## PJ FUNCTIONS
function GR_purejulia(K, R)
    J0 = besselj(0, K*R)
    J1 = besselj(1, K*R)
    K0 = besselk(0, lambda*R)
    K1 = besselk(1, lambda*R)
    GR = -( J0*(lambda*K^2*R*K1 - K^2 - lambda^2) -
            K*R*J1*((K^2+lambda^2)*log(R) + K^2*K0) + lambda^2
            ) / (K^2*(lambda^2+K^2))
    return GR
end    

function GR0_purejulia(R)
    GR = -(lambda^2*R^2 - 2*lambda^2*R^2*log(R) +
           4*lambda*R*besselk(1, R*lambda) - 4)/(4*lambda^2)
    return GR
end

## HELPERS
function k_vectors(M,box)
    function k_vec(M,L)
        if mod(M,2)==0
            MM = M/2;
            k = (2*pi/L)*collect(-MM:(MM-1));
        elseif mod(M-1,2)==0
            MM = (M-1)/2;
            k = (2*pi/L)*collect(-MM:MM);
        else
            error("k-vectors not computed");
        end
        return k
    end
    k1 = k_vec(M[1], box[1])
    k2 = k_vec(M[2], box[2])    
    return k1, k2
end
