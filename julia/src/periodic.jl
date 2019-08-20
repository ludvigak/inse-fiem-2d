function per_modstokes(F1, F2, L, lambda, xnu=[], ynu=[]; zerok0::Bool=false, ifgrad::Bool=false)
    # Check input
    M = size(F1)
    @assert length(M)>1 && M[1]==M[2] "Grid must be square"
    @assert size(F1)==size(F2)
    N = M[1]
    # Do FFTs in-place
    Fhat1, Fhat2 = fft(F1), fft(F2)
    # Convolve with truncated greens function    
    Uhat1, Uhat2 = Fhat1, Fhat2 # Operate in-place
    k1, k2 = k_vectors([N, N], [L, L])
    k1 = ifftshift(k1)
    k2 = ifftshift(k2)
    # Store k==0 values
    Fhat1k0 = Fhat1[1,1]
    Fhat2k0 = Fhat2[1,1]
    for j=1:N
        for i=1:N
            K1 = k1[i]
            K2 = k2[j]
            Ksq = K1*K1 + K2*K2
            G = 1/ (Ksq*(Ksq + lambda^2))
            KdotFhat = K1*Fhat1[i,j] + K2*Fhat2[i,j]            
            Uhat1[i,j] = G*(Ksq*Fhat1[i,j] - K1*KdotFhat)
            Uhat2[i,j] = G*(Ksq*Fhat2[i,j] - K2*KdotFhat)            
        end
    end
    # Correct at k==0 (where just just got div by zero)
    Uhat1[1,1] = Fhat1k0 / lambda^2 
    Uhat2[1,1] = Fhat2k0 / lambda^2 

    if zerok0
        # This improves conditioning, but gets answer wrong by a constant
        # (which we can let integral equation solver pick up)
        # NOTE: It's possible that this is true only because we are using an
        # ill-conditioned test case
        info("Zeroing k=0 mode")        
        Uhat1[1,1], Uhat2[1,1] = 0.0, 0.0
    end
    
    Nnu = length(xnu)    
    if Nnu > 0
        # NUFFT for off-grid points
        @assert length(xnu)==length(ynu)
        scale = 2*pi/L
        # Rescale and shift points from [-L/2,L/2] to [0,2*pi]    
        xnu *= scale
        ynu *= scale
        xnu += pi
        ynu += pi
        # Compute        
        unu1 = Array{Complex128}(Nnu)
        unu2 = Array{Complex128}(Nnu)
        opts = finufft_default_opts()
        opts.modeord = 1
        nufft2d2!(xnu, ynu, unu1, 1, 1e-15, Uhat1, opts)
        nufft2d2!(xnu, ynu, unu2, 1, 1e-15, Uhat2, opts)
        # Get real parts and rescale values
        unu1 = real(unu1)/N^2
        unu2 = real(unu2)/N^2
    else
        unu1, unu2 = [], []
    end

    # Gradient
    if ifgrad
        U1x = zero(Uhat1)
        U1y = zero(Uhat1)
        U2x = zero(Uhat2)
        U2y = zero(Uhat2)
        for j=1:N
            for i=1:N
                U1x[i,j] = 1im*k1[i]*Uhat1[i,j]
                U1y[i,j] = 1im*k2[j]*Uhat1[i,j]
                U2x[i,j] = 1im*k1[i]*Uhat2[i,j]
                U2y[i,j] = 1im*k2[j]*Uhat2[i,j]                            
            end
        end
        ifft!(U1x)
        ifft!(U1y)        
        ifft!(U2x)
        ifft!(U2y)
        U1x, U1y, U2x, U2y = real(U1x), real(U1y), real(U2x), real(U2y)
    else
        U1x, U1y, U2x, U2y = [], [], [], []
    end
    
    # Transform back in-place
    ifft!(Uhat1)
    ifft!(Uhat2)
    # Return real parts
    U1 = real(Uhat1)
    U2 = real(Uhat2)
    return U1, U2, unu1, unu2, U1x, U1y, U2x, U2y
end

# Compute gradient of real periodic function using Fourier transform
function per_gradient(F, L)
    # Check input
    M = size(F)
    @assert length(M)>1 && M[1]==M[2] "Grid must be square"
    N = M[1]
    Fhat = fft(F)    
    k1, k2 = k_vectors([N, N], [L, L])
    k1 = ifftshift(k1)
    k2 = ifftshift(k2)
    Fx = zeros(Fhat)
    Fy = zeros(Fhat)
    for j=1:N
        for i=1:N
            K1 = k1[i]
            K2 = k2[j]
            Fx[i,j] = 1im*K1*Fhat[i,j]
            Fy[i,j] = 1im*K2*Fhat[i,j]            
        end
    end
    # Transform back in-place
    ifft!(Fx)
    ifft!(Fy)
    # Return real parts
    return real(Fx), real(Fy)
end
