# examples/cancellation.jl shows effects of these settings
const POWERSERIES_ZMAX = 1.5
const POWERSERIES_NMAX = 11
const POWERSERIES_ZMAX_PRIME = 2.0
const POWERSERIES_NMAX_PRIME = 13


# GSL implementation faster than SpecialFunctions
import GSL
# This survives large arguments better
besselK0(z) = GSL.sf_bessel_K0_scaled(z)*exp(-z) 
besselK1(z) = GSL.sf_bessel_K1_scaled(z)*exp(-z) 
besselI0(z) = GSL.sf_bessel_I0_scaled(z)*exp(z) 
besselI1(z) = GSL.sf_bessel_I1_scaled(z)*exp(z) 

# besselK0(z) = besselk(0, z)
# besselK1(z) = besselk(1, z)
# besselI0(z) = besseli(0, z)
# besselI1(z) = besseli(1, z)

# ========== STOKESLET
"Computes S_{ij} f_j"
function stokeslet(r, f, alpha)
    rsq = r[1]*r[1] + r[2]*r[2]
    zsq = rsq*alpha*alpha
    z = sqrt(zsq)
    K0 = besselK0(z)
    K1 = besselK1(z)
    S1 = (zsq*K0 + z*K1 - 1) / (2*pi*zsq)
    S2 = -(zsq*K0 + 2*z*K1 - 2) / (2*pi*zsq*zsq)
    rdotf = (r[1]*f[1] + r[2]*f[2])
    return f.*S1 .+ r.*(alpha*alpha*S2*rdotf)
end

function stokeslet_kernel_noarrays(r1, r2, alpha)
    rsq = r1*r1 + r2*r2
    zsq = rsq*alpha*alpha
    z = sqrt(zsq)
    K0 = besselK0(z)
    K1 = besselK1(z)
    S1 = (zsq*K0 + z*K1 - 1) / (2*pi*zsq)
    S2 = -(zsq*K0 + 2*z*K1 - 2) / (2*pi*zsq*zsq)
    A11 = S1 + alpha^2*S2*r1*r1
    A12 = alpha^2*S2*r1*r2
    A21 = A12
    A22 = S1 + alpha^2*S2*r2*r2
    return A11, A12, A21, A22
end

# ========== PRESSURE(LET)
"Computes p_i f_i"
function pressure(r, f, alpha)
    rsq = r[1]*r[1] + r[2]*r[2]
    rdotf = (r[1]*f[1] + r[2]*f[2])
    return rdotf/(2*pi*rsq)
end

# ========= STRESSLET
"Computes T_{ijk} f_i n_k"
function stresslet(r::Array{Float64}, f::Array{Float64},
                   n::Array{Float64}, alpha::Float64)
    rsq = r[1]*r[1] + r[2]*r[2]
    rnorm = sqrt(rsq)
    z = rnorm*alpha
    (T1, T2, T3) = Ti(z)
    
    rdotf = r[1]*f[1] + r[2]*f[2]
    rdotn = r[1]*n[1] + r[2]*n[2]
    fdotn = f[1]*n[1] + f[2]*n[2]

    alpha2 = alpha*alpha
    alpha4 = alpha2*alpha2
    
    u1 = alpha2*T1*(r*fdotn + f*rdotn + n*rdotf)
    u2 = alpha4*T2*rdotf*rdotn*r
    u3 = alpha2*T3*r*fdotn;
    
    return u1 + u2 + u3
end

"Computes K_{ji} = T_{ijk} n_k"
function dblkernel(r::Array{Float64},
                   n::Array{Float64},
                   alpha::Float64 )
    K = Array{Float64}(undef, 2,2)
    dblkernel!(r, n, alpha, K)
    return K
end

function dblkernel!(r::Array{Float64},
                    n::Array{Float64},
                    alpha::Float64,
                    K::Array{Float64})
    K11, K12, K21, K22 = dblkernel_noarrays(r[1], r[2], n[1], n[2], alpha)
    K[1,1] = K11
    K[1,2] = K12
    K[2,1] = K21
    K[2,2] = K22
end

# Compute kernel without passing around arrays
function dblkernel_noarrays(r1, r2, n1, n2, alpha)
    rsq = r1*r1 + r2*r2
    rnorm = sqrt(rsq)
    z = rnorm*alpha
    (T1, T2, T3) = Ti(z)
    alpha2 = alpha*alpha
    alpha4 = alpha2*alpha2
    T1 *= alpha2
    T2 *= alpha4
    T3 *= alpha2
    rdotn = r1*n1 + r2*n2
    K11 = T1*( n1*r1 + n1*r1 + rdotn ) +
        T2*r1*r1*rdotn +
        T3*n1*r1
    K21 = T1*( n2*r1 + n1*r2  ) +
        T2*r1*r2*rdotn +
        T3*n1*r2    
    K12 = T1*( n1*r2 + n2*r1  ) +
        T2*r2*r1*rdotn +
        T3*n2*r1
    K22 = T1*( n2*r2 + n2*r2 + rdotn ) +
        T2*r2*r2*rdotn +
        T3*n2*r2
    K11, K12, K21, K22
end


"""
Computes K_^L{ji}=T_^L{ijk} n_k
# Returns: KL
"""
function dblkernel_log(r::Array{Float64},
                       n::Array{Float64},
                       alpha::Float64 )
    KL = Array{Float64}(undef, 2,2)
    dblkernel_log!(r, n, alpha, KL)
    return KL
end

function dblkernel_log!(r::Array{Float64},
                        n::Array{Float64},
                        alpha::Float64,
                        KL::Array{Float64,2} )
    rsq = r[1]*r[1] + r[2]*r[2]
    rnorm = sqrt(rsq)
    z = rnorm*alpha
    (T1S, T1L, T2S, T2L, T3S, T3L) = Ti_split(z)
    alpha2 = alpha*alpha
    alpha4 = alpha2*alpha2
    T1L *= alpha2
    T2L *= alpha4
    T3L *= alpha2    
    rdotn = r[1]*n[1] + r[2]*n[2]
    for i=1:2
        for j=1:2
            KL[j + 2*(i-1)] = 
                T1L*( n[j]*r[i] + n[i]*r[j] + (i==j)*rdotn ) +
                T2L*r[i]*r[j]*rdotn +
                T3L*n[i]*r[j]
        end
    end
end

"""
Computes T_{ijk} f_i n_k, but with kernel-splitting
# Returns: (uS, uL, uC, uQ)
such that 
u = uS + uL*log(r) + uC*r_k*n_k/r^2 + uQ*f_i*r_i*r_j*r_k*n_k/r^4
"""
function stresslet_split(r, f, n, alpha)
    rnorm = sqrt(r[1]*r[1] + r[2]*r[2])
    z = rnorm*alpha
    (T1S, T1L, T2S, T2L, T3S, T3L) = Ti_split(z)
    rdotf = r[1]*f[1] + r[2]*f[2]
    rdotn = r[1]*n[1] + r[2]*n[2]
    fdotn = f[1]*n[1] + f[2]*n[2]

    vec1 = r*fdotn + f*rdotn + n*rdotf
    vec2 = r*(rdotf*rdotn)
    vec3 = r*fdotn

    alpha2 = alpha*alpha
    alpha4 = alpha2*alpha2    
    
    uS = alpha2*(T1S*vec1 + alpha2*T2S*vec2 + T3S*vec3)
    uL = alpha2*(T1L*vec1 + alpha2*T2L*vec2 + T3L*vec3) 
    uS = uS + log(alpha)*uL # Split on log(r)
    uC = alpha2/(8*pi)*r*rdotf
    uQ = -1/pi
    
    return (uS, uL, uC, uQ)    
end

function stresslet_split_noarrays(r1, r2, f1, f2, n1, n2, alpha)
    rnorm = sqrt(r1*r1 + r2*r2)
    z = rnorm*alpha
    (T1S, T1L, T2S, T2L, T3S, T3L) = Ti_split(z)
    rdotf = r1*f1 + r2*f2
    rdotn = r1*n1 + r2*n2
    fdotn = f1*n1 + f2*n2

    vecA1 = r1*fdotn + f1*rdotn + n1*rdotf
    vecA2 = r2*fdotn + f2*rdotn + n2*rdotf
    vecB1 = r1*(rdotf*rdotn)
    vecB2 = r2*(rdotf*rdotn)
    vecC1 = r1*fdotn
    vecC2 = r2*fdotn

    alpha2 = alpha*alpha
    alpha4 = alpha2*alpha2    
    
    uS1 = alpha2*(T1S*vecA1 + alpha2*T2S*vecB1 + T3S*vecC1)
    uS2 = alpha2*(T1S*vecA2 + alpha2*T2S*vecB2 + T3S*vecC2)    
    uL1 = alpha2*(T1L*vecA1 + alpha2*T2L*vecB1 + T3L*vecC1)
    uL2 = alpha2*(T1L*vecA2 + alpha2*T2L*vecB2 + T3L*vecC2)     
    uS1 += log(alpha)*uL1 # Split on log(r)
    uS2 += log(alpha)*uL2 # Split on log(r)    
    uC1 = alpha2/(8*pi)*r1*rdotf
    uC2 = alpha2/(8*pi)*r2*rdotf    
    uQ = -1/pi
    
    return (uS1, uS2, uL1, uL2, uC1, uC2, uQ)    
end

"""
Stresslet functions
# Returns: (T1, T2, T3)
"""
function Ti(z::Float64)
    if z < POWERSERIES_ZMAX
        (T1S, T1L, T2S, T2L, T3S, T3L) = Ti_split_pow(z)
        T1 = T1S + T1L*log(z)
        T2 = T2S + T2L*log(z) + 1/(8*pi*z^2) - 1/(pi*z^4)
        T3 = T3S + T3L*log(z)
        return (T1, T2, T3)        
    else
        return Ti_direct(z)
    end
end

function Ti_direct(z::Float64)
    zsq = z*z
    z4 = zsq*zsq
    z6 = z4*zsq
    K0 = besselK0(z)
    K1 = besselK1(z)
    T1 = -(2*zsq*K0 + (zsq+4)*z*K1 - 4) / (2*pi*z4)
    T2 = (4*zsq*K0 + (zsq+8)*z*K1 - 8) / (pi*z6)
    T3 = (z*K1 - 1) / (2*pi*zsq)
    return (T1, T2, T3)    
end

function Ti_prime(z::Float64)
    if z < POWERSERIES_ZMAX_PRIME
        (T1S, T1L, T2S, T2L, T3S, T3L) = Ti_split_pow(z, false)
        (DT1S, DT1L, DT2S, DT2L, DT3S, DT3L) = Ti_split_pow(z, true)
        T1p = DT1S + DT1L*log(z) + T1L / z
        T2p = DT2S + DT2L*log(z) + T2L/z - 2/(8*pi*z^3) + 4/(pi*z^5)
        T3p = DT3S + DT3L*log(z) + T3L/z
        return (T1p, T2p, T3p)        
    else
        return Ti_direct_prime(z)
    end
end

function Ti_direct_prime(z::Float64)
    zsq = z*z
    z4 = zsq*zsq
    z6 = z4*zsq
    K0 = besselK0(z)
    K1 = besselK1(z)
    T1p = ((2*z^2 + 8)*K1)/(z^4*pi) - 8/(z^5*pi) + ((z^2 + 8)*K0)/(2*z^3*pi)
    T2p = 48/(z^7*pi) - ((8*z^2 + 48)*K1)/(z^6*pi) - ((z^2 + 24)*K0)/(z^5*pi)
    T3p = 1/(z^3*pi) - K0/(2*z*pi) - K1/(z^2*pi)
    return (T1p, T2p, T3p)    
end

"""
Split of stresslet functions
# Returns: (T1S, T1L, T2S, T2L, T3S, T3L)
"""
function Ti_split(z::Float64)
    if z < POWERSERIES_ZMAX
        return Ti_split_pow(z)
    else
        return Ti_split_direct(z) 
    end    
end

function Ti_split_prime(z::Float64)
    if z < POWERSERIES_ZMAX_PRIME
        return Ti_split_pow(z,true)
    else
        return Ti_split_direct(z,true)
    end    
end

function Ti_split_direct(z::Float64,ifprime::Bool=false)

    if !ifprime

        K0 = besselK0(z)
        K1 = besselK1(z)
        I0 = besselI0(z)
        I1 = besselI1(z)
        
        K0S = K0 + I0*log(z)
        K1S = K1 - I1*log(z) - 1/z    

        T1S = -(2*z*K0S + (z^2+4)*K1S + z) / (2*pi*z^3)
        T1L = (2*z*I0 - (z^2+4)*I1) / (2*pi*z^3)

        T2S = (32*z*K0S + 8*(z^2+8)*K1S  - z*(z^2-16)) / (8*pi*z^5)
        T2L = ((z^2+8)*I1 - 4*z*I0) / (pi*z^5)

        T3S = K1S/(2*pi*z)
        T3L = I1/(2*pi*z)
        return (T1S, T1L, T2S, T2L, T3S, T3L)

    else

        logz = log(z)
        z2 = z*z
        z3 = z*z2
        z4 = z*z3
        z5 = z*z4
        z6 = z*z5
        z7 = z*z6
        
        K0 = besselK0(z)
        K1 = besselK1(z)        
        I0 = besselI0(z)
        I1 = besselI1(z)
        
        DT1S = (-16 +
                z*(z*(8+z2)*K0+4*(4+z2)*K1-(4+z2)*I1*(4*logz-1)+
                   z*I0*((8+z2)*logz-2)))/(2*pi*z5)
        DT2S = (192-16*z2+z4+
                4*z*(-z*(24+z2)*K0-8*(6+z2)*K1+
                     I1*(-8-z2+8*(6+z2)*logz)-
                     z*I0*(-4+(24+z2)*logz)))/(4*pi*z7)
        DT3S = (2-I1*(z-2*z*logz)-
                z2*(2*K1/z+K0+I0*logz))/(2*pi*z3)

        DT1L = (-2*z*(8 + z2)*I0 + 8*(4 + z2)*I1)/(4*pi*z4)
        DT2L = (z*(24 + z2)*I0 - 8*(6 + z2)*I1)/(pi*z6)
        DT3L = (I0-2*I1/z)/(2*pi*z)
        
        return (DT1S, DT1L, DT2S, DT2L, DT3S, DT3L)

    end

end

"""
Direct power series for split stresslet functions
# Returns: (T1S, T1L, T2S, T2L, T3S, T3L)
"""
function Ti_split_pow(z::Float64,ifprime::Bool=false)
    if !ifprime
        # Outputs
        T1S = 0.0
        T1L = 0.0
        T2S = 0.0
        T2L = 0.0
        T3S = 0.0
        T3L = 0.0
        # Add n=0 term
        log2 = log(2)    
        k1L0 = 0.5
        k1S0 = (-2*log2-(-2*eulergamma+1.0)) * 0.25
        T3S += k1S0
        T3L += k1L0
        # Recursion variables, init at n=1
        n=1
        z2 = z*z
        z2nm2 = 1; # z^(2*n-2)    
        z2n = z2; # z^(2*n)    
        digamma_np1 = -eulergamma + 1.0 # digamma(n+1)
        digamma_np2 = -eulergamma + 1.5
        ifact_n = 1.0 # 1/factorial(n)
        ifact_np1 = 1.0/2
        i4n = 1.0/4 # 1/4^n
        i4np1 = 1.0/16
        k1Snm1 = k1S0 # k1S(n-1)
        k1Lnm1 = k1L0 # k1L(n-1)
        ## First iteration unrolled
        # k0, k1
        k0n_denom = i4n*ifact_n*ifact_n
        k1n_denom = i4np1*ifact_n*ifact_np1
        k0Ln = (-1)*k0n_denom
        k0Sn = (digamma_np1 + log2)*k0n_denom
        k1Ln = (2)*k1n_denom
        k1Sn = (-2*log2-digamma_np1-digamma_np2)*k1n_denom
        # T1
        t1Sn = 2*k0Sn + 4*k1Sn + k1Snm1
        t1Ln = 2*k0Ln + 4*k1Ln + k1Lnm1
        T1S += t1Sn*z2nm2
        T1L += t1Ln*z2nm2
        # T3
        T3S += k1Sn*z2n
        T3L += k1Ln*z2n  
        # Update recursions
        n += 1
        inp1 = 1/(n+1)
        z2n *= z2
        z2nm2 = z2
        z2nm4 = 1
        (digamma_np1, digamma_np2) = (digamma_np2, digamma_np2+inp1)
        (ifact_n, ifact_np1) = (ifact_np1, ifact_np1*inp1)
        (i4n, i4np1) = (i4np1, i4np1*0.25)
        k1Snm1 = k1Sn
        k1Lnm1 = k1Ln     
        ## Loop
        while n<POWERSERIES_NMAX
            # k0, k1
            k0n_denom = i4n*ifact_n*ifact_n
            k1n_denom = i4np1*ifact_n*ifact_np1
            k0Ln = (-1)*k0n_denom
            k0Sn = (digamma_np1 + log2)*k0n_denom
            k1Ln = (2)*k1n_denom
            k1Sn = (-2*log2-digamma_np1-digamma_np2)*k1n_denom
            # T1
            t1Sn = 2*k0Sn + 4*k1Sn + k1Snm1
            t1Ln = 2*k0Ln + 4*k1Ln + k1Lnm1
            T1S += t1Sn*z2nm2
            T1L += t1Ln*z2nm2
            # T2
            t2Sn = 4*k0Sn + 8*k1Sn + k1Snm1
            t2Ln = 4*k0Ln + 8*k1Ln + k1Lnm1
            T2S += t2Sn*z2nm4
            T2L += t2Ln*z2nm4
            # T3
            T3S += k1Sn*z2n
            T3L += k1Ln*z2n  
            # Update recursions
            n += 1
            inp1 = 1/(n+1)
            (z2nm4, z2nm2, z2n) = (z2nm2, z2n, z2n*z2)
            (digamma_np1, digamma_np2) = (digamma_np2, digamma_np2+inp1)
            (ifact_n, ifact_np1) = (ifact_np1, ifact_np1*inp1)
            (i4n, i4np1) = (i4np1, i4np1*0.25)
            k1Snm1 = k1Sn
            k1Lnm1 = k1Ln        
        end
        # Multiply with constant factors
        T1S *= -1/(2*pi)
        T1L *= -1/(2*pi)    
        T2S *= 1/pi
        T2L *= 1/pi
        T3S *= 1/(2*pi)
        T3L *= 1/(2*pi)
        return (T1S, T1L, T2S, T2L, T3S, T3L)
    else
        # Outputs
        T1S = 0.0
        T1L = 0.0
        T2S = 0.0
        T2L = 0.0
        T3S = 0.0
        T3L = 0.0
        # Add n=0 term (derivative is zero)
        log2 = log(2)    
        k1L0 = 0.5
        k1S0 = (-2*log2-(-2*eulergamma+1.0)) * 0.25
        #T3S += k1S0
        #T3L += k1L0
        # Recursion variables, init at n=1
        n=1
        z2 = z*z
        z2nm3 = 0; # z^(2*n-3)*(2n-2)
        z2nm1 = z*2; # z^(2*n-1)*2n   
        digamma_np1 = -eulergamma + 1.0 # digamma(n+1)
        digamma_np2 = -eulergamma + 1.5
        ifact_n = 1.0 # 1/factorial(n)
        ifact_np1 = 1.0/2
        i4n = 1.0/4 # 1/4^n
        i4np1 = 1.0/16
        k1Snm1 = k1S0 # k1S(n-1)
        k1Lnm1 = k1L0 # k1L(n-1)
        ## First iteration unrolled
        # k0, k1
        k0n_denom = i4n*ifact_n*ifact_n
        k1n_denom = i4np1*ifact_n*ifact_np1
        k0Ln = (-1)*k0n_denom
        k0Sn = (digamma_np1 + log2)*k0n_denom
        k1Ln = (2)*k1n_denom
        k1Sn = (-2*log2-digamma_np1-digamma_np2)*k1n_denom
        # T1 (constant ... derivative is zero)
        #t1Sn = 2*k0Sn + 4*k1Sn + k1Snm1
        #t1Ln = 2*k0Ln + 4*k1Ln + k1Lnm1
        #T1S += t1Sn*z2nm3
        #T1L += t1Ln*z2nm3
        # T3
        T3S += k1Sn*z2nm1
        T3L += k1Ln*z2nm1  
        # Update recursions
        n += 1
        inp1 = 1/(n+1)
        z2nm1 *= z2*(n/(n-1.0))
        z2nm3 = z*2
        z2nm5 = 0
        (digamma_np1, digamma_np2) = (digamma_np2, digamma_np2+inp1)
        (ifact_n, ifact_np1) = (ifact_np1, ifact_np1*inp1)
        (i4n, i4np1) = (i4np1, i4np1*0.25)
        k1Snm1 = k1Sn
        k1Lnm1 = k1Ln     
        ## Loop
        while n<POWERSERIES_NMAX_PRIME
            # k0, k1
            k0n_denom = i4n*ifact_n*ifact_n
            k1n_denom = i4np1*ifact_n*ifact_np1
            k0Ln = (-1)*k0n_denom
            k0Sn = (digamma_np1 + log2)*k0n_denom
            k1Ln = (2)*k1n_denom
            k1Sn = (-2*log2-digamma_np1-digamma_np2)*k1n_denom
            # T1
            t1Sn = 2*k0Sn + 4*k1Sn + k1Snm1
            t1Ln = 2*k0Ln + 4*k1Ln + k1Lnm1
            T1S += t1Sn*z2nm3
            T1L += t1Ln*z2nm3
            # T2
            t2Sn = 4*k0Sn + 8*k1Sn + k1Snm1
            t2Ln = 4*k0Ln + 8*k1Ln + k1Lnm1
            T2S += t2Sn*z2nm5
            T2L += t2Ln*z2nm5
            # T3
            T3S += k1Sn*z2nm1
            T3L += k1Ln*z2nm1  
            # Update recursions
            n += 1
            inp1 = 1/(n+1)
            (z2nm5, z2nm3, z2nm1) = (z2nm3, z2nm1, z2nm1*z2*(n/(n-1.0)))
            (digamma_np1, digamma_np2) = (digamma_np2, digamma_np2+inp1)
            (ifact_n, ifact_np1) = (ifact_np1, ifact_np1*inp1)
            (i4n, i4np1) = (i4np1, i4np1*0.25)
            k1Snm1 = k1Sn
            k1Lnm1 = k1Ln        
        end
# Multiply with constant factors
        T1S *= -1/(2*pi)
        T1L *= -1/(2*pi)    
        T2S *= 1/pi
        T2L *= 1/pi
        T3S *= 1/(2*pi)
        T3L *= 1/(2*pi)
        return (T1S, T1L, T2S, T2L, T3S, T3L)
    end
end
