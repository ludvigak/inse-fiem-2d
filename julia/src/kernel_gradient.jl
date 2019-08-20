

# Weights for direct gradient eval
function gradient_kernel_direct_weights(tau, dtau, nu, z, alpha)
    N = length(tau)
    DUDX = complex(zeros(2, N))
    DUDY = complex(zeros(2, N))    
    for i=1:N
        # This term is done in primary vars
        r1 = real(tau[i]-z)
        r2 = imag(tau[i]-z)
        n1 = real(nu[i])
        n2 = imag(nu[i])        
        # Kernel functions
        rnorm = sqrt(r1^2 + r2^2)
        T1S, T2S, T3S = Ti(alpha*rnorm)
        D1S, D2S, D3S = Ti_prime(alpha*rnorm)        
        t1s = T1S*alpha^2
        t2s = T2S*alpha^4
        t3s = T3S*alpha^2        
        t1sp = D1S*alpha^3
        t2sp = D2S*alpha^5
        t3sp = D3S*alpha^3                        
        # Compute weights with hack
        d1 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 1.0, 0.0,
                                        t1s, t1sp, t2s, t2sp, t3s, t3sp)
        d2 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 0.0, 1.0,
                                        t1s, t1sp, t2s, t2sp, t3s, t3sp)
        dS = abs(dtau[i])
        DUDX[1, i] = (d1[1] + 1im*d1[3])*dS 
        DUDX[2, i] = (d2[1] + 1im*d2[3])*dS 
        DUDY[1, i] = (d1[2] + 1im*d1[4])*dS 
        DUDY[2, i] = (d2[2] + 1im*d2[4])*dS                 
    end
    DUDX = transpose(vec(DUDX))
    DUDY = transpose(vec(DUDY))
    return DUDX, DUDY
end


# Compute weights for gradient special quadrature
function gradient_kernel_weights(tau, dtau, nu, kappa, tau1, tau2, nu1, nu2,
                                 omega, omega1, omega2, 
                                 z, alpha, wL, wC, wQ, wT, DD, partint::Bool)
    tauminusz = tau - z
    # Complex panel differentiation matrix
    DDtau = diagm(1 ./ (DD*tau))*DD    
    domegadtau = DDtau*omega
    domegabardtau = DDtau*conj(omega)    
    DUDX, DUDY = _smoothlog_term(tauminusz, nu, dtau, wL, wC, alpha)
    Q = _stresslet_term(omega, domegadtau, domegabardtau,
                        tau, z, tauminusz, nu, kappa, dtau, wC, wQ, wT,
                        tau1, tau2, omega1, omega2, nu1, nu2, partint)
    C = _cauchy_term(omega, dtau, z, tauminusz, wC, wQ, domegadtau, nu, alpha, tau1, tau2, omega1, omega2, partint)    
    @. DUDX += Q[1]+C[1]
    @. DUDY += Q[2]+C[2]        
    return DUDX, DUDY
end

function _smoothlog_term(tauminusz, nu, dtau, wL, wC, alpha)
    N = length(tauminusz)
    DUDX = complex(zeros(2, N))
    DUDY = complex(zeros(2, N))    
    for i=1:N
        # This term is done in primary vars
        r1 = real(tauminusz[i])
        r2 = imag(tauminusz[i])
        n1 = real(nu[i])
        n2 = imag(nu[i])        
        # Kernel functions
        rnorm = sqrt(r1^2 + r2^2)
        T1S, T1L, T2S, T2L, T3S, T3L = Ti_split(alpha*rnorm)
        D1S, D1L, D2S, D2L, D3S, D3L = Ti_split_prime(alpha*rnorm)        
        t1s = (T1S + T1L*log(alpha))*alpha^2
        t2s = (T2S + T2L*log(alpha))*alpha^4
        t3s = (T3S + T3L*log(alpha))*alpha^2        
        t1l = T1L*alpha^2
        t2l = T2L*alpha^4
        t3l = T3L*alpha^2
        t1sp = (D1S + D1L*log(alpha))*alpha^3
        t2sp = (D2S + D2L*log(alpha))*alpha^5
        t3sp = (D3S + D3L*log(alpha))*alpha^3                        
        t1lp = D1L*alpha^3
        t2lp = D2L*alpha^5
        t3lp = D3L*alpha^3                
        ##### LOG PART
        # Unit vector responses
        e1 = gen_mod_stresslet(r1, r2, n1, n2, 1.0, 0.0, t1l, t2l, t3l)
        e2 = gen_mod_stresslet(r1, r2, n1, n2, 0.0, 1.0, t1l, t2l, t3l)
        d1 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 1.0, 0.0,
                                        t1l, t1lp, t2l, t2lp, t3l, t3lp)
        d2 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 0.0, 1.0,
                                        t1l, t1lp, t2l, t2lp, t3l, t3lp)
        tmp = conj( 1im./nu[i].*wC[i] ) # Integration weight for |dtau|/(tau-z)
        DUDX[1, i] = -(e1[1] + 1im*e1[2]) * real(tmp)
        DUDX[2, i] = -(e2[1] + 1im*e2[2]) * real(tmp)
        DUDY[1, i] = -(e1[1] + 1im*e1[2]) * imag(tmp)
        DUDY[2, i] = -(e2[1] + 1im*e2[2]) * imag(tmp)
        DUDX[1, i] += (d1[1] + 1im*d1[3])*wL[i]
        DUDX[2, i] += (d2[1] + 1im*d2[3])*wL[i]
        DUDY[1, i] += (d1[2] + 1im*d1[4])*wL[i]
        DUDY[2, i] += (d2[2] + 1im*d2[4])*wL[i]
        ##### SMOOTH PART
        d1 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 1.0, 0.0,
                                        t1s, t1sp, t2s, t2sp, t3s, t3sp)
        d2 = gen_mod_stresslet_gradient(r1, r2, n1, n2, 0.0, 1.0,
                                        t1s, t1sp, t2s, t2sp, t3s, t3sp)
        dS = abs(dtau[i])
        DUDX[1, i] += (d1[1] + 1im*d1[3])*dS 
        DUDX[2, i] += (d2[1] + 1im*d2[3])*dS 
        DUDY[1, i] += (d1[2] + 1im*d1[4])*dS 
        DUDY[2, i] += (d2[2] + 1im*d2[4])*dS                 
    end
    DUDX = transpose(vec(DUDX))
    DUDY = transpose(vec(DUDY))
    return DUDX, DUDY
end

function _stresslet_term(omega, domegadtau, domegabardtau,
                         tau, z, tauminusz, nu, kappa, dtau, wC, wQ, wT,
                         tau1, tau2, omega1, omega2, nu1, nu2, partint)    
    ##### STRESSLET PART
    if partint
        # Partial integration of as much as possible
        # dfdz = (        
        #     -wC.'*domegadtau
        #     - conj(wC).'*conj(domegadtau)
        # )
        # dfdz_corr = 2*real(omega2/(tau2-z) - omega1/(tau1-z))       
        # tmp = 2 * kappa ./ nu.^3        
        # dfdzbar = conj(
        #     + (wC.'*domegabardtau)
        #     + (wC.*conj(nu.^2)).'*domegadtau
        #     + (wC.*tmp).' * omega
        #     + (wQ.*conj(tauminusz)).'*domegadtau
        #     - (wQ.*conj(nu.^2)).'*omega
        # )
        # fG(tau, omega, nu) = -conj( (conj(omega) + omega.*conj(nu).^2)./(tau-z)
        #                             +  omega.*conj(tau-z)./(tau-z).^2 )    
        # dfdzbar_corr = fG(tau2,omega2,nu2)-fG(tau1,omega1,nu1)

        # Only partial integration of worst term
        dfdz = (
            -(transpose(wQ)*omega + transpose(conj(wQ))*conj(omega))
        )
        dfdz_corr = 0.0        
        dfdzbar = (
            + conj(wQ).'*omega
            + conj(conj(nu).^2.*wQ).'*conj(omega)
            + conj((wQ.*conj(tauminusz)).'*domegadtau
                   - (wQ.*conj(nu.^2)).'*omega)
        )
        fG(tau, omega, nu) = -conj( omega.*conj(tau-z)./(tau-z).^2 )
        dfdzbar_corr = fG(tau2,omega2,nu2)-fG(tau1,omega1,nu1)        
    else
        # Direct derivatives
        dfdz = (
            -(wQ.'*omega + conj(wQ).'*conj(omega))
        )
        dfdz_corr = 0.0        
        dfdzbar = (
            + conj(wQ).'*omega
            + conj(conj(nu).^2 .* wQ).'*conj(omega)
            + 2*conj(conj(tauminusz).*wT).'*conj(omega)
        )
        dfdzbar_corr = 0.0
    end
    
    dudx = (dfdz + dfdzbar)/(-4*1im*pi)
    dudy = 1im*(dfdz - dfdzbar)/(-4*1im*pi)

    dudx_corr = (dfdz_corr + dfdzbar_corr)/(-4*1im*pi)
    dudy_corr = 1im*(dfdz_corr - dfdzbar_corr)/(-4*1im*pi)    
    return (dudx + dudx_corr,
            dudy + dudy_corr)
end

function _cauchy_term(omega, dtau, z, tauminusz, wC, wQ, domegadtau, nu, alpha, tau1, tau2, omega1, omega2, partint)
    ##### CAUCHY PART
    # Direct derivatives
    dfdz_dir = (
        + dtau.'*conj(omega)
        - conj(dtau).'*omega
        - 2*(tauminusz.*conj(wC)).'*conj(omega)
    )
    
    dfdzbar_dir = (
        + dtau.'*omega
        + conj(conj(tauminusz.^2).*wQ).'*conj(omega)
    )

    partint = false # Partial integration of this terms has no benefit
    if partint    
        # Partial integration:
        fG(tau,omega) = -conj( omega.*conj(tau-z).^2 ./ (tau-z) )
        dfdzbar_pi = (dtau.'*omega
                      + conj(
                          (conj(tauminusz.^2).*wC).'*domegadtau
                          +
                          (2 .* conj(tauminusz).*(-conj(nu.^2)).*wC).'*omega
                      )
                      ) + fG(tau2,omega2)-fG(tau1,omega1)
        # Put it together
        dfdzbar = dfdzbar_pi
    else
        dfdzbar = dfdzbar_dir
    end
    dfdz = dfdz_dir    
    dudx = (dfdz + dfdzbar)/(4*1im)*(alpha^2/(8*pi))
    dudy = 1im*(dfdz - dfdzbar)/(4*1im)*(alpha^2/(8*pi))
    return dudx, dudy
end


##### Helpers for structure that appears in S and L terms:

function gen_mod_stresslet(r1, r2, n1, n2, q1, q2, T1, T2, T3)
    rdotn = r1*n1 + r2*n2
    rdotq = r1*q1 + r2*q2
    ndotq = n1*q1 + n2*q2        
    u1 = (
        T1*(rdotq*n1 + ndotq*r1 + rdotn*q1) +
        T2*(rdotq*rdotn*r1) +
        T3*(ndotq*r1)
    )
    u2 = (
        T1*(rdotq*n2 + ndotq*r2 + rdotn*q2) +
        T2*(rdotq*rdotn*r2) +
        T3*(ndotq*r2)
    )
    return u1, u2    
end

function gen_mod_stresslet_gradient(r1, r2, n1, n2, q1, q2, T1, T1p, T2, T2p, T3, T3p)
    rnorm = sqrt(r1^2 + r2^2)
    rdotn = r1*n1 + r2*n2
    rdotq = r1*q1 + r2*q2
    ndotq = n1*q1 + n2*q2       
    du1dx = -(
        T1p/rnorm*r1*(rdotq*n1 + ndotq*r1 + rdotn*q1) +
        T1*(q1*n1 + ndotq + n1*q1) +
        T2p/rnorm*r1*(rdotq*rdotn*r1) +
        T2*(q1*rdotn*r1 + rdotq*n1*r1 + rdotq*rdotn) +
        T3p/rnorm*r1*ndotq*r1 +
        T3*(ndotq)
    )
    du1dy = -(
        T1p/rnorm*r2*(rdotq*n1 + ndotq*r1 + rdotn*q1) +
        T1*(q2*n1 + n2*q1) +
        T2p/rnorm*r2*(rdotq*rdotn*r1) +        
        T2*(q2*rdotn*r1 + rdotq*n2*r1) + 
        T3p/rnorm*r2*ndotq*r1
    )
    du2dy = -(
        T1p/rnorm*r2*(rdotq*n2 + ndotq*r2 + rdotn*q2) +
        T1*(q2*n2 + ndotq + n2*q2) +
        T2p/rnorm*r2*(rdotq*rdotn*r2) +
        T2*(q2*rdotn*r2 + rdotq*n2*r2 + rdotq*rdotn) +
        T3p/rnorm*r2*ndotq*r2 +
        T3*(ndotq)
    )
    du2dx = -(
        T1p/rnorm*r1*(rdotq*n2 + ndotq*r2 + rdotn*q2) +
        T1*(q1*n2 + n1*q2) +
        T2p/rnorm*r1*(rdotq*rdotn*r2) +        
        T2*(q1*rdotn*r2 + rdotq*n1*r2) + 
        T3p/rnorm*r1*ndotq*r2
    )
    return du1dx, du1dy, du2dx, du2dy
end
