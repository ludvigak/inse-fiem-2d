using FastGaussQuadrature
import GSL

"""
    U = upsampling_matrix(n, m)

um = U*un upsamples a function u(x) from n to m GL nodes on [-1,1] \\
Can be used to create a downsampling matrix as well
"""
function upsampling_matrix(n, m)
    L = legendre_matrix(n)
    B = Array{Float64}(m, n)
    glpoints, glweights = gausslegendre(m)        
    for i=1:m
        P = GSL.sf_legendre_Pl_array(n-1, glpoints[i])
        for j=1:n
            B[i,j] = P[j]
        end
    end
    return B*L
end

"""
    L = legendre_matrix(order)

L computes Legendre expansion coefficients for a function
defined at the Gauss-Legendre quadrature nodes on [-1, 1] \\
``c = L*f  ->  c_l = (2l+1)/2 \\sum_n P_l(x_n) f_n w_n``
"""
function legendre_matrix(n)
    glpoints, glweights = gausslegendre(n)    
    L = Array{Float64}(n, n)
    for i=1:n
        P = GSL.sf_legendre_Pl_array(n-1, glpoints[i])
        for j=1:n
            l = j-1
            L[j, i] = P[j]*glweights[i]*(2*l+1)/2
        end
    end
    return L
end

""" 
    D = legendre_diff_matrix(order)

D approximates the derivative of a function defined at 
the Gauss-Legendre quadrature nodes on [-1, 1] \\
``D*f =~ f'``
"""
function legendre_diff_matrix(n)
    glpoints, glweights = gausslegendre(n)    
    dPt = Array{Float64}(n, n)
    for i=1:n
        _, dP = legendre_Pl_deriv_array(n-1, glpoints[i])
        for j=1:n
            dPt[i, j] = dP[j]
        end
    end
    L = legendre_matrix(n)
    return dPt*L
end


"""
P, D = legendre_Pl_deriv_array(n, x) \\
Recurrence relation used: [http://dlmf.nist.gov/18.9]

As fast as the GSL implementation, but also takes x outside [-1, 1].
    """
function legendre_Pl_deriv_array(n, x)
    P = Array{Complex{Float64}}(n+1)
    D = Array{Complex{Float64}}(n+1)
    P[1] = 1.0 # l=0
    D[1] = 0.0
    P[2] = x # l=1
    D[2] = 1.0
    for l=1:n-1
        # Compute l+1
        P[l+1+1] = ( (2*l+1)*x * P[l+1] - l*P[l-1+1] ) / (l+1)
        D[l+1+1] = ( (2*l+1)*(P[l+1] + x * D[l+1]) - l*D[l-1+1] ) / (l+1)
    end
    P, D
end

"""
    t, relres, iter = newton_legendre(c, zr, t0, maxiter, tol)

Find root t of Legendre expansion c such that
``\\sum_n c_n*P_n(t) = zr``
"""
function newton_legendre(c, zr, t0, maxiter, tol)
    t = t0
    iter = 0
    dt = Inf
    lmax = length(c)-1
    relres = tol
    for iter=1:maxiter
        P, D = legendre_Pl_deriv_array(lmax, t)
        f = sum(P.*c) - zr
        fp = sum(D.*c)
        dt = -f / fp
        t = t + dt
        relres = abs(dt/t)
        if relres < tol
            break
        end
    end
    return t, relres, iter
end
