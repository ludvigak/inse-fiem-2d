__precompile__()
module BarycentricLagrange
# Barycentric Lagrange Interpolation.
# Berrut, J.-P., & Trefethen, L. N. (2004).
# SIAM Review, 46(3), 501–517. doi:10.1137/S36144502417716

export bclag_interp_eval
export bclag_interp_matrix
export bclag_interp_weights

function bclag_interp_eval(x, f, xi, w)
    # fi = bclag_interp_eval(x, f, xi, w)
    # xi scalar, not in x
    tmp = w ./ (xi - x)
    fi = sum(f.*tmp)/sum(tmp)
    return fi
end
bclag_interp_eval(x, f, xi) = bclag_interp_eval(x, f, xi, bclag_interp_weights(x))

function bclag_interp_matrix(x, xi, w)
    n = length(x)
    N = length(xi)
    @assert length(w)==n

    B = zeros(N, n)
    denom = zeros(N)
    exact = zeros(Int64, N)
    for j=1:n
        for k=1:N
            xdiff = xi[k]-x[j]
            if xdiff != 0
                temp = w[j]/xdiff
                B[k,j] = temp
                denom[k] += temp
            else
                exact[k] = j
            end
        end
    end

    B ./= denom
    for jj=1:N
        if exact[jj] != 0
            B[jj,:] = 0.0
            B[jj + N*(exact[jj]-1)] = 1.0
        end
    end
    return B
end
bclag_interp_matrix(x, xi) = bclag_interp_matrix(x, xi, bclag_interp_weights(x))


function bclag_interp_weights(x)
    n = length(x);
    w = zeros(n);
    for j=1:n
        w[j] = 1/prod(x[j]-x[1:n.!=j]);
    end
    return w
end


end # module
