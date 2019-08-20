# Machinery for complex plane mappings

"zr = rotate_and_scale(za, zb, z)"
function rotate_and_scale(za, zb, z)
    # computes
    # zr = M(z)
    # where
    # M(za) = -1
    # M(zb) = 1
    c_map = ((za+zb)/(za-zb), 2/(zb-za))
    zr = c_map[1] + z*c_map[2]
    return zr
end

function z_edges(grid, i)
    za = grid.edges[1,i] + 1im*grid.edges[2,i]
    zb = grid.edges[3,i] + 1im*grid.edges[4,i]
    return za, zb
end

function map_panels(grid::DiscreteCurve)
    L = legendre_matrix(grid.panelorder)
    ffit(z) = L*z
    coeffs = Array{Complex{Float64}}(grid.panelorder, grid.numpanels)
    for i=1:grid.numpanels
        za, zb = z_edges(grid, i)
        idx = (1:grid.panelorder) + grid.panelorder*(i-1)
        zpanel = grid.points[1,idx] + 1im*grid.points[2,idx]
        # First rescale [za,zb] to [-1,1], which helps Newton search later on
        zscaled = rotate_and_scale(za, zb, zpanel)
        c = ffit(zscaled)
        coeffs[:,i] = c
    end
    return coeffs
end

function invert_map(grid, coeffs, panel_idx, z; fall_back_to_initial::Bool=true)
    # Step 1: Rotate and scale
    za, zb = z_edges(grid, panel_idx)
    zr = rotate_and_scale(za, zb, z)
    # Step 2: Get the right expansion and solve with Newton
    c = coeffs[:, panel_idx]
    c = c[1 : min(length(c), 16)]     # Limit to 16 modes
    t0 = zr
    maxiter = 30
    converged = true
    t, relres, iter = newton_legendre(c, zr, t0, maxiter, 1e-13) # Try to hit 1e-13
    if iter == maxiter && relres > 1e-10 # ...but consider things converged if only 1e-10
        info("Newton did not converge, relres=", @sprintf("%.2e", relres), ". t=$t, t0=$t0")
        converged = false
        if fall_back_to_initial
            t = t0 # Fall back to initial guess
        end
    end                    
    tlocal = t
    # Step 3: Rescale t back to global parametrization
    ta = grid.t_edges[1, panel_idx]
    tb = grid.t_edges[2, panel_idx]
    tglobal = ta + (t+1)/2*(tb-ta)
    return tlocal, tglobal, converged
end
