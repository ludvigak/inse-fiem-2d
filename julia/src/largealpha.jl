# Routines used when dealing with quadrature for large alpha

LARGEALPHA_LIMIT_16 = 4.5
LARGEALPHA_LIMIT_32 = 4.5

function set_largealpha_limit_16(x)
    global LARGEALPHA_LIMIT_16 = x
end

function set_largealpha_limit_32(x)
    global LARGEALPHA_LIMIT_32 = x
end

"""
    subdivide_interval(hmin, target)

Create a fixed level-restricted subdivision of [-1,1] starting from a panel of
length hmin centered around real(target).
"""
function subdivide_interval(hmin, target)
    intervals = [real(target)-hmin/2, real(target)+hmin/2]
    h = hmin
    while intervals[1] > -1 || intervals[end] < 1
        h = h*2 # This factor can be large, possibly up to 5
        push!(intervals, intervals[end]+h)
        unshift!(intervals, intervals[1]-h)
    end
    intervals = intervals[(intervals .> -1.0) .& (intervals .< 1.0)]
    push!(intervals, 1.0)
    unshift!(intervals, -1.0)
    return intervals
end    


"""
    subdivide_interval_with_bisection(dt_max, troot, Nsec)


Subdivide [-1,1] by placing a panel of size dt_max around real(troot), and then
subdividing neighboring panels until they satisfy quadrature requirements.

Generally produces fewer subdivisions than `subdivide_interval`.
"""
function subdivide_interval_with_bisection(dt_max, troot, Nsec)
    # First check: bail out if this algorithm is not needed (avoids a bit of computation)
    R = NEARLIM_RADIUS_16^(16/Nsec)  
    if !(bernstein_radius(troot) < R && 2.0 > dt_max)
        return [-1.0, 1.0]
    end
    # Second check: Vanilla bisection if target point outside [-1,1]
    if abs(real(troot)) > 1
        return recursive_bisection(-1.0, 1.0, R, dt_max, troot)
    end
    # Compute length of longest center panel s.t. specquad not needed
    d = abs(imag(troot))
    # bernstein_radius(d*1im) = d + sqrt(d^2+1)
    # So the closest center point not requiring specquad, d_nospec*1im, satisfies
    # d_nospec + sqrt(d_nospec^2+1) = R
    # Solving this:
    d_nospec = (R^2-1)/(2*R)
    # For this to be true, panel has to have length
    # (fudge factor to avoid roundoff kicking rho below R on subinterval)
    dt_nospec = 2*d/d_nospec * (1 - 1e-8)
    # Also figure out longest center panel that doesn't spill over
    distance_to_edge = 1-abs(real(troot))
    dt_nospill = 2*distance_to_edge
    # Pick center panel to be longest of dt_max and dt_nospec,
    # but not so long that it doesn't spill over
    dt_center = min(dt_nospill, max(dt_nospec, dt_max))
    # Create initial subdivision
    ta = real(troot)-dt_center/2
    tb = real(troot)+dt_center/2
    # Recursively bisect subintervals
    sub1 = recursive_bisection(-1.0, ta, R, dt_max, troot)
    sub2 = [ta, tb]
    sub3 = recursive_bisection(tb, 1.0, R, dt_max, troot)
    intervals = union(sub1, sub2, sub3)
    @assert intervals[1] == -1.0 && intervals[end] == 1.0
    return intervals
end

function recursive_bisection(ta, tb, R, dt_max, troot)
    if ta < tb
        dt = tb-ta
        troot_sec = (troot-ta)*2/(tb-ta) - 1.0                
        rho_sec = bernstein_radius(troot_sec)
        needs_specquad = rho_sec < R
        too_long = dt > dt_max
        if needs_specquad && too_long
            tc = ta + dt/2
            sub1 = recursive_bisection(ta, tc, R, dt_max, troot)
            sub2 = recursive_bisection(tc, tb, R, dt_max, troot)
            return union(sub1, sub2)
        end
    end
    return [ta, tb]
end
