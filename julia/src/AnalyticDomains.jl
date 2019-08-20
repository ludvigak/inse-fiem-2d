#__precompile__()
module AnalyticDomains

# Curve parametrization
struct CurveDesc
    tau
    dtau
    d2tau
    normal
    tangent
end

export starfish, CurveDesc

using SymPy

# Use any curve parametrized in t in [0, 2 pi), described as string
# Ex. desc = paramcurve("(1 + 0.3*cos(5*t))*exp(I*t)")
# !Note that imaginary unit is I
function paramcurve(str::String)
    desc = Sym(str)
    t = symbols("t")
    # Differentiate curve
    diff1 = diff(desc, t)
    diff2 = diff(diff1, t)
    tang = diff1 / abs(diff1)
    normal = 1im*tang    
    return CurveDesc(
        lambdify(desc),
        lambdify(diff1),
        lambdify(diff2),
        lambdify(normal),
        lambdify(tang)
    )
end

function starfish(;n_arms = 5, amplitude = 0.3, radius = 1.0, center=(0,0), rotation=0.0, interior=false)
    if interior
        s = -1.0
    else
        s = 1.0
    end
    R = exp(1im*rotation)
    desc(t) = R*radius*((1 + amplitude*cos(n_arms*t)).*exp(1im*s*t)) + (center[1] + 1im*center[2])
    diff1(t) = R*radius*(exp(t*1im*s)*(amplitude*cos(n_arms*t) + 1)*1im*s - amplitude*n_arms*exp(t*1im*s)*sin(n_arms*t))
    diff2(t) = R*radius*(- exp(t*1im*s)*(amplitude*cos(n_arms*t) + 1) - amplitude*n_arms*exp(t*1im*s)*sin(n_arms*t)*2*1im*s - amplitude*n_arms^2*exp(t*1im*s)*cos(n_arms*t))
    tang(t) = diff1(t) / abs(diff1(t))
    normal(t) = 1im*tang(t)
    return CurveDesc(
        desc,
        diff1,
        diff2,
        normal,
        tang
    )
end

function rectangle(;width = 1, height = 1, order = 10, interior=false)
    @assert iseven(order)
    a = width
    b = height
    p = order;
    desc(t) = (a*cos(t) + b*sin(t)*1im)/(cos(t)^p + sin(t)^p)^(1/p)
    _diff1(t) = -(a*cos(t)*sin(t)^p - b*cos(t)^p*sin(t)*1im)/(cos(t)*sin(t)*(cos(t)^p + sin(t)^p)^((p + 1)/p))    
    _diff2(t) = (-(a*cos(t)^(p + 1)*sin(t)^p + b*cos(t)^(2*p + 2)*sin(t)*2im - 2*a*cos(t)^3*sin(t)^(2*p) - 2*a*cos(t)^(p + 3)*sin(t)^p - b*cos(t)^p*sin(t)^(p + 1)*1im - b*cos(t)^(2*p)*sin(t)*2im + b*cos(t)^(p + 2)*sin(t)^(p + 1)*2im + a*p*cos(t)^(p + 1)*sin(t)^p + b*p*cos(t)^p*sin(t)^(p + 1)*1im)/(cos(t)^2*sin(t)^2*(cos(t)^p + sin(t)^p)^((2*p + 1)/p)))

    # Fudge because derivative expression is NaN at t=0
    diff1(t) = _diff1(t + eps(t))
    diff2(t) = _diff2(t + eps(t))

    
    tang(t) = diff1(t) / abs(diff1(t))
    normal(t) = 1im*tang(t)
    
    return CurveDesc(
        desc,
        diff1,
        diff2,
        normal,
        tang
    )
end

end
