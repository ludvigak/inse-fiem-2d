# Test PUX by extending a function and computing gradient using Fourier method

push!(LOAD_PATH, string(pwd(),"/src"))
using ModifiedStokesSolver
using CurveDiscretization
import AnalyticDomains
using PUX
using Base.Test

@testset "PUX" begin
    srand(0)
    # Params
    alpha = 2.0
    numpanels = 10
    panelorder = 16
    Ngrid_int = 150
    pux_ep = 2.0
    pux_P = 40
    # Discretize
    curve = AnalyticDomains.starfish(n_arms=3, amplitude=0.1)
    dcurve = discretize(curve, numpanels, panelorder,
                        equal_arclength=true)
    ## PUX prep
    PUXStor, X, Y, Lgrid, Ngrid, interior = PUX.pux_precompute(curve, dcurve, Ngrid_int, pux_ep; P=pux_P)
    # Test function
    f(x, y) = x.^2 + y.^2 + x.*y
    fx(x, y) = 2*x + y
    fy(x, y) = 2*y + x   
    F = f.(X,Y)
    Fxref = fx.(X,Y)
    Fyref = fy.(X,Y)        
    # Extend
    FEXT = PUX.pux_eval(F, PUXStor)
    # Differentiate
    Fx, Fy = per_gradient(FEXT, Lgrid)
    Ex = abs.(Fx-Fxref)    
    Ey = abs.(Fy-Fyref)
    E = max.(Ex,Ey)
    logE = log10.(E)
    Emax = norm(vec(E[interior]), Inf)
    @test Emax < 1e-6
end
