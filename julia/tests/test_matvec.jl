# Test FMM matvec vs full matrix
push!(LOAD_PATH, string(pwd(),"/src"))
using Base.Test
using ModifiedStokesSolver
import AnalyticDomains

@testset "Matvec" begin
    # Params
    alpha = 10.0
    # Discretize
    numpanels = 20
    panelorder = 16
    curve = AnalyticDomains.starfish(n_arms = 3, amplitude=0.2)
    dcurve = discretize(curve, numpanels, panelorder)
    N = dcurve.numpoints
    # Compute
    x = rand(2*N)
    matvec = system_matvec(dcurve, alpha)
    yfmm = matvec(x)
    A = system_matrix(dcurve, alpha)
    ydir = A*x
    erel = norm(yfmm-ydir, Inf) / norm(ydir, Inf)
    @test erel < 5e-14
end
