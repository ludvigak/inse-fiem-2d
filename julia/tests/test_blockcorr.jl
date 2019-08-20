push!(LOAD_PATH, string(pwd(),"/src"))

using Base.Test
using ModifiedStokesSolver
import AnalyticDomains

@testset "BlockCorrections" begin
    # Params
    alpha = 10.0
    # Discretize
    numpanels = 20
    panelorder = 16
    curve = AnalyticDomains.starfish(n_arms = 3, amplitude=0.2)
    dcurve = discretize(curve, numpanels, panelorder)
    N = dcurve.numpoints

    # Assemble complete double layer matrix
    Dmat = doublelayer_matrix(dcurve, alpha)

    # Assemble block correction matrix
    BCmat = ModifiedStokesSolver.correction_block_matrix(dcurve, alpha)

    tovec(x) = reshape(x, :)
    fromvec(x) = reshape(x, 2, :)  
    matvec(x) = tovec(doublelayer_self(dcurve, fromvec(x), alpha, zerolimit=false)) + BCmat*x

    x = rand(2*N)
    v1 = Dmat*x
    v2 = matvec(x)
    relerr = norm(v1-v2, Inf) / norm(v1, Inf)

    @test relerr < 1e-14
end
