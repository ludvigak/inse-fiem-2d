push!(LOAD_PATH, string(pwd(),"/src"))

import AnalyticDomains
using ModifiedStokesSolver
using Base.Test

@testset "Discretization" begin
    numpanels = 100
    panelorder = 16
    curve = AnalyticDomains.starfish(n_arms = 5, amplitude=0.3);
    dcurve = discretize(curve, numpanels, panelorder);
    # Test equal arclength discretization
    h = sum(reshape(dcurve.dS, 16, :), 1)
    err = norm(h-h[1], Inf)
    @show err
    @test err < 1e-13
end
