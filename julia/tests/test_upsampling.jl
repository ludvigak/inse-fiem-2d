push!(LOAD_PATH, string(pwd(),"/src"))

using ModifiedStokesSolver
using FastGaussQuadrature
using Base.Test

panelorder1 = 32
panelorder2 = 160

glpoints1, glweights1 = gausslegendre(panelorder1)
glpoints2, glweights2 = gausslegendre(panelorder2)

@testset "Upsampling" begin
    U = ModifiedStokesSolver.upsampling_matrix(panelorder1, panelorder2)
    f(x) = @. 1/(x^2+2)
    fu = U * f(glpoints1)
    err = fu - f(glpoints2)
    relerr_inf = norm(err, Inf) / norm(fu, Inf)
    @test relerr_inf < 2e-14
end

