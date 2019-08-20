push!(LOAD_PATH, string(pwd(),"/src"))

using Compat.LinearAlgebra
using Compat.Test
using Compat.Random

using ModifiedStokesSolver

@testset "Kernels" begin
    srand(1)
    r = rand(2)
    f = rand(2)
    n = rand(2)
    alpha = rand()

    # Stokeslet functions are same
    u1 = stokeslet(r, f, alpha)
    K11, K12, K21, K22 = ModifiedStokesSolver.stokeslet_kernel_noarrays(r[1], r[2], alpha)
    K = [K11 K12; K21 K22]
    u2 = K*f
    @test norm(u1-u2, Inf) / norm(u1, Inf) < 1e-15   
    
    # Stresslet functions are same
    u1 = stresslet(r, f, n, alpha)
    u2 = dblkernel(r, n, alpha) * f
    @test norm(u1-u2, Inf) / norm(u1, Inf) < 1e-15

    # Test log split
    (uS, uL, uC, uQ) = stresslet_split(r, f, n, alpha)
    KL = dblkernel_log(r, n, alpha)
    uL2 = KL * f    
    @test norm(uL-uL2, Inf) / norm(uL, Inf) < 1e-15

    @testset "Ti derivatives" begin
        # Finite difference test
        h = 1e-8
        z = 2.0
        f = ModifiedStokesSolver.Ti_direct
        dFD = (f(z+h).-f(z-h))./(2*h)
        d = ModifiedStokesSolver.Ti_direct_prime(z)
        @test maximum(abs.((dFD .- d)./d)) < 1e-7
    end
    
    
end

    
