push!(LOAD_PATH, string(pwd(),"/src"))

import ModifiedStokesSolver
using Base.Test

mss = ModifiedStokesSolver

@testset "Split" begin

    srand(1)
    r = 1 .- 2.*(rand(), rand())
    f = 1 .- 2.*(rand(), rand())
    n = 1 .- 2.*(rand(), rand())
    rnorm = sqrt(r[1]^2+r[2]^2)

    cutoff = mss.POWERSERIES_ZMAX

    @testset "T1-T3" begin
        # Test split of stresslet kernel functions
        z = 0.987 # This value should work well for both
        (T1,T2,T3) = mss.Ti_direct(z)
        (T1S, T1L, T2S, T2L, T3S, T3L) =
            mss.Ti_split(z)
        E1 = (T1S + T1L*log(z) - T1)/T1
        E2 = (T2S + T2L*log(z) + 1/(8*pi*z^2) - 1/(pi*z^4)
              - T2)/T2        
        E3 = (T3S + T3L*log(z) - T3)/T3
        
        @test abs(E1) < 1e-14
        @test abs(E2) < 1e-14        
        @test abs(E3) < 1e-14    

        if cutoff>0
            # Test with dispatcher for large and small arguments
            # 3.5 is the largest argument that works with explicit split
            for z=[cutoff/2 cutoff*2 1e-5 3.5]
                @show z
                (T1,T2,T3) = mss.Ti(z)
                (T1S, T1L, T2S, T2L, T3S, T3L) =
                    mss.Ti_split(z)
                E1 = (T1S + T1L*log(z) - T1)/T1
                E2 = (T2S + T2L*log(z) + 1/(8*pi*z^2) - 1/(pi*z^4)
                      - T2)/T2        
                E3 = (T3S + T3L*log(z) - T3)/T3
                
                @test abs(E1) < 2e-14
                @test abs(E2) < 1e-14        
                @test abs(E3) < 1e-14

            end
        end
    end

    @testset "Derivatives" begin
        @testset "Split" begin
            z = 0.9999*ModifiedStokesSolver.POWERSERIES_ZMAX_PRIME
            err = abs.( (mss.Ti_split_pow(z, true) .- mss.Ti_split_direct(z, true)) ./ mss.Ti_split_direct(z, true) )
            @test maximum(err) < 1e-13
        end
        @testset "Recombine" begin
            z = 0.99*ModifiedStokesSolver.POWERSERIES_ZMAX_PRIME
            err = abs.( (mss.Ti_prime(z) .- mss.Ti_direct_prime(z)) ./ mss.Ti_direct_prime(z) )
            @test maximum(err) < 3e-14
        end
    end

    @testset "Stresslet" begin
        # Test split stresslet
        r = rand(2)
        f = rand(2)
        n = rand(2)
        alpha = rand()
        u1 = mss.stresslet(r, f, n, alpha)
        (uS, uL, uC, uQ) = mss.stresslet_split(r, f, n, alpha)
        rnorm = norm(r)
        rdotf = dot(r,f)
        rdotn = dot(r,n)    
        u2 = uS + uL*log(rnorm) + uC*rdotn/rnorm^2 +
            uQ*rdotf*r*rdotn/rnorm^4
        E = norm(u1-u2,Inf)/norm(u2,Inf)
        @test E < 1e-14
    end
end
