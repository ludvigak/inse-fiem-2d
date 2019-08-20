push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using SpecialFunctions
using ModifiedStokesSolver
using Base.Test

@testset "NewPowerSeriesTests" begin

    z = 0.99
    #(T1Sp, T1Lp, T2Sp, T2Lp, T3Sp, T3Lp) =
    Tp = ModifiedStokesSolver.Ti_split_pow(z)    
    #(T1Sd, T1Ld, T2Sd, T2Ld, T3Sd, T3Ld) =
    Td = ModifiedStokesSolver.Ti_split_direct(z)    
    Ep = Tp .- Td
    @test maximum(abs.(Ep./Td)) < 9e-14

    z = 1.8
    #(T1Sp, T1Lp, T2Sp, T2Lp, T3Sp, T3Lp) =
    Tp = ModifiedStokesSolver.Ti_split_pow(z,true)    
    #(T1Sd, T1Ld, T2Sd, T2Ld, T3Sd, T3Ld) =
    Td = ModifiedStokesSolver.Ti_split_direct(z,true)    
    Ep = Tp .- Td
    @test maximum(abs.(Ep./Td)) < 9e-14

    
end



# @testset "PowerSeries" begin

# # Series evaluator
# function evalseries(fn, z, nstart)
#     y = 0.0
#     dy = 1.0
#     n = nstart
#     while abs(dy)>=eps(y)
#         dy = fn(z, n)*z^(2*n)
#         y += dy
#         n += 1
#     end
#     return y
# end
# # Definitions for modified Bessel functions
# i0(z,n) = 1/(4^n*factorial(n)^2)
# i1(z, n) = 1/(2*4^n*factorial(n)*factorial(n+1))
# k0(z,n) = ((digamma(n+1) + log(2) - log(z)) / (4^n * factorial(n)^2))
# k1(z, n) = (2*log(z) - 2*log(2) - digamma(n+1) - digamma(n+2)) /
#     (4^(n+1)*factorial(n)*factorial(n+1))
# # Splits
# k0S(z,n) = ((digamma(n+1) + log(2)) / (4^n * factorial(n)^2))
# k0L(z,n) = (-1) / (4^n * factorial(n)^2)
# k1S(z, n) = ( - 2*log(2) - digamma(n+1) - digamma(n+2)) /
#     (4^(n+1)*factorial(n)*factorial(n+1))
# k1L(z, n) = (2) / (4^(n+1)*factorial(n)*factorial(n+1))

# # Direct power series for split K0, K1
# # Returns: (K0S, K0L, K1S, K1L, K1C)
# function besselpowsplit(z)
#     # Outputs
#     K0S = 0
#     K0L = 0
#     K1S = 0
#     K1L = 0
#     K1C = 1
#     # Recursion variables
#     n = 0
#     z2 = z*z
#     z2n = 1; # z^(2*n)
#     digamma_np1 = -eulergamma # digamma(n+1)
#     digamma_np2 = -eulergamma + 1
#     ifact_n = 1 # 1/factorial(n)
#     ifact_np1 = 1
#     i4n = 1 # 1/4^n
#     i4np1 = 1/4
#     log2 = log(2)
#     # Loop
#     maxit = 20
#     while n<maxit
#         k0n_denom = i4n*ifact_n*ifact_n
#         k1n_denom = i4np1*ifact_n*ifact_np1
#         k0Ln = (-1)*k0n_denom
#         k0Sn = (digamma_np1 + log2)*k0n_denom
#         k1Ln = (2)*k1n_denom
#         k1Sn = (-2*log2-digamma_np1-digamma_np2)*k1n_denom
#         K0S += k0Sn*z2n
#         K0L += k0Ln*z2n
#         K1S += k1Sn*z2n*z
#         K1L += k1Ln*z2n*z
#         # Update recursions
#         n += 1
#         z2n *= z2;
#         inp1 = 1/(n+1)
#         (digamma_np1, digamma_np2) = (digamma_np2, digamma_np2+inp1)
#         (ifact_n, ifact_np1) = (ifact_np1, ifact_np1*inp1)
#         (i4n, i4np1) = (i4np1, i4np1*0.25)
#     end
#     return (K0S, K0L, K1S, K1L, K1C)
# end

# # Tests
# z = 1/10

# ref = besseli(0, z)
# I0test = evalseries(i0, z, 0)
# @test abs((ref-I0test)/ref) < 1e-14
# ref = besseli(1, z)
# I1test = z*evalseries(i1, z, 0)
# @test abs((ref-I1test)/ref) < 1e-14
# ref = besselk(0, z)
# K0test = evalseries(k0, z, 0)
# K0split = evalseries(k0S, z, 0) + log(z)*evalseries(k0L, z, 0)
# @test abs((ref-K0test)/ref) < 1e-14
# @test abs((ref-K0split)/ref) < 1e-14
# ref = besselk(1, z)
# K1test = 1/z + z*evalseries(k1, z, 0)
# K1split = 1/z + z*evalseries(k1S, z, 0) + z*log(z)*evalseries(k1L, z, 0)
# @test abs((ref-K1test)/ref) < 1e-14
# @test abs((ref-K1split)/ref) < 1e-14
# # Test direct split power series
# ref = (besselk(0, z), besselk(1, z))
# (K0S, K0L, K1S, K1L, K1C) = besselpowsplit(z)
# split = (log(z)*K0L+K0S, K1S+log(z)*K1L+K1C/z)
# @test maximum(@. abs((ref-split)/ref) ) < 1e-14


# # Evaluate T, high precision reference from Mathematica
# for n=1:2
#     if n==1
#         z = 0.001
#         T1exact = 0.00028941085628129279453492498775927
#         T2exact = -318.30984639510755467178629422316
#         T3exact = -0.00059871605412101857968611091116825
#     end
#     if n==2
#         z = 0.99
#         T1exact = 0.023658849757364503594228422585464
#         T2exact = -0.28952030153253490875314823065265
#         T3exact = -0.063320117528153811164290571635020
#     end
#     # z = 1
#     # (T1exact, T2exact, T3exact) =
#     #     (0.023620976507380996090769874829016,
#     #      -0.28607692796680641340443170447384,
#     #      -0.063358432123254121248207660793628)

#     # z = 2
#     # (T1exact, T2exact, T3exact) =
#     #     (0.016930005787967333888517613314308,
#     #      -0.11224071607388210664375507421808,
#     #      -0.057317125084941282339599571205834)


#     println("z = ", z)
#     println()

#     K0 = besselk(0, z)
#     K1 = besselk(1, z)
#     K0S = K0 + besseli(0, z)*log(z)
#     K1S = K1 - besseli(1, z)*log(z) - 1/z    
#     I0 = besseli(0,z)
#     I1 = besseli(1,z)

#     ( T1Spow, T1Lpow,
#       T2Spow, T2Lpow, T2Cpow,
#       T3Spow, T3Lpow) = ModifiedStokesSolver.Ti_split_pow(z)

#     ## T1
#     T1 = -( z^2*(2*K0 + z*K1) + 4*(z*K1 - 1) ) / (2*pi*z*z^2)
#     t1(z, n) = 2*k0(z, n) + 4*k1(z, n) + k1(z, n-1)
#     t1S(z, n) = 2*k0S(z, n) + 4*k1S(z, n) + k1S(z, n-1)
#     t1L(z, n) = 2*k0L(z, n) + 4*k1L(z, n) + k1L(z, n-1)
#     T1series = -evalseries(t1,z,1)/(2*pi*z)
#     T1Sseries = -evalseries(t1S,z,1)/(2*pi*z)
#     T1Lseries = -evalseries(t1L,z,1)/(2*pi*z)
#     T1splitseries = (T1Sseries + T1Lseries*log(z))
#     T1pow = (T1Spow + T1Lpow*log(z))

#     #quit()

#     T1S = -(2*z*K0S + (z^2+4)*K1S + z) / (2*pi*z^2)
#     T1L = (2*z*I0 - (z^2+4)*I1) / (2*pi*z^2)
#     T1split = (T1S + T1L*log(z))

#     println("== T1")
#     println("exact  = ", T1exact)
#     println("series = ", T1series, "  relerr = ", (T1series-T1exact)/T1exact)
#     T1spliterr = (T1splitseries-T1exact)/T1exact
#     println("spl+ser= ", T1splitseries, "  relerr = ", T1spliterr)
#     #@test abs(T1spliterr)<1e-14
#     println("regular= ", T1, "  relerr = ", (T1-T1exact)/T1exact)
#     println("reg+spl= ", T1split, "  relerr = ", (T1split-T1exact)/T1exact)
#     #@test abs( (T1series-T1exact)/T1exact )<1e-14
#     T1powerr = (T1pow-T1exact)/T1exact
#     println("final  = ", T1pow, "  relerr = ", T1powerr)
#     #@test abs(T1powerr)<1e-14
#     println()

#     ## T2
#     T2 = (4*z^2*K0 + (z^2+8)*z*K1 - 8) / (pi*z^3)
#     t2(z, n) = 4*k0(z, n) + 8*k1(z, n) + k1(z, n-1)
#     T2series = -1/(pi*z) + evalseries(t2,z,1)/(pi*z)
#     println("== T2")
#     println("exact  = ", T2exact)
#     println("series = ", T2series, "  relerr = ", (T2series-T2exact)/T2exact)
#     println("regular= ", T2, "  relerr = ", (T2-T2exact)/T2exact)
#     @test abs( (T2series-T2exact)/T2exact )<1e-14
#     T2pow = (T2Spow + T2Lpow*log(z) + T2Cpow/z)
#     T2powerr = (T2pow-T2exact)/T2exact
#     println("final  = ", T2pow, "  relerr = ", T2powerr)
#     @test abs(T2powerr)<1e-14
#     println()

#     ## T3
#     T3 = (z*K1 - 1) / (2*pi*z)    
#     t3(z, n) = k1(z, n)
#     T3series = z*evalseries(t3,z,0)/(2*pi)
#     println("== T3")
#     println("exact  = ", T3exact)
#     println("series = ", T3series, "  relerr = ", (T3series-T3exact)/T3exact)
#     println("regular= ", T3, "  relerr = ", (T3-T3exact)/T3exact)
#     @test abs( (T3series-T3exact)/T3exact )<1e-14
#     T3pow = (T3Spow + T3Lpow*log(z))
#     T3powerr = (T3pow-T3exact)/T3exact
#     println("final  = ", T3pow, "  relerr = ", T3powerr)
#     @test abs(T3powerr)<1e-14

#     println()
# end


#end # End testset
