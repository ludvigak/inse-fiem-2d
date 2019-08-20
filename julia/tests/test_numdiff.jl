push!(LOAD_PATH, string(pwd(),"/src"))

using ModifiedStokesSolver
using Base.Test

@testset "NumDiff" begin

    function numdiv(F, r, h)
        hx = [h, 0]
        hy = [0, h]
        Fx = (F(r.+hx) .- F(r.-hx)) ./ (2*h)
        Fy = (F(r.+hy) .- F(r.-hy)) ./ (2*h)
        return Fx[1] + Fy[2]
    end

    function numgrad(F, r, h)
        hx = [h, 0]
        hy = [0, h]
        Fx = (F(r.+hx) - F(r.-hx)) ./ (2*h)
        Fy = (F(r.+hy) - F(r.-hy)) ./ (2*h)
        return (Fx, Fy)
    end

    function numlap(F, r, h)
        hx = [h, 0]
        hy = [0, h]
        Fxx = (F(r.+hx) .-2.*F(r) .+ F(r.-hx)) ./ (h^2)
        Fyy = (F(r.+hy) .-2.*F(r) .+ F(r.-hy)) ./ (h^2)
        return Fxx .+ Fyy
    end


    h = 1e-4
    srand(1)
    r = 1 .- 2.*rand(2)
    f = 1 .- 2.*rand(2)
    n = 1 .- 2.*rand(2)
    a = rand()
    rnorm = sqrt(sum(r.^2))
    println("a|r| = ", a*rnorm)

    
    p = pressure(r, f, a)
    u = stokeslet(r, f, a)
    divu = numdiv(r -> stokeslet(r, f, a), r, h)
    Lp = numlap(r -> pressure(r, f, a), r, h)
    Lu = numlap(r -> stokeslet(r, f, a), r, h)
    gradp = numgrad(r -> pressure(r, f, a), r, h)
    pde = a^2.*u .- Lu .+ gradp

    println("u = ",u)
    println("div u = ",divu)
    println("L p = ", Lp)
    println("a^2 u - Lu + gradp = ", pde)

    # Test streslet with a|r| <1 and >1    
    T = stresslet(r, f, n, a)
    divT_zlt1 = numdiv(r -> stresslet(r, f, n, a), r, h)
    println("T = ", T)
    println("div T = ", divT_zlt1)

    a = 2/rnorm
    println("a|r| = ", a*rnorm)    
    T = stresslet(r, f, n, a)
    divT_zgt1 = numdiv(r -> stresslet(r, f, n, a), r, h)
    println("T = ", T)
    println("div T = ", divT_zgt1)

    
    @test abs(divu) < 1e-7
    @test abs(Lp) < 1e-7
    @test maximum(abs.(pde)) < 1e-7
    @test abs(divT_zlt1) < 1e-7
    @test abs(divT_zgt1) < 1e-7

end
