push!(LOAD_PATH, string(pwd(),"/src"))
using ModifiedStokesSolver
using Base.Test

@testset "FMM" begin
    ## define source and target points
    ns = 100
    nt = 102

    srand(1)
    src = 1 - 2*rand(2,ns)
    fvec = 1 - 2*rand(2,ns)
    nvec = 1 - 2*rand(2,ns)
    targ = 1 - 2*rand(2,nt)

    # Compute using the kernels from ModifiedStokesSolver
    function mss_stokeslet_direct(src, targ, fvec, alpha)
        ns = size(src, 2)
        nt = size(targ, 2)
        u = zeros(2, nt)
        for i=1:ns
            for j=1:nt
                f = fvec[:, i]
                r = src[:,i]-targ[:,j]
                u[:,j] += stokeslet(r, f, alpha)
            end
        end
        return u
    end
    function mss_stresslet_direct(src, targ, fvec, nvec, alpha)
        ns = size(src, 2)
        nt = size(targ, 2)
        u = zeros(2, nt)
        for i=1:ns
            for j=1:nt
                r = src[:,i]-targ[:,j]
                u[:,j] += stresslet(r, fvec[:,i], nvec[:,i], alpha)
            end
        end
        return u
    end
    function mss_stresslet_direct_self(src, fvec, nvec, alpha)
        ns = size(src, 2)
        u = zeros(2, ns)
        for i=1:ns
            for j=1:ns
                i==j && continue
                r = src[:,i]-src[:,j]
                u[:,j] += stresslet(r, fvec[:,i], nvec[:,i], alpha)
            end
        end
        return u
    end

    
    ## Test with large alpha
    alpha = 20.0
    
    ## Stokeslet
    u1 = mss_stokeslet_direct(src, targ, fvec, alpha)
    u2 = fmm_stokeslet_direct(src, targ, fvec, alpha)
    u3 = fmm_stokeslet_targ(src, targ, fvec, alpha)    
    println(" * Stokeslet, alpha=", alpha)
    Gerrdirect = norm(u1[:]-u2[:]) / norm(u1[:], Inf)
    Gerrfmm = norm(u1[:]-u3[:]) / norm(u1[:], Inf)
    @test Gerrdirect < 5e-14
    @test Gerrfmm < 5e-14

    ## Stresslet
    ## Test also with small alpha
    ## (Stokeslet has no special eval for small args)
    for alpha=[alpha, 1e-15]
        println(" * Stresslet, alpha=", alpha)        
        # Stresslet targ
        d1 = mss_stresslet_direct(src, targ, fvec, nvec, alpha)
        d2 = fmm_stresslet_direct(src, targ, fvec, nvec, alpha)
        d3 = fmm_stresslet_targ(src, targ, fvec, nvec, alpha)
        # Stresslet self
        ds1 = mss_stresslet_direct_self(src, fvec, nvec, alpha)
        ds2 = fmm_stresslet_direct(src, src, fvec, nvec, alpha; self=true)    
        ds3 = fmm_stresslet_self(src, fvec, nvec, alpha)
        # Test
        Terrdirect = norm(d1[:]-d2[:]) / norm(d1[:], Inf)
        Terrdirect_self = norm(ds1[:]-ds2[:]) / norm(d1[:], Inf)    
        Terrfmm = norm(d1[:]-d3[:]) / norm(d1[:], Inf)
        Terrfmm_self = norm(ds1[:]-ds3[:]) / norm(ds1[:], Inf)    
        @test Terrdirect < 5e-14
        @test Terrdirect_self < 5e-14    
        @test Terrfmm < 5e-14
        @test Terrfmm_self < 5e-14    
    end
end
