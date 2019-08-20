push!(LOAD_PATH, string(pwd(),"/../../code/src"))
using ModifiedStokesSolver
using Base.Test

include("ModStokesFMM2D.jl")

alpha = 20.0

## define source and target points

ns = 500
nt = 600

srand(1)
stokeslet_str = 1 - 2*rand(2,ns)
src = 1 - 2*rand(2,ns)
fvec = 1 - 2*rand(2,ns)
nvec = 1 - 2*rand(2,ns)
targ = 1 - 2*rand(2,nt)

# Compute using the kernels from ModifiedStokesSolver
function mss_stokeslet_direct(src, targ, stokeslet_str, alpha)
    ns = size(src, 2)
    nt = size(targ, 2)
    u = zeros(2, nt)
    for i=1:ns
        for j=1:nt
            f = stokeslet_str[:, i]
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

## Stokeslet
u1 = mss_stokeslet_direct(src, targ, stokeslet_str, alpha)
u2 = fmm_stokeslet_direct(src, targ, stokeslet_str, alpha)
u3 = fmm_stokeslet_targ(src, targ, stokeslet_str, alpha)

## Stresslet
d1 = mss_stresslet_direct(src, targ, fvec, nvec, alpha)
d2 = fmm_stresslet_direct(src, targ, fvec, nvec, alpha)
d3 = fmm_stresslet_targ(src, targ, fvec, nvec, alpha)


## Print + test
println(" * Stokeslet")
Gerrdirect = norm(u1[:]-u2[:]) / norm(u1[:], Inf)
Gerrfmm = norm(u1[:]-u3[:]) / norm(u1[:], Inf)
@show Gerrdirect
@show Gerrfmm

println(" * Stresslet")
Terrdirect = norm(d1[:]-d2[:]) / norm(d1[:], Inf)
Terrfmm = norm(d1[:]-d3[:]) / norm(d1[:], Inf)
@show Terrdirect
@show Terrfmm

