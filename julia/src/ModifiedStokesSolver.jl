__precompile__()

module ModifiedStokesSolver

using Compat
using Compat.Libdl
using Compat.LinearAlgebra
using Compat.MathConstants
using Compat.@info


using PyCall

using CurveDiscretization
using BarycentricLagrange

export DiscreteCurve
export discretize, multiply_connected
export doublelayer, doublelayer_fast, doublelayer_self, doublelayer_matrix
export system_matrix, system_matvec
export correction_block_matrix
export stokeslet
export pressure
export stresslet, stresslet_split
export dblkernel, dblkernel_log, dblkernel_block!
export fmm_stokeslet_direct, fmm_stokeslet_targ
export fmm_stresslet_direct, fmm_stresslet_prep
export fmm_stresslet_targ, fmm_stresslet_self
export ndgrid, interior_points
export fs_modstokes, per_modstokes, per_gradient

include("../../external/mbh2dfmm/julia/ModStokesFMM2D.jl")
include("types.jl")
include("helsing.jl")
include("layerpotential.jl")
include("integralequation.jl")
include("kernels.jl")
include("kernel_gradient.jl")
include("legendre.jl")
include("sparse.jl")
include("freespace.jl")
include("periodic.jl")
include("matlab_compat.jl")
include("largealpha.jl")

function __init__()
    # Load libgomp for FMM
    Libdl.dlopen("libgomp", Libdl.RTLD_GLOBAL)
end

end
