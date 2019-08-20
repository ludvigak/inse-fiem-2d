using Compat
using Compat.Printf

include("BoxFMMType.jl")
include("BoxFMMUtil.jl")
include("MBHFMM2D.jl")

using PyPlot


## define source and target points

ns = 100
nt = 20

# sources on an ellipse, targets inside

src = Array{Float64}(undef,2,ns)
srcsort = copy(src)
isrcsort = Array{Int32}(undef,ns)
targ = Array{Float64}(undef,2,nt)
targsort = copy(targ)
itargsort = Array{Int32}(undef,nt)

center = [6.1;-1.3]

h = 2.0*pi/ns
a = 1.123
b = 0.7

for i = 1:ns
    t = (i-1)*h
    src[:,i] = center + [a*cos(t);b*sin(t)] 
end

for i = 1:nt
    xt = randn(); yt = randn();
    elnorm = sqrt((xt/a)^2 + (yt/b)^2)
    sc = rand();
    xt = sc*xt/elnorm; yt = sc*yt/elnorm
    targ[:,i] = center + [xt;yt]
end

## build tree based on points

maxboxes = 100000 # length of arrays allocated

# no more than 30 sources and targets per box
tree, sorted_pts, ier = BoxTree2DMaxBoxesST(src,targ,maxboxes,
                                 maxnodes=30)

@printf "built tree with %d boxes and %d levels\n" tree.nboxes tree.nlev

println("sorted_pts ",sorted_pts.srcsort)

@printf "plotting tree and points ...\n"

boxfmm2d_plot(tree)

plot(src[1,:],src[2,:],"bo")
plot(targ[1,:],targ[2,:],"rx")
