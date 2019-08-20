# Plots of quadrature geometry

push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using PyPlot
import IterativeSolvers
import LinearMaps

import ModifiedStokesSolver
import AnalyticDomains

mss = ModifiedStokesSolver

# Params
alpha = 200

# Discretize
numpanels = 20
panelorder = 16
curve = AnalyticDomains.starfish(n_arms = 3, amplitude=0.1)
dcurve = mss.discretize(curve, numpanels, panelorder)

N = dcurve.numpoints
panel_length = sum(dcurve.dS[1:panelorder])

n_sub = 32
dt_max = mss.LARGEALPHA_LIMIT_16*2/(alpha*panel_length)    
glpoints1, glweights1 = mss.gausslegendre(panelorder)
glpoints2, glweights2 = mss.gausslegendre(n_sub)


######################## Near boundary point
panel_idx = 1
point_idx = 6

t1 = dcurve.t_edges[:,panel_idx]
idx = panelorder*(panel_idx-1) + (1:panelorder)  
xpanel = dcurve.points[1,idx]
ypanel = dcurve.points[2,idx]
zj = xpanel + 1im*ypanel
nj = dcurve.normals[1,idx] + 1im*dcurve.normals[2,idx]

h1 = t1[2]-t1[1]
@show t = glpoints1[point_idx] + 0.08im;
ztmp = curve.tau( t1[1] + (t+1)*h1/2 )
zt = [real(ztmp), imag(ztmp)]

coeffs = CurveDiscretization.map_panels(dcurve)
zcpx = zt[1] + 1im*zt[2]
troot, _, converged =
    CurveDiscretization.invert_map(dcurve, coeffs, panel_idx, zcpx,
                                   fall_back_to_initial=false)
@show troot

paperplots = true
if paperplots
    close("all")
    figure(1, figsize=(3,3))
else
    figure(1)
    clf()
end

plot(dcurve.points[1,:], dcurve.points[2,:], "-k", label=L"\Gamma")

cpx2vec(z) = [real(z),imag(z)]
function plot_edge(x, n, h, style="b")
    plot(x[1]+n[1]*[-1,1]*h, x[2]+n[2]*[-1,1]*h, style)
end

h1 = 0.06*panel_length
h2 = 0.03*panel_length

for i=1:numpanels
    x = dcurve.edges[1:2,i]
    n = dcurve.n_edges[1:2,i]
    plot_edge(x, n, h1)
end

plot(zt[1], zt[2], ".", label=L"x \in \Omega")

# Plot subintervals
bcweights = mss.bclag_interp_weights(glpoints1)    
@show intervals = mss.subdivide_interval_with_bisection(dt_max, troot, n_sub)
#intervals = mss.subdivide_interval(dt_max, troot)
for idx = 1:length(intervals)-1
    ta = intervals[idx]
    tb = intervals[idx+1]
    troot_sec = (troot-ta)*2/(tb-ta) - 1.0
    # New source points, in original frame
    tsrc = ta + (1+glpoints2)*(tb-ta)/2
    # Interpolation to new source points
    Psec = mss.bclag_interp_matrix(glpoints1, tsrc, bcweights)
    zj_sec = Psec*zj
    za_sec = mss.bclag_interp_eval(glpoints1, zj, ta, bcweights)
    zb_sec = mss.bclag_interp_eval(glpoints1, zj, tb, bcweights)
    nj_sec = Psec*nj
    na_sec = mss.bclag_interp_eval(glpoints1, nj, ta, bcweights)
    nb_sec = mss.bclag_interp_eval(glpoints1, nj, tb, bcweights)

    if idx>1
        plot_edge(cpx2vec(za_sec), cpx2vec(na_sec), h2, "r")
    end
    if idx<length(intervals)-1
        plot_edge(cpx2vec(zb_sec), cpx2vec(nb_sec), h2, "r")
    end
end   

xticks([])
yticks([])
minorticks_off()
axis("equal")
axis((0.9153349805587189, 1.2063410174030949, -0.02194527272134522, 0.3516612032549363))
legend()

function write_fig(num, name)
    figure(num)
    savefig(name)
    println("Saved $name")
end

filename = joinpath(dirname(@__DIR__),"figs", "quad_geo_near.pdf")
write_fig(1, filename)
run(`cp $filename ../docs/paper/figs`)

######################## On boundary point
zt = dcurve.points[:,point_idx]

coeffs = CurveDiscretization.map_panels(dcurve)
zcpx = zt[1] + 1im*zt[2]
troot = glpoints1[point_idx]

paperplots = true
if paperplots
    figure(2, figsize=(3,3))
else
    figure(2)
    clf()
end

plot(dcurve.points[1,:], dcurve.points[2,:], "-k", label=L"\Gamma")

for i=1:numpanels
    x = dcurve.edges[1:2,i]
    n = dcurve.n_edges[1:2,i]
    plot_edge(x, n, h1)
end

plot(zt[1], zt[2], ".", label=L"x \in \Gamma")

# Plot subintervals
intervals = mss.subdivide_interval_with_bisection(dt_max, troot, 32)
#intervals = mss.subdivide_interval(dt_max, troot)
for idx = 2:length(intervals)-2
    ta = intervals[idx]
    tb = intervals[idx+1]
    troot_sec = (troot-ta)*2/(tb-ta) - 1.0
    # New source points, in original frame
    tsrc = ta + (1+glpoints2)*(tb-ta)/2
    # Interpolation to new source points
    Psec = mss.bclag_interp_matrix(glpoints1, tsrc, bcweights)
    zj_sec = Psec*zj
    za_sec = mss.bclag_interp_eval(glpoints1, zj, ta, bcweights)
    zb_sec = mss.bclag_interp_eval(glpoints1, zj, tb, bcweights)
    nj_sec = Psec*nj
    na_sec = mss.bclag_interp_eval(glpoints1, nj, ta, bcweights)
    nb_sec = mss.bclag_interp_eval(glpoints1, nj, tb, bcweights)

    plot_edge(cpx2vec(za_sec), cpx2vec(na_sec), h2, "r")
    plot_edge(cpx2vec(zb_sec), cpx2vec(nb_sec), h2, "r")    
end   

xticks([])
yticks([])
minorticks_off()
axis("equal")
axis((0.9153349805587189, 1.2063410174030949, -0.02194527272134522, 0.3516612032549363))
legend()

filename = joinpath(dirname(@__DIR__),"figs", "quad_geo_on.pdf")
write_fig(2, filename)
run(`cp $filename ../docs/paper/figs`)
