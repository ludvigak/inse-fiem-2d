using NullableArrays
using PyCall
@pyimport numpy.ma as ma
PyObject(a::NullableArray) =  pycall(ma.array, Any, a.values, mask=a.isnull)

function interior_pcolor(xx, yy, ff, interior_mask; vmin=nothing, vmax=nothing, cmap=nothing)
    P = NullableArray(ff, .!reshape(interior_mask, size(ff)))
    pcolormesh(xx, yy, PyObject(P), vmin=vmin, vmax=vmax, cmap=cmap)
end

function plot_boundary(dcurve, style="-k")
    Ncurves = dcurve.curvenum[end]
    for n=1:Ncurves
        idx = findin(dcurve.curvenum, n)    
        plot(dcurve.points[1,idx], dcurve.points[2,idx], style)
    end
    axis("equal")
end

function zoom_to_boundary(dcurve)
    xmin, xmax = extrema(dcurve.points[1,:]) .+ (-0.01, 0.01)
    ymin, ymax = extrema(dcurve.points[2,:]) .+ (-0.01, 0.01)
    axis((xmin, xmax, ymin, ymax))
end
