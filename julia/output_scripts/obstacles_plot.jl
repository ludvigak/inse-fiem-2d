#
# Plot results for demo_unsteady.jl with obstacles flow case
#

# Speed up loading
module ModifiedStokesSolver
include("../src/types.jl")
end

using Compat
using Compat.@info
using JLD
using HDF5
using PyPlot
using Interpolations

include("../src/plottools.jl")

# Geometry
boxW, boxH = 5, 3

basedir = "data/obstacles"
geoname = "$basedir/geo.jld"
dataname = i -> @sprintf("%s/vstep_%06d.jld", basedir, i)

function load_geo()
    datafile = geoname
    info("Loading ", datafile)        
    D = load(datafile)    
    X = D["X"]
    Y = D["Y"]
    dcurve = D["dcurve"]
    interior = D["interior"]    
    return X, Y, dcurve, interior
end

function vorticity()
    X, Y, dcurve, interior = load_geo()

    ilist = [10000]

    stepnum = 0
    for i=ilist
        stepnum += 1
        datafile = dataname(i)
        info("Loading ", datafile)
        U1, U2, vorticity = [], [], []
        h5open(datafile, "r") do io
            U1 = read(io, "U1")
            U2 = read(io, "U2")                        
            dU2dx = read(io, "dU2dx")
            dU1dy = read(io, "dU1dy")
            vorticity = dU2dx - dU1dy            
        end                
        figure(1, figsize=[boxW, boxH].*2)
        clf()
        s = 1
        axis("off")        
        plot_boundary(dcurve)
        zoom_to_boundary(dcurve)
        @. vorticity[ (U1==0) & (U2==0) ] = NaN
        pcolormesh(X[1:s:end,1:s:end], Y[1:s:end,1:s:end], vorticity[1:s:end,1:s:end],
                   cmap=PyPlot.cm[:coolwarm], shading="gouraud")
        clim(-6,6)
        colorbar(label="Vorticity", shrink=0.8, pad=0.03)        
        zoom_to_boundary(dcurve)
        ystart = collect(linspace(-boxH+0.05,boxH-0.05,40))    
        xstart = -boxW+0.7 + zeros(ystart)

        #xstart, ystart = [],[]
        
        function addpt(r1,r2)
            #plot(r1, r2, "ok")
            push!(xstart, r1)
            push!(ystart, r2)
        end

        # top-bottom, left-right
        
        addpt(-2.4, 0.7)
        addpt(-2.7, 1.0)
       
        addpt(-2.0, -1.0)
        addpt(-2.1, -0.4)

        addpt(0.0, 1.6)

        addpt(-0.2, 0.2)

        addpt(-0.2, -1.5)
        addpt(0.6, -2)

        addpt(3.7, 1.3)
        addpt(3.9, 2.0)
        addpt(1.8, -0.35)

        addpt(2.5, -1.5)        
        addpt(2.5, -2.3)        


        # Random
        Nrand = 2000
        xstart = boxW*(2*rand(Nrand)-1)
        ystart = boxH*(2*rand(Nrand)-1)

        s = 3
        xstart = vec(X[1:s:end,1:s:end])
        ystart = vec(Y[1:s:end,1:s:end])
        
        vel = @. sqrt(U1^2+U2^2)        
        start_points = [xstart ystart]
        streamplot(y=Y[1,:],
                   x=X[:,1],
                   u=U1',
                   v=U2',
                   density=4,
                   start_points=start_points,
                   integration_direction="both",
                   maxlength=10,
                   color="k"
                   )
        #tight_layout(0)
        show()
        pause(0.1/2)
        savefig("../docs/paper/figs/demo.png", dpi=200)
        run(`mogrify -trim ../docs/paper/figs/demo.png`)
    end
end

vorticity()
