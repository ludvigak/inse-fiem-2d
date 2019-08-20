# Generate plots based on data from dev/demo_unsteady.jl

# Speed up loading
module ModifiedStokesSolver
include("../src/types.jl")
end

using FileIO
using HDF5
using PyPlot
using Interpolations

include("../src/plottools.jl")
include("../src/streakplot.jl")

function load_geo(datafile)
    info("Loading ", datafile)        
    D = load(datafile)    
    X = D["X"]
    Y = D["Y"]
    dcurve = D["dcurve"]
    interior = D["interior"]    
    return X, Y, dcurve, interior
end

function make_vortex_plot(;Re=100, vorticity=true, stepnum=15000,
                          streaklines=[], streamlines=false,
                          plotcolorbar=false)
    basedir = "data/vortex_Re$(Re)"
    geoname = "$basedir/geo.jld"
    dataname = i -> @sprintf("%s/vstep_%06d.jld", basedir, i)
    datafile = dataname(stepnum)
    X, Y, dcurve, interior = load_geo(geoname)
    vorticity, U1, U2 = [], [], []
    t = 0.0;
    h5open(datafile, "r") do io
        U1 = read(io, "U1")
        U2 = read(io, "U2")                
        dU2dx = read(io, "dU2dx")
        dU1dy = read(io, "dU1dy")
        t = read(io, "ti")
        vorticity = dU2dx - dU1dy            
    end
    @. vorticity[interior==false] = NaN

    h5open(datafile, "r") do io            
        dU2dx = read(io, "dU2dx")
        dU1dy = read(io, "dU1dy")
        t = read(io, "ti")
        vorticity = dU2dx - dU1dy            
    end
    @. vorticity[interior==false] = NaN        
    clf()
    s = 1        
    plot_boundary(dcurve)
    zoom_to_boundary(dcurve)
    axis("off")
    pcolormesh(X, Y, vorticity, cmap="coolwarm", shading="gouraud")
    clim(-3, 3)
    if plotcolorbar
        colorbar(label="Vorticity", shrink=0.8, pad=0.01)        
    end        
    if streamlines
        # Start points at inflow and in wake
        ystart = 0.1:0.5:3
        ystart = vcat(-reverse(ystart), ystart)
        xstart = zeros(length(ystart)) .- 4.9
        push!(xstart, -2.5)
        push!(ystart, 0.03)
        push!(xstart, -2.5)
        push!(ystart, -0.03)
        
        streamplot(y=Y[1,:],
                   x=X[:,1],
                   u=U1',
                   v=U2',
                   density=40,
                   start_points=hcat( xstart, ystart-0.1),
                   #integration_direction="both",
                   maxlength=4,
                   color="black"
                   
                   )  
    end

    if !isempty(streaklines)
        starty = 0.1:0.5:3
        starty = vcat(reverse(-starty), starty) .- 0.1    
        startx = ones(size(starty)) .* (-5)        
        streakplot(X=X, Y=Y, interior=interior,
                   startx=startx,
                   starty=starty,
                   dataname=dataname,
                   ilist=streaklines)
    end
    
    zoom_to_boundary(dcurve)
    tight_layout(0)

    #pngname = "../docs/paper/figs/vortex_Re$(Re).png"
    #savefig(pngname, dpi=200)
    #run(`mogrify -trim $pngname`)
end

close("all")
figsize = (5.5, 3)

# figure(1, figsize=figsize)
# make_vortex_plot(Re=100, stepnum=6250)

# figure(2, figsize=figsize)
# make_vortex_plot(Re=200, stepnum=14500, plotcolorbar=true)

#figure(3, figsize=figsize)
#sn = 14600
#make_vortex_plot(Re=50, stepnum=sn, plotcolorbar=true, streaklines=sn-4000:1:sn)

figure(4, figsize=figsize)
make_vortex_plot(Re=25, stepnum=10, streamlines=true)
