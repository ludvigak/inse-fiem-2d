#
# Demo program for timestepping unsteady flows
#
# Currently there's a memory leak in HDF5.jl, so
# this runner will shut down after a certain number
# of steps. To restart automatically, run it as
# while true; do julia output_scripts/demo_unsteady.jl; done


using Revise
using Compat.@info
using JLD
using FileIO
using HDF5

include("timestep_base.jl")

# Bounding box
boxW, boxH = 5, 3

# flow_case = :obstacles
# Re = 30

flow_case = :vortex
Re = 25
offset = true

if flow_case == :obstacles
    basedir = "data/obstacles"
elseif flow_case == :vortex
    if offset
        basedir = "data/vortex_Re$Re"
    else
        basedir = "data/vortex_Re$(Re)_center"
    end
else
    error("Unknown: $(flow_case)")
end
mkpath(basedir)

# Cache for multiple runs
if !isdefined(:solvers)
    solvers = Dict()
end

# Resume timestepping
geoname = "$basedir/geo.jld"
dataname = i -> @sprintf("%s/vstep_%06d.jld", basedir, i)
resume = 0
while isfile(dataname(resume+1))
    resume += 1
end
@show resume

function hook_quiver(solver, sol, ti, i)
    # Unpack
    U1 = sol.U1
    U2 = sol.U2    
    X = solver.X
    Y = solver.Y
    dcurve = solver.dcurve
    interior = solver.interior

    figure(1)
    clf()
    k = 10 # Point skip in quiver plot
    quiver(X[1:k:end, 1:k:end], Y[1:k:end, 1:k:end], U1[1:k:end, 1:k:end], U2[1:k:end, 1:k:end])    
    axis("equal")
end

function hook_save(solver, sol, ti, i)
    @time begin        
        if i==1
            filename = geoname
            println("Writing ", filename)
            save(filename, 
                 "X", solver.X, "Y", solver.Y,
                 "interior", solver.interior,
                 "dcurve", solver.dcurve
                 )
        end
        filename = dataname(i)
        println("Writing ", filename)
        U1    = sol.U1   
        U2    = sol.U2   
        dU1dx = sol.dU1dx
        dU1dy = sol.dU1dy
        dU2dx = sol.dU2dx
        dU2dy = sol.dU2dy
        save(File(format"JLD", filename),
             "ti", ti,
             "U1"    , U1   ,
             "U2"    , U2   ,
             "dU1dx" , dU1dx,
             "dU1dy" , dU1dy,
             "dU2dx" , dU2dx,
             "dU2dy" , dU2dy,
             "RHS1"  , sol.RHS1,
             "RHS2"  , sol.RHS2
             )
    end
end
    
function hook_vortex(solver, sol, ti, i)
    hook_save(solver, sol, ti, i)    
    if mod(i, 1)==0
        println("Vorticity plot:")
        @time hook_vorticity(solver, sol, ti, i)
    end
    if mod(i, 20)==0
        hook_streamlines(solver, sol, ti, i)
    end
end

function hook_vorticity(solver, sol, ti, i)
    figure(1)
    clf()
    X = solver.X
    Y = solver.Y        
    vorticity = sol.dU2dx - sol.dU1dy
    s = 1
    pcolormesh(X[1:s:end,1:s:end], Y[1:s:end,1:s:end], vorticity[1:s:end,1:s:end],
           cmap=PyPlot.cm[:coolwarm])
    clim(-3,3)
    plot_boundary(solver.dcurve)        
    axis("image")
    axis("off")    
end

function hook_streamlines(solver, sol, ti, i)
    figure(2)
    clf()
    U1 = sol.U1
    U2 = sol.U2    
    X = solver.X
    Y = solver.Y
    dcurve = solver.dcurve
    interior = solver.interior    
    ystart = collect(linspace(-boxH+0.05,boxH-0.05,40))    
    xstart = -boxW+0.1 + zeros(ystart)
    start_points = [xstart ystart]
    streamplot(y=Y[1,:],
               x=X[:,1],
               u=U1',
               v=U2',
               density=20,
               start_points=start_points,
               integration_direction="both",
               maxlength=6)       
    
    plot_boundary(dcurve)    
    axis("image")
    axis("off")
end    

    

function demo_unsteady(flow_case, Re, resume)
    # Defaults
    fastdirect_TOL = 1e-8
    
    # Geometry setup
    if flow_case == :obstacles        
        numpanels1 = 500
        numpanels2 = 50
        Ngrid = 500
        dt = 0.01
        sdc_order = 2
        pux_R = 0.4
        Nstep = 10000

        curve1 = AnalyticDomains.rectangle(width=boxW, height=boxH)
        # Setup a bunch of starfish
        srand(1)
        curve2 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = (-3.0, -1.0), rotation=2*pi*rand(), interior=true)
        curve3 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = (-1.0, -0.1), rotation=2*pi*rand(), interior=true)
        curve4 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = (-0.9,  1.5), rotation=2*pi*rand(), interior=true)
        curve5 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = ( 1.0,  0.0), rotation=2*pi*rand(), interior=true)
        curve6 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = (-3.5,  0.7), rotation=2*pi*rand(), interior=true)
        curve7 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = ( 1.7, -2.0), rotation=2*pi*rand(), interior=true)
        curve8 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = (-0.7, -1.9), rotation=2*pi*rand(), interior=true)
        curve9 = AnalyticDomains.starfish(;amplitude = 0.2, radius = 0.5, center = ( 2.9,  1.7), rotation=2*pi*rand(), interior=true)
        
        curve = [curve1, curve2, curve3, curve4, curve5, curve6, curve7, curve8, curve9]
        numpanels = [numpanels1, numpanels2, numpanels2, numpanels2, numpanels2, numpanels2, numpanels2, numpanels2, numpanels2]
    elseif flow_case == :vortex       
        if Re==200
            # Stable, but slow
            # Oscillations around 2500
            Re = 200
            dt = 0.01/2
            Ngrid = 500*2
            pux_R = 0.5/2
        elseif Re==100
            # Separates after ~2000 steps with offset
            Re = 100
            dt = 0.01
            Ngrid = 500
            pux_R = 0.5
        elseif Re==50
            # offset and SDC=1: 3500 steps oscillations visible
            Re = 50
            dt = 0.02
            Ngrid = 500
            pux_R = 0.5
        elseif Re==25
            # No oscillations after 10000 steps
            Re = 25
            dt = 1/Re
            Ngrid = 500
            pux_R = 0.5
        end        
        
        sdc_order = 2
        Nstep = 20000
        
        numpanels1 = 500
        numpanels2 = 100
        if offset
            cylinder = (-boxW+1.5, -0.1) # Offset triggers faster separation
        else
            cylinder = (-boxW+1.5, 0.0) # No offset
        end
        
        curve1 = AnalyticDomains.rectangle(width=boxW, height=boxH)
        curve2 = AnalyticDomains.starfish(;amplitude = 0.0, radius = 0.5, center = cylinder, interior=true)
        curve = [curve1, curve2]
        numpanels = [numpanels1, numpanels2]
    else
        error()
    end

    # Restart every NN steps (HDF5.jl memory leak!)        
    restart = Ngrid > 500 ? 300 : 1000
    Nstep = min(Nstep, resume+restart)
    @info("Will run steps $resume:$Nstep on flow case $flow_case")    

    # Solver precomputation
    dt_substep = sdc_substep(dt, sdc_order)
    alpha = sqrt(Re/dt_substep)
    @info("alpha = $alpha")
    if haskey(solvers, alpha)
        println(" * USING CACHE")
        solver = solvers[alpha]
    else
        solver = precompute_solver(alpha=alpha,
                                   pux_R=pux_R,
                                   Ngrid=Ngrid,
                                   curve = curve,
                                   numpanels = numpanels,
                                   fastdirect_TOL = fastdirect_TOL=1e-8
                                   )
        solvers[alpha] = solver
    end

    # Set BCs
    dcurve = solver.dcurve
    ubc1 = zeros(dcurve.numpoints)
    ubc2 = zeros(dcurve.numpoints)    
    if flow_case == :obstacles
        # Flow on outer boundary
        @. ubc1[dcurve.curvenum==1] = 1.0
    elseif flow_case == :vortex
        # # Set flow rate on left/right box boundaries
        # idx1 = 1:dcurve.panelorder*numpanels1
        # xbdry = dcurve.points[1,idx1]
        # ybdry = dcurve.points[2,idx1]
        # sigmax = boxW/10
        # u1bdry(x,y) = (exp(-((x-boxW)/sigmax).^2) + exp(-((x+boxW)/sigmax).^2))
        # ubc1[idx1] = u1bdry.(xbdry,ybdry)
        # @. ubc1[idx1] = 1.0

        # Flow on outer boundary
        @. ubc1[dcurve.curvenum==1] = 1.0
    else
        error("Don't know what BCs to set")
    end

    # Estimate BC resolution
    panelorder = dcurve.panelorder
    L = mss.legendre_matrix(panelorder)
    bc_resolution = 0.0
    for i=1:dcurve.numpanels
        xcoeff = L*ubc1[(1:panelorder) + panelorder*(i-1)]
        ycoeff = L*ubc2[(1:panelorder) + panelorder*(i-1)]
        xres = max(abs(xcoeff[end]), abs(xcoeff[end-1]))
        yres = max(abs(ycoeff[end]), abs(ycoeff[end-1]))
        res = max(xres, yres)
        bc_resolution = max(bc_resolution, res)
    end
    info("BC resolution: ", bc_resolution)    

    # Start from rest
    X = solver.X
    U1, U2 = zeros(X), zeros(X)
    RHS1, RHS2 = zeros(X), zeros(X)
    dU1dx, dU1dy = zeros(X), zeros(X)
    dU2dx, dU2dy = zeros(X), zeros(X)    
    sol0 = Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy, RHS1, RHS2, [], [])
    
    
    t0 = 0.0
    istart = 0

    hook = hook_vortex

    if resume > 0
        info("Resuming from ", resume)
        istart = resume
        filename = dataname(istart)
        h5open(filename, "r") do io
            U1 = read(io, "U1")            
            U2 = read(io, "U2")
            dU1dx = read(io, "dU1dx")
            dU1dy = read(io, "dU1dy")            
            dU2dx = read(io, "dU2dx")
            dU2dy = read(io, "dU2dy")            
            RHS1 = read(io, "RHS1")
            RHS2 = read(io, "RHS2")
            sol0 = Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy, RHS1, RHS2, [], [])            
        end   
    end

    t, sol = timestep(Nstep, dt, sdc_order, hook, solver, t0, sol0, ubc1, ubc2; istart=istart)
    
end

demo_unsteady(flow_case, Re, resume)
