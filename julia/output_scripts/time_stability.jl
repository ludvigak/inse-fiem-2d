#
# Demo for unsteady flows
#
using Revise
using Compat.@info
using JLD
using FileIO
using HDF5

include("timestep_base.jl")

size_with_colorbar = [4,3]
size_no_colorbar = [3,3]

function hook_plot(solver, sol, ti, i)
    @time begin
        figure(2)
        clf()
        hook_vorticity(solver, sol, ti, i)
        hook_quiver(solver, sol, ti, i)        
        println("* Plotting done")
    end
    show()
end

function hook_quiver(solver, sol, ti, i)
    # Unpack
    U1 = sol.U1
    U2 = sol.U2    
    X = solver.X
    Y = solver.Y
    dcurve = solver.dcurve
    interior = solver.interior

    @. U1[interior==false] = NaN
    @. U2[interior==false] = NaN
    
    k = 20*2 # Point skip in quiver plot
    quiver(X[1:k:end, 1:k:end], Y[1:k:end, 1:k:end], U1[1:k:end, 1:k:end], U2[1:k:end, 1:k:end])    
    plot_boundary(solver.dcurve)        
    axis("image")
    axis("off")
end

function hook_vorticity(solver, sol, ti, i)
    X = solver.X
    Y = solver.Y        
    vorticity = sol.dU2dx - sol.dU1dy
    interior = solver.interior
    @. vorticity[interior==false] = NaN
    s = 1
    pcolormesh(X[1:s:end,1:s:end], Y[1:s:end,1:s:end], vorticity[1:s:end,1:s:end],
               cmap=PyPlot.cm[:coolwarm])
    clim(-3,3)
    if figure(2).o[:get_size_inches]() == size_with_colorbar
        colorbar(shrink=0.8, pad=0, label="Vorticity")
    end
    plot_boundary(solver.dcurve)        
    axis("image")
    axis("off")
end

function precompute(Re = 20, dt=0.1; R=0.5, Ngrid=150, pux_R_over_R=1/4, numpanels=100, sdc_order=1)
    dt_substep = sdc_substep(dt, sdc_order)
    alpha = sqrt(Re/dt_substep)
    @info("alpha = $alpha")
    figure(1)
    circle = AnalyticDomains.starfish(n_arms=0, amplitude=0.0, radius=R)        
    solver = precompute_solver(alpha=alpha,
                               pux_R=R*pux_R_over_R,
                               Ngrid=Ngrid,
                               curve = circle,
                               numpanels = numpanels
                               )
    return solver, dt, sdc_order
end

function timestep(pre, Nstep=10; velocity=1)
    solver, dt, sdc_order = pre
    # Set BCs
    dcurve = solver.dcurve
    xbdry = dcurve.points[1,:]
    ybdry = dcurve.points[2,:]
    n1 = dcurve.normals[1,:]
    n2 = dcurve.normals[2,:]    
    phi = (x, y) -> atan2(y, x)
    k = 1
    ubc1 = @. -n1*sin(k*phi(xbdry, ybdry))*velocity
    ubc2 = @. -n2*sin(k*phi(xbdry, ybdry))*velocity
    
    compat = @. (n1*ubc1 + n2*ubc2)*dcurve.dS
    info("Compability condition=", sum(compat))
    assert(sum(compat) < 1e-10)
    
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
    Y = solver.Y    
    U1, U2 = zeros(X), zeros(X)
    RHS1, RHS2 = zeros(X), zeros(X)
    dU1dx, dU1dy = zeros(X), zeros(X)
    dU2dx, dU2dy = zeros(X), zeros(X)  
    
    sol0 = Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy, RHS1, RHS2, [], [])
    # Run   
    t0 = 0.0
    istart = 0
    hook = hook_plot
    t, sol = timestep(Nstep, dt, sdc_order, hook, solver, t0, sol0, ubc1, ubc2; istart=istart)    
end

#  pre = precompute()
# timestep(pre)



function instab1()
    close("all")
    figure(2, figsize=size_no_colorbar)
    pre = precompute(50, 0.04; R=0.5, Ngrid=150*2, pux_R_over_R=1/4, numpanels=100*2, sdc_order=1)
    timestep(pre, 45, velocity=1)
end

function instab2()
    close("all")
    figure(2, figsize=size_no_colorbar)
    pre = precompute(50, 0.02; R=0.5, Ngrid=150*2, pux_R_over_R=1/4, numpanels=100*2, sdc_order=1)
    timestep(pre, 200, velocity=1)
end

function instab3()
    close("all")
    figure(2, figsize=size_no_colorbar)
    pre = precompute(100, 0.02; R=0.5, Ngrid=150*2, pux_R_over_R=1/4, numpanels=100*2, sdc_order=1)
    timestep(pre, 45, velocity=1)
end

function instab4()
    close("all")
    figure(2, figsize=size_with_colorbar)
    pre = precompute(100, 0.01; R=0.5, Ngrid=150*2, pux_R_over_R=1/4, numpanels=100*2, sdc_order=1)
    timestep(pre, 200, velocity=1)
end

function write_figure(num)
    pngname = "../docs/paper/figs/time_stability_$(num).png"
    savefig(pngname, dpi=200)
    Base.run(`mogrify -trim $pngname`)
end    

function paper_scripts()
    instab1()
    write_figure(1)
    instab2()
    write_figure(2)
    instab3()
    write_figure(3)
    instab4()
    write_figure(4)
end
