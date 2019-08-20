#
# Demo and convergences test for the spindown problem
#

# Repeated runner:
# while true; do julia dev/demo_spindown.jl; done

using Compat.@info
using LaTeXStrings
using JLD
using Glob

include("timestep_base.jl")
include("spindown_solution.jl")

# Cache
if !isdefined(:solvers)
    solvers = Dict()
end


global TOTAL_STEPS_TAKEN = 0 
global MAX_TOTAL_STEPS = 1000 # HDF5 bug forces us to limit this

function precompute_spindown(dt, sdc_order, Re; dummy=false)
    dt_substep = sdc_substep(dt, sdc_order)
    alpha = sqrt(Re/dt_substep)
    @info("alpha = $alpha")
    if dummy
        method = :none
    else
        method = :heuristic
    end
    return precompute_solver(alpha=alpha, pux_R=0.3, Ngrid=500, numpanels=200, method=method)
end

function spindown_ubc(solver)
    dcurve = solver.dcurve
    ubc1 = zeros(dcurve.numpoints)
    ubc2 = zeros(dcurve.numpoints)
    return ubc1, ubc2
end
    
function spindown_initialize(t0, solver, Re)
    # Unpack
    X = solver.X
    Y = solver.Y
    dcurve = solver.dcurve
    interior = solver.interior
    PUXStor = solver.PUXStor
    # Initial information 
    U1, U2 = zeros(X), zeros(X)
    L1, L2 = zeros(X), zeros(X)
    dU1dx, dU1dy = zeros(X), zeros(X)
    dU2dx, dU2dy = zeros(X), zeros(X)    
    @time for i=1:length(X)
        if interior[i]
            r = sqrt(X[i]^2+Y[i]^2)
            if r > 0.0
                u = u_spindown(r, t0, Re)
                Du = Du_spindown(r, t0, Re)                
                Lu = Lu_spindown(r, t0, Re)            
                U1[i] = -Y[i]/r*u
                U2[i] = X[i]/r*u
                L1[i] = -Y[i]/r*Lu
                L2[i] = X[i]/r*Lu
                dU1dx[i] = X[i]*Y[i]/r^2*(u/r-Du)
                dU1dy[i] = Y[i]^2/r^2*(u/r-Du)-u/r
                dU2dx[i] = -X[i]^2/r^2*(u/r-Du)+u/r                
                dU2dy[i] = -X[i]*Y[i]/r^2*(u/r-Du)                
            end
        end
    end
    #dU1dx, dU1dy, dU2dx, dU2dy = pux_gradient(U1, U2, PUXStor, solver.    
    # SDC needs initial Laplacian, which it recovers from RHS
    # -> Reverse engineer an approx of RHS (pressure missing, but thats ok)
    RHS1 = solver.alpha^2*U1 - L1
    RHS2 = solver.alpha^2*U2 - L2
    sol = Solution(U1, U2, dU1dx, dU1dy, dU2dx, dU2dy, RHS1, RHS2, [], [])
    return sol
end

function hook_spindown(solver, sol, ti, i)
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
    
    # Plot velocities
    figure(3)
    clf()
    
    #mid = Int64(round(Ngrid/2)+1)
    #plot(X[:, mid], U1[:, mid], ".", label="U1")
    #plot(X[:, mid], U2[:, mid], ".", label="U2")

    r = sqrt.(vec(X[interior]).^2 + vec(Y[interior]).^2)
    r[r .== 0] = 1
    u = sqrt.(vec(U1[interior]).^2 + vec(U2[interior]).^2)

    u = (vec(U1[interior]).*vec(-Y[interior]) + vec(U2[interior]).*vec(X[interior]))./r
    
    plot(r, u, ".", label="Solution")
    
    #plot_spindown(ti, Re)    
    
    xlim(0,1)
    ylim(0,1)
    #legend()
    grid("on")

    # Plot error
    figure(4)
    clf()
    uref = map(x -> u_spindown(x, ti, Re), r)
    plot(r, (u-uref), ".", label="Error")
    legend()
    grid("on")
    show()
    
    maxerr = norm(u-uref, Inf)
    @show maxerr

    l2err = sqrt( sum( (u-uref).^2)/length(u) )
    @show l2err
end

function spindown_error(solver, sol, ti, Re)
    println(" * Computing error at t=$ti")
    U1 = sol.U1
    U2 = sol.U2    
    X = solver.X
    Y = solver.Y
    dcurve = solver.dcurve
    interior = solver.interior
    r = sqrt.(vec(X[interior]).^2 + vec(Y[interior]).^2)
    r[r .== 0] = 1
    #u = sqrt.(vec(U1[interior]).^2 + vec(U2[interior]).^2)
    u = (vec(U1[interior]).*vec(-Y[interior]) + vec(U2[interior]).*vec(X[interior]))./r
    uref = map(x -> u_spindown(x, ti, Re), r)
    unorm = norm(uref, Inf)

    maxnorm(v) = maximum(abs.(v))
    l2norm(v) = sqrt(sum(vec(v.^2)))

    maxerr = maxnorm(u-uref) / maxnorm(uref)
    l2err = l2norm(u-uref) / l2norm(uref)
    
    @show maxerr
    @show l2err

    return maxerr, l2err
end

function spindown_basename(Re, t0, dt, sdc_order)
    basepath = @sprintf("data/spindown/Re%.g/init%.g/dt%.g/sdc%d",
                        Re, t0, dt, sdc_order)
    return basepath    
end

function spindown_stepname(i)
    return @sprintf("step_%06d.jld", i)
end

function spindown_save(Re, t0, dt, sdc_order)
    basepath = spindown_basename(Re, t0, dt, sdc_order)
    mkpath(basepath)
    function save_hook(solver, sol, ti, i)
        filename = joinpath(basepath, spindown_stepname(i))
        save(filename, "sol", sol, "ti", ti)
        println("Saved $filename")
    end
    return save_hook
end

function spindown_resume(Nstep, Re, t0, dt, sdc_order)
    basepath = spindown_basename(Re, t0, dt, sdc_order)
    if isdir(basepath)
        stepnum = 0
        filename = joinpath(spindown_basename(Re, t0, dt, sdc_order),
                            spindown_stepname(Nstep))
        if isfile(filename)
            stepnum = Nstep
        else
            files = glob("step_*.jld", basepath)
            if !isempty(files)
                sort!(files)
                lastfile = files[end]
                m = match(r"^step_(\d+)\.jld$", basename(lastfile))
                stepnum = parse(Int64, m.captures[1])
                filename = lastfile
            end
        end
        if stepnum > 0
            info("Resuming at step $stepnum")
            println("Loading $filename")
            D = load(filename)
            return true, stepnum, D["ti"], D["sol"]
        end            
    end
    return false, 0, 0, 0
end


function demo_spindown(Re, t0, iterations, Nstep, dt, sdc_order)
    # Convergence test

    #hook = hook_null

    @show iterations
    @show l2err = zeros(iterations)
    maxerr = zeros(iterations)
    dtlist = zeros(iterations)

    dummy = precompute_spindown(dt, sdc_order, Re, dummy=true)

    figure(1)
    for k=1:iterations
        println("=================================")
        @show k

        # Check if resume possible
        do_resume, i, t_i, sol_i = spindown_resume(Nstep, Re, t0, dt, sdc_order)
        if do_resume
            istart = i
            sol = sol_i
        else
            istart = 0
        end

        t = t0 + istart*dt

        if istart < Nstep
            # Setup solver
            alpha = sqrt(Re/sdc_substep(dt, sdc_order))
            @show alpha   
            if haskey(solvers, alpha)
                println(" * USING CACHE")
                solver = solvers[alpha]
            else
                solver = precompute_spindown(dt, sdc_order, Re)
                solvers[alpha] = solver
            end
            ubc1, ubc2 = spindown_ubc(solver)
            
            # Initial conditions
            if istart==0
                println("* Compute initial conditions")
                tic()            
                sol = spindown_initialize(t0, solver, Re)
                toq()
            end

            # Save hook
            hook = spindown_save(Re, t0, dt, sdc_order)        

            # Check if we're past step limit
            global TOTAL_STEPS_TAKEN
            global MAX_TOTAL_STEPS
            steps_left = Nstep - istart
            if TOTAL_STEPS_TAKEN + steps_left > MAX_TOTAL_STEPS
                # Truncate how many steps we're taking on this pass
                Nstep = istart + MAX_TOTAL_STEPS - TOTAL_STEPS_TAKEN
            end
            
            # Do the actual timestepping
            t, sol = timestep(Nstep, dt, sdc_order, hook, solver, t0, sol, ubc1, ubc2; istart=istart)

            TOTAL_STEPS_TAKEN += Nstep - istart
            @show TOTAL_STEPS_TAKEN
            if TOTAL_STEPS_TAKEN >= MAX_TOTAL_STEPS
                println("* Out of steps!")
                quit()
            end
            
        else
            # Dummy solver
            solver = dummy
        end

        @show spindown_error(solver, sol, t, Re)
        @show maxerr
        @show l2err
        
        maxerr[k], l2err[k] = spindown_error(solver, sol, t, Re)
        dtlist[k] = dt

        figure(99)
        clf()
        loglog(dtlist[1:k], maxerr[1:k], ".-")
        loglog(dtlist[1:k], l2err[1:k], ".-")
        loglog(dtlist[1:k], (dtlist[1:k]./dtlist[1]).^sdc_order * maxerr[1], "--")
        axis("tight")
        figure(1)
        
        Nstep *= 2
        dt /= 2
    end

    return dtlist, l2err, maxerr
end

# Convergence test
Re = 1
t0 = 0.05*Re # Very smooth start
dt = 0.004*Re
Nstep = 5
iterations = 8

dtlist1, l2err1, maxerr1 = demo_spindown(Re, t0, iterations, Nstep, dt, 1)
dtlist2, l2err2, maxerr2 = demo_spindown(Re, t0, iterations, Nstep, dt, 2)
dtlist3, l2err3, maxerr3 = demo_spindown(Re, t0, iterations, Nstep, dt, 3)
dtlist4, l2err4, maxerr4 = demo_spindown(Re, t0, iterations, Nstep, dt, 4)

include("demo_spindown_plot.jl")
