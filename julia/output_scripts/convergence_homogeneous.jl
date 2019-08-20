# Dev code for simple convergence test, copied out from dev.jl

push!(LOAD_PATH, string(pwd(),"/src"))

using Revise
using PyPlot
import IterativeSolvers
import LinearMaps

using ModifiedStokesSolver
import AnalyticDomains
using PUX

# Params
alpha = 10.0
numpanels_list = vcat(25,30:10:120)
ticks = [30, 60, 120]

panelorder = 16
curve = AnalyticDomains.starfish()

# Setup a problem
srand(1)
Nsrc = 5
Rsrc = 1.4

theta = 2*pi*rand(Nsrc)
src = Rsrc*[cos.(theta) sin.(theta)]'
str = 1-2*rand(2, Nsrc)
ubc(targ) = fmm_stokeslet_targ(src, targ, str, alpha)
ubc_grad(targ) = fmm_stokeslet_targ(src, targ, str, alpha, ifgrad=true)

# Ntrg = 50
# theta = 2*pi*rand(Ntrg)
# zt_interior = 0.65*(rand(Ntrg).*[cos.(theta) sin.(theta)])'

Ngrid = 30
dcurve = discretize(curve, 500, panelorder)
x = linspace(-1.4, 1.4, Ngrid)
X, Y = ModifiedStokesSolver.ndgrid(x, x)
zt = [vec(X) vec(Y)]'
interior = ModifiedStokesSolver.interior_points(dcurve, zt)
zt_interior = zt[:, interior]

function solve(numpanels, method=:dense)
    # Discretize
    @show numpanels
    dcurve = discretize(curve, numpanels, panelorder)
    panel_length = sum(dcurve.dS[1:panelorder])
    @show panel_length*alpha
    
    N = dcurve.numpoints
    # Right hand side
    rhs = ubc(dcurve.points)
    # Left hand side
    println("* Assembly")
    if method==:dense
        LHS = system_matrix(dcurve, alpha)
    elseif method==:fmm
        flhs = system_matvec(dcurve, alpha)
        LHS = LinearMaps.LinearMap(flhs, 2*N)
    else
        error(method)
    end
    println("* Solve")
    if method==:dense
        sol = LHS\vec(rhs)
    else
        # GMRES solve
        maxiter = 150
        tol = 1e-15
        sol, gmlog = IterativeSolvers.gmres(LHS, rhs; tol=tol, restart=maxiter, maxiter=maxiter, verbose=false, log=true)
        println( (gmlog.isconverged ? "  Converged" : "Did NOT converge"),
                 " in $(gmlog.iters) iterations, residual=",
                 gmlog.data[:resnorm][end])
    end
    density = reshape(sol, 2, N)

    # Estimate density resolution
    L = ModifiedStokesSolver.legendre_matrix(panelorder)
    density_resolution = 0.0;
    for i=1:numpanels
        xcoeff = L*density[1, (1:panelorder) + panelorder*(i-1)]
        ycoeff = L*density[2, (1:panelorder) + panelorder*(i-1)]
        xres = max(abs(xcoeff[end]), abs(xcoeff[end-1]))
        yres = max(abs(ycoeff[end]), abs(ycoeff[end-1]))
        res = max(xres, yres)
        density_resolution = max(density_resolution, res)
    end
    @show density_resolution
    
    ## VOLUME ERROR EVAL 
    dlp, dlpgrad1, dlpgrad2 = doublelayer_fast(dcurve, density, zt_interior, alpha, ifgrad=true)
    ref, refgrad1, refgrad2 = ubc_grad(zt_interior)
    @show relerr = maximum(abs.(ref.-dlp)) / maximum(abs.(ref))
    unorm = maximum(abs.(ref))
    gradnorm = max( maximum(abs.(refgrad1)), maximum(abs.(refgrad2)) )
    @show rmserr = sqrt( sum(vec((dlp.-ref).^2))/length(ref) ) / unorm

    @show graderr = max( maximum(abs.(refgrad1.-dlpgrad1)),
                         maximum(abs.(refgrad2.-dlpgrad2)) ) / gradnorm

    # plot
    err = abs.(dlp-ref)
    E = zero(X)
    E[interior] = max.(err[1,:],err[2,:]) / unorm
    figure(3)
    clf()
    pcolor(X, Y, log10.(E))
    clim(-15, -2)
    axis("equal")
    colorbar()


    
    return relerr, rmserr, graderr, density_resolution
end

ntests = length(numpanels_list)
max_errors = zeros(ntests)
rms_errors = zeros(ntests)
grad_errors = zeros(ntests)
density_resolution = zeros(ntests)
for i=1:ntests
    max_errors[i], rms_errors[i], grad_errors[i], density_resolution[i] = solve(numpanels_list[i])
end

close("all")

## Convergence plot
xquant = numpanels_list
figure(1, figsize=(4,3))
clf()
loglog(xquant, max_errors, ".-", label=L"||\tilde u - u||_\infty")
loglog(xquant, grad_errors, "*--", label=L"||\nabla\tilde u - \nabla u||_\infty")
#loglog(xquant, rms_errors, ".-", label="RMS err")
#loglog(xquant, density_resolution, ".--", label="Density")
# loglog(xquant, density_resolution[end]*(xquant/xquant[end]).^-(panelorder-2), "--k", label="$(panelorder-2)th order")
# loglog(xquant, density_resolution[end]*(xquant/xquant[end]).^-(panelorder-1), "--k", label="$(panelorder-1)th order")
# loglog(xquant, density_resolution[end]*(xquant/xquant[end]).^-panelorder, "--k", label="$(panelorder)th order")
mid = div(length(xquant),2)
loglog(xquant, (max_errors[mid]+grad_errors[mid])/2*(xquant/xquant[mid]).^-panelorder, "--k", label="$(panelorder)th order")
grid()
legend()
axis("auto")
ylim(1e-15, 1)
xticks(ticks, ticks)
minorticks_off()
xlabel("Number of panels")
ylabel("Relative error")
#title(L"\alpha="*string(alpha))
tight_layout(0)

## Domain plot
figure(2, figsize=(4,3))
clf()
t = linspace(0,2*pi,1000)
tau = curve.tau.(t)
plot(real(tau), imag(tau), label=L"\partial\Omega")
plot(src[1,:], src[2,:], "*", label="sources")
plot(zt_interior[1,:], zt_interior[2,:], ".", label="targets")
legend(loc="lower left")
axis("off")
axis("tight")
axis("equal")
tight_layout(0)

function write_fig(num, name)
    figure(num)
    savefig(name)
    println("Saved $name")
end

function save_conv_figs()
    write_fig(1, "../docs/paper/figs/convergence_homogeneous.pdf")
    write_fig(2, "../docs/paper/figs/convergence_homogeneous_geometry.pdf")
end


println("Done.\n")
println("Number of interior points: ", size(zt_interior,2),"\n")
println("Call save_conv_figs() to write plots to disk");
