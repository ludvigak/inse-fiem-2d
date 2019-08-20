function plot_results(;paperplots=false)
    if paperplots
        close(3)
        figure(3, figsize=(4,3))
    end

    figure(3)
    clf()

    #N = 1./hgrid
    N = Ngrid_list

    
    region = (l2_rel_err .> 1e-13) .& (l2_rel_err .< 1e-5)
    p = -10
    c = match_convergence(N[region], l2_rel_err[region], p)
    loglog(N, c*N.^(p*1.0), "-k", label=L"\mathcal{O}(N^{-10})", zorder=-1)

    loglog(N, grad_max_rel_err, "*-", label=L"||\nabla \tilde u-\nabla u||_{\infty}")
    loglog(N, grad_l2_rel_err, "+-", label=L"||\nabla \tilde u-\nabla u||_{2}")    
    loglog(N, max_rel_err, ".-", label=L"||\tilde u-u||_{\infty}")
    loglog(N, l2_rel_err, "x-", label=L"||\tilde u-u||_{2}")


    region = (grad_l2_rel_err .> 1e-10) .& (grad_l2_rel_err .< 1e-4)
    p = -10
    c = match_convergence(N[region], grad_l2_rel_err[region], p)
    loglog(N, c*N.^(p*1.0), "-k", zorder=-1)
    
    legend()        
    xlabel("N")
    grid("on")
    minorticks_off()
    ticks = [50,100,200,400,800]
    xticks(ticks, round.(Int, ticks))
    ylim(1e-15, 1)
    ylabel("Relative error")

    tight_layout(0)
end
