function plot_results(;paperplots=false)
    if paperplots
        close("all")
        figure(3, figsize=(4,3))
    end

    figure(3)
    clf()

    #N = 1./hgrid
    N = alpha_list
    
    loglog(N, grad_max_rel_err, "*-", label=L"||\nabla \tilde u-\nabla u||_{\infty}")
    loglog(N, grad_l2_rel_err, "+-", label=L"||\nabla \tilde u-\nabla u||_{2}")    
    loglog(N, max_rel_err, ".-", label=L"||\tilde u-u||_{\infty}")
    loglog(N, l2_rel_err, "x-", label=L"||\tilde u-u||_{2}")
    
    legend()        
    xlabel(L"\alpha")
    grid("on")   
    ylim(1e-15, 1)
    ylabel("Relative error")

    tight_layout(0)
end
