
function spindown_convergence_plot(dtlist, l2err, maxerr, sdc_order)
    orderstring = latexstring("\\mathcal{O}(\\Delta t^$sdc_order)")
    h = loglog(dtlist, maxerr, ".-", label=LaTeXString("\$||\\tilde u^{($sdc_order)} - u||_\\infty\$"))
    color = h[1][:get_color]()
    loglog(dtlist, l2err, "*-", label=LaTeXString("\$||\\tilde u^{($sdc_order)} - u||_2\$"), color=color)    
    matchidx = 4
    C = ( maxerr[matchidx] + l2err[matchidx] )/2
    @show dtlist_ext = vcat(dtlist[1]*2, dtlist, dtlist[end]/2)
    loglog(dtlist_ext, C*(dtlist_ext/dtlist[matchidx]).^sdc_order, "--", label=orderstring, color=color)
    #axis("equal")
    grid("on")    
    #xlim(left=3e-5, right=1e-2)
    axis("tight")
    ylim([1e-13, 1e-2])    
    xlabel(L"Time step $\Delta t$")
    ylabel("Relative error")
    legend(bbox_to_anchor=(1.0,1.03), ncol=2)
    tight_layout(0.1)
end


close("all")
figure(figsize=(9,3))
spindown_convergence_plot(dtlist1, l2err1, maxerr1, 1)
spindown_convergence_plot(dtlist2, l2err2, maxerr2, 2)
spindown_convergence_plot(dtlist3, l2err3, maxerr3, 3)
spindown_convergence_plot(dtlist4, l2err4, maxerr4, 4)

savefig("../docs/paper/figs/convergence_time.pdf")
