# Show cancellation errors in stresslet when evaluated directly

push!(LOAD_PATH, string(pwd(),"/src"))
using PyPlot
using ModifiedStokesSolver

for ifprime = [false, true]
    N = 100
    z = linspace(eps(), 3, N)
    E = zeros(N, 6)
    for i=1:N
        Tp = ModifiedStokesSolver.Ti_split_pow(z[i], ifprime)    
        Td = ModifiedStokesSolver.Ti_split_direct(z[i], ifprime)
        for j=1:6
            E[i, j] = abs(Tp[j]-Td[j]) / maximum(abs.(Tp))
        end
    end

    figure(ifprime ? 1 : 2)    
    clf()
    semilogy(z, E[:,1], label="T1S")
    semilogy(z, E[:,2], label="T1L")
    semilogy(z, E[:,3], label="T2S")
    semilogy(z, E[:,4], label="T2L")
    semilogy(z, E[:,5], label="T3S")
    semilogy(z, E[:,6], label="T3L")
    ylim(1e-16, 1)

    if ifprime
        stem(ModifiedStokesSolver.POWERSERIES_ZMAX_PRIME*[1.0, 1.0], [1e-16, 1], label="cutoff")
    else    
        stem(ModifiedStokesSolver.POWERSERIES_ZMAX*[1.0, 1.0], [1e-16, 1], label="cutoff")
    end
    title("ifprime=$ifprime")

    legend()
    grid("on")
end
show()
