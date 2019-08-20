using Compat


function advance_trace(s1, s2, U1interp, U2interp, dt)
    u1 = U1interp.(s1,s2)
    u2 = U2interp.(s1,s2)
    @. s1 += dt*u1
    @. s2 += dt*u2    
end
    
function streakplot(;X=[], Y=[], interior=[], dataname=nothing, ilist=[],
                    startx=[],
                    starty=[])
    x = X[1:end,1]
    y = Y[1,1:end]
    traces = length(startx)
    s1 = Array{Array{Float64}}(undef, traces)
    s2 = Array{Array{Float64}}(undef, traces)    
    for k=1:traces
        s1[k] = [startx[k]]
        s2[k] = [starty[k]]       
    end
    # Iterate
    t = load(dataname(ilist[1]))["ti"] # Get initial time
    stepnum = 0
    for i=ilist
        stepnum += 1
        datafile = dataname(i)
        info("Loading ", datafile)
        U1,U2 = [],[]
        dt = 0.0
        @time begin
            U1 = h5read(datafile, "U1")
            U2 = h5read(datafile, "U2")
            dt = h5read(datafile, "ti") - t
        end
        # Set free stream speed to get traces coming in and out
        @. U1[ (U1==0.0) & (U2==0.0) ] = 1.0        
        U1interp = LinearInterpolation((x,y), U1, extrapolation_bc=0)
        U2interp = LinearInterpolation((x,y), U2, extrapolation_bc=0)        
        t += dt
        for k=1:traces
            advance_trace(s1[k], s2[k], U1interp, U2interp, dt)
        end
        # Add points, fill in gaps and trim
        if mod(stepnum, 10) == 0
            for k=1:traces
                # Add points to beginning
                unshift!(s1[k], startx[k])
                unshift!(s2[k], starty[k])                
                # Interpolate into line
                kx, ky = s1[k], s2[k]
                m = 1
                while m <= length(kx)-1
                    dx = kx[m+1]-kx[m]
                    dy = ky[m+1]-ky[m]
                    if sqrt(dx^2+dy^2) > 0.1
                        newx = kx[m]+dx/2
                        newy = ky[m]+dy/2
                        splice!(kx, m+1:m, newx)
                        splice!(ky, m+1:m, newy)                            
                    else
                        m += 1                            
                    end
                end
                xmax = maximum(x)-0.1
                inside = @. kx < xmax
                s1[k] = kx[inside]
                s2[k] = ky[inside]
            end
        end
    end
    # Plot
    begin
        marker = zeros(size(X))
        marker[interior] = 1.0
        markerInterp = LinearInterpolation((x,y), marker, extrapolation_bc=0)
        for k=1:traces
            kx, ky = s1[k], s2[k]
            #inside = @. kx < 5
            inside = markerInterp.(kx, ky) .> 0
            plot(kx[inside], ky[inside], "-k")
        end
    end    
end
