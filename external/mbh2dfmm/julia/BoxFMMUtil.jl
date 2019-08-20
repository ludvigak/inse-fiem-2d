using PyPlot
using PyCall

function boxfmm2d_formquadvec(der1::Array{Float64},
                              der2::Array{Float64})
    quadvec = zeros(Float64,3)
    quadvec[1] = der1[1]*der2[1]
    quadvec[2] = der1[1]*der2[2] + der1[2]*der2[1]
    quadvec[3] = der1[2]*der2[2]
    return quadvec
end

function boxfmm2d_formoctvec(der1::Array{Float64},
                          der2::Array{Float64},
                          der3::Array{Float64})

    octvec = zeros(Float64,4)
    octvec[1] = der1[1]*der2[1]*der3[1] # xxx derivative
    octvec[2] = der1[1]*der2[1]*der3[2] +
        der1[1]*der2[2]*der3[1] +
        der1[2]*der2[1]*der3[1] # xxy derivative
    octvec[3] = der1[1]*der2[2]*der3[2] +
        der1[2]*der2[1]*der3[2] +
        der1[2]*der2[2]*der3[1] # xyy derivative
    octvec[4] = der1[2]*der2[2]*der3[2] # yyy derivative

    return octvec

end

function boxfmm2d_booltoint32(in::Bool)
    out = Int32[0]
    if (in)
        out[1] = 1
    end
    return out
end

function boxfmm2d_plot(tree::BoxTree2D)
    patch = pyimport("matplotlib.patches")

    
    nlev = tree.nlev
    nblevel = tree.nblevel
    istartlev = tree.istartlev
    blength = tree.blength
    zll = tree.zll
    icolbox = tree.icolbox
    irowbox = tree.irowbox
    iboxlev = tree.iboxlev

    fig = figure()
    ax = gca()
    axis([zll[1],zll[1]+blength,zll[2],zll[2]+blength])
    
    for i = 1:nlev+1
        h = blength*0.5^(i-1)
        for j = 1:nblevel[i]
            ind = istartlev[i]+j-1
            ibox = iboxlev[ind]
            icol = icolbox[ibox]
            irow = irowbox[ibox]

            rect = patch[:Rectangle]([zll[1]+h*(icol-1),zll[2]+h*(irow-1)],h,h,
                                     linewidth=2,facecolor="none",edgecolor="black")
            ax[:add_patch](rect)
        end
    end
            
    return
end
