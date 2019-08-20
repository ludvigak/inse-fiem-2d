using Compat

const LIBMBHFMM2D = string(Base.Filesystem.dirname(Base.source_path()), "/../bin/libmbhfmm2d")

function modbh_qn(x::Array{Float64,1},
                   lambda::Float64,rscale::Float64,n::Int)

    ifders = zero(Int32)
    ders = zeros(Float64,1)
    n = convert(Int32,n)

    n2 = length(x)

    qtemp = zeros(Float64,n+2)
    qn = zeros(Float64,n+1,n2)

    kvec = zeros(Float64,n+5)
    
    for i = 1:n2
        ccall( (:diffslogbk_fast_,LIBMBHFMM2D), Cvoid,
               (Ref{Float64},Ref{Float64},Ref{Float64},
                Ref{Float64},Ref{Int32},Ref{Float64},
                Ref{Float64},Ref{Int32}),
               x[i],lambda,rscale,qtemp,ifders,ders,
               kvec,n)
        qn[:,i] = qtemp[1:n+1]
    end
    
    return qn
end

function besselq(n::Int,x::Float64,lambda::Float64;
                 rscale::Float64=1.0)
    xvec = zeros(Float64,1)
    xvec[1] = x
    qn = modbh_qn(xvec,lambda,rscale,n)

    return qn[n+1,1]
end

function modbh_qn_test()
    xtest1 = [0.1]
    lamtest1 = 1.0
    q0test1 = 0.1244839317079709285005145657

    rscale = 1.0
    n = 1
    
    qn = modbh_qn(xtest1,lamtest1,rscale,n)

    err = abs(q0test1-qn[1,1])
    relerr = err/abs(q0test1)

    return err, relerr
end
    
function modbh_pn(x::Array{Float64,1},
                   lambda::Float64,rscale::Float64,n::Int)

    ifders = zero(Int32)
    ders = zeros(Float64,1)
    n = convert(Int32,n)

    n2 = length(x)

    ptemp = zeros(Float64,n+2)
    pn = zeros(Float64,n+1,n2)

    ivec = zeros(Float64,n+5)
    
    for i = 1:n2
        ccall( (:diffszkik_fast_,LIBMBHFMM2D), Cvoid,
               (Ref{Float64},Ref{Float64},Ref{Float64},
                Ref{Float64},Ref{Int32},Ref{Float64},
                Ref{Float64},Ref{Int32}),
               x[i],lambda,rscale,ptemp,ifders,ders,
               ivec,n)
        pn[:,i] = ptemp[1:n+1]
    end
    
    return pn
end

function modbh_pn_test()
    xtest1 = [0.1]
    lamtest1 = 1.0
    p0test1 = 0.00250156293409560140021055

    rscale = 1.0
    n = 1
    
    pn = modbh_pn(xtest1,lamtest1,rscale,n)

    err = abs(p0test1-pn[1,1])
    relerr = err/abs(p0test1)

    return err, relerr
end
    
function besselp(n::Int,x::Float64,lambda::Float64;
                 rscale::Float64=1.0)
    xvec = zeros(Float64,1)
    xvec[1] = x
    pn = modbh_pn(xvec,lambda,rscale,n)

    return pn[n+1,1]
end

function besselj0m1(x::Float64)
    # evaluate J_0(x) - 1 stably
    
    xh2 = (0.5*x)^2
    if (xh2 <= 1)
        tot = 0.0
        xterm = 1.0
        for i = 1:12
            xterm = -xterm*xh2/(i*i)
            tot = tot + xterm
        end
    else
        tot = besselj(0,x)-1
    end

    return tot
end

function besselj0m1_test()
    testin = 0.1
    testval = -0.002498437933959967718713101525
    funval = besselj0m1(testin)
    err = abs(testval-funval)
    relerr = err/abs(testval)

    return err, relerr
end
