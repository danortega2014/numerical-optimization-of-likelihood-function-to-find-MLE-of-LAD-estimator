function bisection(fun,a,b,tol)
    sa = sign(fun(a))
    sb = sign(fun(b))
    if sa == sb
        error("Interval is not valid")
    end

    while abs(b-a)>tol
        c = (a+b)/2
        sc = sign(fun(c))
        if sc == 0
            a,b = c,c
        elseif sa==sc
            a,b = c,b
        else
            a,b = a,c
        end
    end
    return (a+b)/2
end



using Random, Distributions


x= rand(Normal(),100)
y= rand(Normal(),100)

function loglike(Beta)
    llbeta  = -length(x)*log(2)-sum(abs.(y .-x*Beta))
    return -llbeta
end


function derivat(Beta)
    h = .00001
    g2 = (loglike(Beta + h) - loglike(Beta) ) / h
    return g2
end


mle = bisection(derivat,.0001, 100, 1e-15)

using Optim
opt = optimize(loglike, .0001, 100, Brent())
