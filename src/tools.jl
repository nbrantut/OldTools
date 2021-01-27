"""
    logspace(n1, n2, N=50)

Generate a logarithmically spaced array ranging from 10^n1 to 10^n2 with N points. Same as old Matlab function.
"""
function logspace(n1, n2, N=50)
    return exp10.(range(n1; length=N, stop=n2))
end

"""
    linspace(n1, n2, N=50)

Generate a linearly spaced Range from n1 to n2 with N points. Same as old Matlab function. This is just another name for `range(n1; length=N, stop=n2)` but useful for Matlab-heads like me.
"""
function linspace(n1, n2, N=50)
    return collect(range(n1; length=N, stop=n2))
end

"""
    linfit(x,y)

The most basic least-square regression algorithm. Take two arrays `x` and `y` and fit a straight line through it as `y = ax+b`. Return the tuple `(a,b)`.
"""
function linfit(x, y)
    if length(x)!=length(y)
        error("x and y should be of the same size")
    end
    
    A = [x  ones(size(x,1))]

    coefs = A\y

    return coefs[1], coefs[2]
end

"""
    trapz(x,y)

Compute numerical integral of y(x) with trapezoidal rule.
"""
trapz(x, y) = integrate(x,y,Trapezoidal())

"""
    cumtrapz(x,y)

Compute cumulative numerical integral of y(x) with trapezoidal rule.
"""
cumtrapz(x,y) = cumul_integrate(x,y,Trapezoidal())

"""
	cleancurve(x0,y0,step)

Generate a pair of vectors (x,y) such that the values of y are separated by the given step.
"""
function cleancurve(x0,y0,step)
    x=eltype(x0)[]
    y=eltype(y0)[]

    push!(x, x0[1])
    push!(y, y0[1])

    for k in eachindex(y0)        
        if isfinite(y0[k]) && abs(y0[k]-y[end])>=step
            push!(x, x0[k])
            push!(y, y0[k])
        end
    end
    return (x,y)
end

"""
	cleancurvei(x0,y0,step)

Generate a pair of vectors (x,y) such that the values of y are separated by the given step. Forces the sequence of y's to be monotonically increasing.
"""
function cleancurvei(x0,y0,step)
    x=eltype(x0)[]
    y=eltype(y0)[]

    push!(x, x0[1])
    push!(y, y0[1])

    for k in eachindex(y0)        
        if isfinite(y0[k]) && (y0[k]-y[end])>=step
            push!(x, x0[k])
            push!(y, y0[k])
        end
    end
    return (x,y)
end

"""
	pow10latexstring(x,n=3;symbol="\\times")

    Make a valid LaTeX string to render a number `x` into something like a×10ᵇ.
"""
function pow10latexstring(x,n=3;symbol="\\times")
    sn = sign(x)
    x *= sn
    exponent = round(Int, log10(x))
    pf = sn*x /( 10.0^exponent)
    rpf = round(pf, sigdigits=n)
    return "\$ $(rpf) $(symbol) 10^{$(exponent)}\$"
end


"""
	peaks(n)

Return arrays x, y, and z of size, n, n, and n-by-n, respectively, where (x,y) is a regular grid and z is a nicely peaked function. Similar to Matlab's peaks function.
"""
function peaks(n=50)
    x = collect(range(-3,stop=3,length=n))
    y = collect(range(-3,stop=3,length=n))
    z = zeros(n,n)
    for j=1:n, i=1:n
        z[i,j] =  3*(1-x[i])^2*exp(-(x[i]^2) - (y[j]+1)^2) -
            10*(x[i]/5 - x[i]^3 - y[j]^5)*exp(-x[i]^2-y[j]^2) -
            1/3*exp(-(x[i]+1)^2 - y[j]^2)
    end

    return x,y,z
end
