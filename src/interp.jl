function lininterp(x::AbstractVector,y::AbstractVector,xi::AbstractVector)
    yi = zeros(eltype(y), length(xi))

    # first check if x is sorted and there are no repeated values
    if ~allunique(x)
        error("Values in x should be distinct")
    end
    
    if ~issorted(x)
        I = sortperm(x)
        x = x[I] #this apparently creates a copy of x, so output not modifed
        y = y[I]
    end

    # sort xi
    if ~issorted(xi)
        J = sortperm(xi)
        xi = xi[J]
    else
        J = 1:length(xi)
    end
    
    
    N = length(x)

    for (k,xk) in enumerate(xi)
        m = searchsortedlast(x, xk)
        n = min(max(1, m), N-1)
        yi[J[k]] = y[n] + (xk - x[n])*(y[n+1] - y[n])/(x[n+1]- x[n])
    end

    return yi

end

function lininterp(x::AbstractVector,y::AbstractVector,xi::Number)

    # first check if x is sorted and there are no repeated values
    if ~allunique(x)
        error("Values in x should be distinct")
    end
    
    if ~issorted(x)
        I = sortperm(x)
        x = x[I] #this apparently creates a copy of x, so output not modifed
        y = y[I]
    end

    N = length(x)

    m = searchsortedlast(x, xi)
    n = min(max(1, m), N-1)
    yi = y[n] + (xi - x[n])*(y[n+1] - y[n])/(x[n+1]- x[n])

    return yi

end

function cubicinterp(x::AbstractVector,y::AbstractVector,xi::AbstractVector)
    yi = zeros(eltype(y), length(xi))

    # first check if x is sorted and there are no repeated values
    if ~allunique(x)
        error("Values in x should be distinct")
    end

    if length(x)<3
        error("Input must be of size >=3")
    end
    
    if ~issorted(x)
        I = sortperm(x)
        x = x[I] #this apparently creates a copy of x, so output not modifed
        y = y[I]
    end

    # sort xi
    if ~issorted(xi)
        J = sortperm(xi)
        xi = xi[J]
    else
        J = 1:length(xi)
    end

    if (xi[1]<x[1]) || (xi[end]>x[end])
        error("Cannot extrapolate")
    end
    
    N = length(x)

    # initialise m, D, h for first set of three points
    k = 1
    h = x[2]-x[1]
    hkp = x[3]-x[2]
    Dk = (y[2]-y[1])/h
    Dkp = (y[3]-y[2])/hkp
    if abs(Dk)<=eps()
        mk = 0.0
        mkp = 0.0
    else
        mk = Dk
        mkp = Dk*Dkp<0.0 ? 0.0 : 0.5*(Dk+Dkp)
        r2 = (mk^2+mkp^2)/(Dk^2)
        if r2>9.0
            mk *= 3/sqrt(r2)
            mkp *= 3/sqrt(r2)
        end
    end
    
    for (i,s) in enumerate(xi)
        if s>x[k+1]
            # then search for next index
            k = searchsortedlast(x, xi[i]; lt=<=)

            # and recompute what is needed
            Dkm = (y[k]-y[k-1])/(x[k]-x[k-1])
            
            h = x[k+1]-x[k]
            Dk = (y[k+1]-y[k])/h

            Dkp = k<N-1 ? (y[k+2]-y[k+1])/(x[k+2]-x[k+1]) : Dk
            
            mk = Dkm*Dk<0.0 ? 0.0 : 0.5*(Dkm+Dk)
            mkp = Dk*Dkp<0.0 ? 0.0 : 0.5*(Dk+Dkp)
            
            if abs(Dk)<=eps()
                mk = 0.0
                mkp = 0.0
            else
                r2 = (mk^2+mkp^2)/(Dk^2)
                if r2>9.0
                    mk *= 3/sqrt(r2)
                    mkp *= 3/sqrt(r2)
                end
            end
        end

        # compute nterpolation
        t = (s-x[k])/h
        t2 = t^2
        t3 = t^3
        h00 = 2*t3-3*t2+1
        h10 = t3-2*t2+t
        h01 = -2*t3+3*t2
        h11 = t3-t2
        yi[J[i]] = y[k]*h00 + h*mk*h10 + y[k+1]*h01 + h*mkp*h11
    end

    return yi
end

