"""
	frequenciesdft(fs, N)

Compute the array of frequencies ordered as in output of fftw functions. see http://www.fftw.org/doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
"""
function frequenciesdft(fs, N)
    f = zeros(N)
    for i=1:N
        if i<=(N/2+1)
            f[i] = i-1
        else
            f[i] = i-1-N
        end
    end
    f *= 2π*fs/N

    return f
end

"""
	movingaverage(X, N)

Compute moving average for vector X with N points. 
"""
function movingaverage(X,N)
    backwin = div(N,2) 
    fwdwin = isodd(N) ? div(N,2) : div(N,2) - 1
    
    Nx = length(X)
    Y = zeros(eltype(X),size(X))
    
    # first points
    for n in 1:backwin+1
        for i in 1:n+fwdwin
            Y[n] += X[i]
        end
        Y[n] /= (n+fwdwin)
    end
    
    @inbounds for n in backwin+2:Nx-fwdwin
        i0 = n - backwin - 1
        i1 = n + fwdwin

        Y[n] = Y[n-1] + (X[i1] - X[i0])/N
    end

    # last points
    for n in Nx-fwdwin+1:Nx
        for i in n-backwin:Nx
            Y[n] += X[i]
        end
        Y[n] /= (Nx-n+backwin+1)
    end
    
    return Y
end
