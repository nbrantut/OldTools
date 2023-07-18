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

Compute moving average for vector X with N points. This is modified from https://stackoverflow.com/questions/59562325/moving-average-in-julia
"""
function movingaverage(X::Vector, N::Int)
    backwin = div(N,2) 
    fwdwin = isodd(N) ? div(N,2) : div(N,2) - 1
    Nx = length(X)
    Y = similar(X)
    @inbounds for n in eachindex(X)
        i0 = max(1,n - backwin)
        i1 = min(Nx,n + fwdwin)
        Y[n] = zero(eltype(Y))
        for i in i0:i1
            Y[n] += X[i]
        end
        Y[n] /= (i1-i0+1)
    end
    return Y
end
