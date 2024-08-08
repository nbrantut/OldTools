module OldTools

using NumericalIntegration: Trapezoidal, integrate, cumul_integrate

export linspace, logspace, linfit, trapz, cumtrapz, lininterp, cubicinterp, subsample, subsample_grow, subsample_decr, frequenciesdft, movingaverage, pow10latexstring, peaks

include("tools.jl")
include("interp.jl")
include("signalproc.jl")

end
