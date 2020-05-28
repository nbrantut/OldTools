module OldTools

using NumericalIntegration: integrate, cumul_integrate

export linspace, logspace, linfit, trapz, cumtrapz, lininterp, cubicinterp, cleancurve, cleancurvei, frequenciesdft, movingaverage, pow10latexstring

include("tools.jl")
include("interp.jl")
include("signalproc.jl")

end
