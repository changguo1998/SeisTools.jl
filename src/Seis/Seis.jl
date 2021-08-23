module Seis
using DSP, Statistics, Printf, Dates

export Frame, WaveFrame, SACFrame, SEGYFrame, HEADER, readsac, readsachead, writesac, newsachead,
detrend, taper!, taper, filt!, filt

include("types.jl")
include("process.jl")
include("sac.jl")

end
