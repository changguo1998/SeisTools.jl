module Seis
using DSP, Statistics, Printf, Dates

export Frame, WaveFrame, SACFrame, SEGYFrame, HEADER, readsac, readsachead, writesac, newsachead,
detrend, taper!, taper, filt!, filt

include("types.jl")
include("sac.jl")
include("segy.jl")
include("process.jl")

end
