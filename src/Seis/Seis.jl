module Seis
using DSP, Statistics, Printf

include("types.jl")
export Frame, WaveFrame, SACFrame, SEGYFrame, HEADER

include("sac.jl")
export readsac, readsachead, writesac, newsachead

include("process.jl")
export detrend, taper!, taper, filt!, filt
end