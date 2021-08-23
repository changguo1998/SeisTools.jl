module Seis
using DSP, Statistics, Printf, Dates

export Frame, WaveFrame, SACFrame, SEGYFrame, HEADER, readsac, readsachead, writesac, newsachead,
detrend, taper!, taper, filt!, filt

include("types.jl")
<<<<<<< HEAD
=======
include("sac.jl")
include("segy.jl")
>>>>>>> e688cec2c3f59bba732f955f5a0825338594c82d
include("process.jl")
include("sac.jl")

end
