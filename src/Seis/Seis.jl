module Seis
export SAC, detrend, taper!, taper, filt!, filt

include("types.jl")
include("sac.jl")
include("segy.jl")
include("process.jl")

end
