module SeisTools

export Seis, Source, Frame, WaveFrame, SACFrame, SEGYFrame, HEADER, readsac, readsachead, writesac, newsachead, detrend,
       taper!, taper, filt!, filt, Mechanism, DoubleCouple, MomentTensor, normal2DC, normalvector, beachball

include("Seis/Seis.jl")
include("Source/Source.jl")

end # module
