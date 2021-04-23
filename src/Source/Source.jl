module Source
using LinearAlgebra, Plots
pyplot()
include("types.jl")
export Mechanism, DoubleCouple, MomentTensor

include("DC.jl")
export normal2DC, normalvector

include("Beachball.jl")
export beachball
end