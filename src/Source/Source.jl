module Source
using LinearAlgebra
include("types.jl")
export Mechanism, DoubleCouple, MomentTensor, normal2DC, normalvector, beachball

include("DC.jl")
include("Beachball.jl")

end
