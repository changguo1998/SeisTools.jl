module VelocityModel

include("macros.jl")

function _read_f64_vector(io::IO)
    n = read!(io, Int64)
    v = zeros(Float64, n)
    read!(io, v)
    return v
end

function _write_f64_vector!(io::IO, v::Vector{<:Real})
    b = zeros(Float64, length(v))
    b .= Float64.(v)
    n = Int64(length(b))
    write(io, n)
    write(io, b)
    return nothing
end

function _read_f64_array(io::IO)
    dm = read(io, Int64)
    sz = zeros(Int64, dm)
    read!(io, sz)
    arr = zeros(Float64, sz)
    read!(io, arr)
    return arr
end

function _write_f64_array(io::IO, arr::Array{<:Real})
    dm = ndims(arr)
    sz = zeros(Int64, dm)
    sz .= size(arr)
    buf = zeros(Float64, sz)
    buf .= Float64.(arr)
    write(io, Int64(dm))
    write(io, sz)
    write(io, buf)
    return nothing
end


"""
`AbstractVelocityModel` including velocity model types as below:

- FlatLayer
- SphereLayer
- PsudoRegular3D
- PsudoIrregular3D
- Regular3D
- Irregular3D
"""
abstract type AbstractVelocityModel <: Any end

struct FlatLayer <: AbstractVelocityModel
    d::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    ρ::Vector{Float64}
end

function FlatLayer(d::Vector{<:Real}=Float64[],
    α::Vector{<:Real}=Float64[],
    β::Vector{<:Real}=Float64[],
    ρ::Vector{<:Real}=Float64[])
    return FlatLayer(Float64.(d), Float64.(α), Float64.(β), Float64.(ρ))
end

function FlatLayer(vel::Matrix{<:Real}; α::Real=5.0, β::Real=3.0, ρ::Real=2.7)
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return FlatLayer(vel[:,1], fill(α, nlayer), fill(β, nlayer), fill(ρ, nlayer))
    elseif size(vel, 2) == 3
        return FlatLayer(vel[:,1], vel[:, 2], vel[:, 3], 1.743.*(vel[:, 2].^0.25))
    elseif size(vel, 2) == 4
        return FlatLayer(vel[:,1], vel[:, 2], vel[:, 3], vel[:, 4])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

struct SphereLayer <: AbstractVelocityModel
    r::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    ρ::Vector{Float64}
end

function SphereLayer(name::String)
end

struct PsudoRegular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Array{Float64}
    α::Array{Float64}
    β::Array{Float64}
    ρ::Array{Float64}
end

function PsudoRegular3D(fp_or_mn::String)
    if isfile(fp_or_mn)
        io = open(fp_or_mn, "r")
        x = _read_f64_vector(io)
        y = _read_f64_vector(io)
        d = _read_f64_array(io)
        α = _read_f64_array(io)
        β = _read_f64_array(io)
        ρ = _read_f64_array(io)
        close(io)
        return PsudoRegular3D(x, y, d, α, β, ρ)
    end
    available_velocity_models = String[]
    for md in ("CRUST1.0")
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "data.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_mn ∈ available_velocity_models
        return PsudoRegular3D(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_mn, "data.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

struct PsudoIrregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Matrix{Float64}
    α::Matrix{Float64}
    β::Matrix{Float64}
    ρ::Matrix{Float64}
end

struct Regular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    α::Array{Float64}
    β::Array{Float64}
    ρ::Array{Float64}
end

struct Irregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    α::Matrix{Float64}
    β::Matrix{Float64}
    ρ::Matrix{Float64}
end

end