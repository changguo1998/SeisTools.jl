module VelocityModel

include("basic.jl")

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
    Qk::Vector{Float64}
    Qm::Vector{Float64}
end

function FlatLayer(d::Vector{<:Real} = Float64[],
                   α::Vector{<:Real} = Float64[],
                   β::Vector{<:Real} = Float64[],
                   ρ::Vector{<:Real} = Float64[],
                   Qk::Vector{<:Real} = Float64[],
                   Qm::Vector{<:Real} = Float64[])
    return FlatLayer(Float64.(d), Float64.(α), Float64.(β), Float64.(ρ), Float64.(Qk), Float64.(Qm))
end

function FlatLayer(vel::Matrix{<:Real}; α::Real = 5.0, β::Real = 3.0, ρ::Real = 2.7, Qk::Real = 300.0, Qm::Real = 150.0)
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return FlatLayer(vel[:, 1], fill(α, nlayer), fill(β, nlayer), fill(ρ, nlayer), fill(Qk, nlayer),
                         fill(Qm, nlayer))
    elseif size(vel, 2) == 3
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], 1.743 .* (vel[:, 2] .^ 0.25), fill(Qk, nlayer),
                         fill(Qm, nlayer))
    elseif size(vel, 2) == 4
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], fill(Qk, nlayer), fill(Qm, nlayer))
    elseif size(vel, 2) == 5
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], fill(Qm, nlayer))
    elseif size(vel, 2) == 6
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], vel[:, 6])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

function FlatLayer(fp_or_nm::String)
    if isfile(fp_or_nm)
        io = open(fp_or_nm, "r")
        d = _read_f64_vector(io)
        α = _read_f64_vector(io)
        β = _read_f64_vector(io)
        ρ = _read_f64_vector(io)
        Qk = _read_f64_vector(io)
        Qm = _read_f64_vector(io)
        close(io)
        return FlatLayer(d, α, β, ρ, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("AK135", "IASP91", "PREM")
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "flat.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return FlatLayer(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "flat.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

struct SphereLayer <: AbstractVelocityModel
    r::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    ρ::Vector{Float64}
    Qk::Vector{Float64}
    Qm::Vector{Float64}
end

function SphereLayer(r::Vector{<:Real} = Float64[],
                     α::Vector{<:Real} = Float64[],
                     β::Vector{<:Real} = Float64[],
                     ρ::Vector{<:Real} = Float64[],
                     Qk::Vector{<:Real} = Float64[],
                     Qm::Vector{<:Real} = Float64[])
    return SphereLayer(Float64.(r), Float64.(α), Float64.(β), Float64.(ρ), Float64.(Qk), Float64.(Qm))
end

function SphereLayer(vel::Matrix{<:Real}; α::Real = 5.0, β::Real = 3.0, ρ::Real = 2.7, Qk::Real = 300.0,
                     Qm::Real = 150.0)
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return SphereLayer(vel[:, 1], fill(α, nlayer), fill(β, nlayer), fill(ρ, nlayer), fill(Qk, nlayer),
                           fill(Qm, nlayer))
    elseif size(vel, 2) == 3
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], 1.743 .* (vel[:, 2] .^ 0.25), fill(Qk, nlayer),
                           fill(Qm, nlayer))
    elseif size(vel, 2) == 4
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], fill(Qk, nlayer), fill(Qm, nlayer))
    elseif size(vel, 2) == 5
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], fill(Qm, nlayer))
    elseif size(vel, 2) == 6
        return SphereLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], vel[:, 6])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

function SphereLayer(fp_or_nm::String)
    if isfile(fp_or_nm)
        io = open(fp_or_nm, "r")
        r = _read_f64_vector(io)
        α = _read_f64_vector(io)
        β = _read_f64_vector(io)
        ρ = _read_f64_vector(io)
        Qk = _read_f64_vector(io)
        Qm = _read_f64_vector(io)
        close(io)
        return SphereLayer(r, α, β, ρ, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("AK135", "IASP91", "PREM")
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "sphere.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return SphereLayer(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "sphere.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

struct PsudoRegular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Array{Float64}
    α::Array{Float64}
    β::Array{Float64}
    ρ::Array{Float64}
end

function PsudoRegular3D(fp_or_nm::String)
    if isfile(fp_or_nm)
        io = open(fp_or_nm, "r")
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
    if fp_or_nm ∈ available_velocity_models
        return PsudoRegular3D(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "data.bin"))
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