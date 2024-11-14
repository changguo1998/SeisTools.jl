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
    vp::Vector{Float64}
    vs::Vector{Float64}
    dens::Vector{Float64}
    Qk::Vector{Float64}
    Qm::Vector{Float64}

    function FlatLayer(d::Vector{Float64},
                       vp::Vector{Float64},
                       vs::Vector{Float64},
                       dens::Vector{Float64},
                       Qk::Vector{Float64},
                       Qm::Vector{Float64})
        flag = map(eachindex(d)) do i
            if i < length(d)
                if d[i] == d[i+1]
                    return false
                end
            end
            local tf = true
            if i == 1
                return tf
            end
            if !isempty(vp)
                tf &= vp[i-1] == vp[i]
            end
            if !isempty(vs)
                tf &= vs[i-1] == vs[i]
            end
            if !isempty(dens)
                tf &= dens[i-1] == dens[i]
            end
            if !isempty(Qk)
                tf &= Qk[i-1] == Qk[i]
            end
            if !isempty(Qm)
                tf &= Qm[i-1] == Qm[i]
            end
            return !tf
        end
        _d = isempty(d) ? d : d[flag]
        _vp = isempty(vp) ? vp : vp[flag]
        _vs = isempty(vs) ? vs : vs[flag]
        _dens = isempty(dens) ? dens : dens[flag]
        _Qk = isempty(Qk) ? Qk : Qk[flag]
        _Qm = isempty(Qm) ? Qm : Qm[flag]
        return new(_d, _vp, _vs, _dens, _Qk, _Qm)
    end
end

function FlatLayer(d::Vector{<:Real} = Float64[],
                   vp::Vector{<:Real} = Float64[],
                   vs::Vector{<:Real} = Float64[],
                   dens::Vector{<:Real} = Float64[],
                   Qk::Vector{<:Real} = Float64[],
                   Qm::Vector{<:Real} = Float64[])
    return FlatLayer(Float64.(_d), Float64.(_vp), Float64.(_vs), Float64.(_dens), Float64.(_Qk), Float64.(_Qm))
end

function FlatLayer(vel::Matrix{<:Real})
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return FlatLayer(vel[:, 1])
    elseif size(vel, 2) == 3
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], 1.743 .* (vel[:, 2] .^ 0.25))
    elseif size(vel, 2) == 4
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4])
    elseif size(vel, 2) == 6
        return FlatLayer(vel[:, 1], vel[:, 2], vel[:, 3], vel[:, 4], vel[:, 5], vel[:, 6])
    else
        error("Invalid shape of velocity matrix. see document for more information.")
    end
end

function FlatLayer(fp_or_nm::String)
    if isfile(fp_or_nm)
        io = open(fp_or_nm, "r")
        n = read(io, Int64)
        if n > 0
            d = _read_f64_vector(io)
        else
            d = Float64[]
        end
        if n > 1
            vp = _read_f64_vector(io)
        else
            vp = Float64[]
        end
        if n > 2
            vs = _read_f64_vector(io)
        else
            vs = Float64[]
        end
        if n > 3
            dens = _read_f64_vector(io)
        else
            dens = Float64[]
        end
        if n > 4
            Qk = _read_f64_vector(io)
        else
            Qk = Float64[]
        end
        if n > 5
            Qm = _read_f64_vector(io)
        else
            Qm = Float64[]
        end
        close(io)
        return FlatLayer(d, vp, vs, dens, Qk, Qm)
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

function Matrix(m::FlatLayer)
    lfieldnames = [:d, :vp, :vs, :dens, :Qk, :Qm]
    lf = filter(_f -> !isempty(getfield(m, _f)), lfieldnames)
    return reduce(hcat, map(_f -> getfield(m, _f), lf))
end

function _slice_index_or_empty(v::Vector, x::Vector{<:Integer})
    if isempty(v)
        return eltype(v)[]
    else
        return v[x]
    end
end

function update(m::FlatLayer; d::Vector{<:Real} = Float64[],
                vp::Vector{<:Real} = Float64[],
                vs::Vector{<:Real} = Float64[],
                dens::Vector{<:Real} = Float64[],
                Qk::Vector{<:Real} = Float64[],
                Qm::Vector{<:Real} = Float64[])
    _d = isempty(d) ? m.d : d
    _vp = isempty(vp) ? m.vp : vp
    _vs = isempty(vs) ? m.vs : vs
    _dens = isempty(dens) ? m.dens : dens
    _Qk = isempty(Qk) ? m.Qk : Qk
    _Qm = isempty(Qm) ? m.Qm : Qm
    return FlatLayer(_d, _vp, _vs, _dens, _Qk, _Qm)
end

function cut(m::FlatLayer, idx::AbstractVector{<:Integer})
    _d = isempty(m.d) ? m.d : m.d[idx]
    _vp = isempty(m.vp) ? m.vp : m.vp[idx]
    _vs = isempty(m.vs) ? m.vs : m.vs[idx]
    _dens = isempty(m.dens) ? m.dens : m.dens[idx]
    _Qk = isempty(m.Qk) ? m.Qk : m.Qk[idx]
    _Qm = isempty(m.Qm) ? m.Qm : m.Qm[idx]
    return FlatLayer(_d, _vp, _vs, _dens, _Qk, _Qm)
end

function cut(m::FlatLayer,
             min_depth::Union{<:Real,Nothing} = nothing,
             max_depth::Union{<:Real,Nothing} = nothing)
    maxi = isnothing(max_depth) ? length(m.d) : findlast(m.d .< max_depth)
    mini = isnothing(min_depth) ? 1 : findlast(m.d .<= min_depth)
    tm = cut(m, mini:maxi)
    if isnothing(min_depth)
        return tm
    else
        td = deepcopy(tm.d)
        td[1] = min_depth
        return update(tm; d = td)
    end
end

function overlay(m1::FlatLayer, m2::FlatLayer)
    if m1.d[1] > m2.d[1]
        return m2
    elseif m1.d[end] < m2.d[1]
        return FlatLayer([m1.d; m2.d], [m1.vp; m2.vp], [m1.vs; m2.vs], [m1.dens; m2.dens], [m1.Qk; m2.Qk],
                         [m1.Qm; m2.Qm])
    else
        m3 = cut(m1, nothing, m2.d[1])
        return overlay(m3, m2)
    end
end

struct SphereLayer <: AbstractVelocityModel
    r::Vector{Float64}
    vp::Vector{Float64}
    vs::Vector{Float64}
    dens::Vector{Float64}
    Qk::Vector{Float64}
    Qm::Vector{Float64}
end

function SphereLayer(r::Vector{<:Real} = Float64[],
                     vp::Vector{<:Real} = Float64[],
                     vs::Vector{<:Real} = Float64[],
                     dens::Vector{<:Real} = Float64[],
                     Qk::Vector{<:Real} = Float64[],
                     Qm::Vector{<:Real} = Float64[])
    return SphereLayer(Float64.(r), Float64.(vp), Float64.(vs), Float64.(dens), Float64.(Qk), Float64.(Qm))
end

function SphereLayer(vel::Matrix{<:Real}; vp::Real = 5.0, vs::Real = 3.0, dens::Real = 2.7, Qk::Real = 300.0,
                     Qm::Real = 150.0)
    nlayer = size(vel, 1)
    if size(vel, 2) == 1
        return SphereLayer(vel[:, 1], fill(vp, nlayer), fill(vs, nlayer), fill(dens, nlayer), fill(Qk, nlayer),
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
        n = read(io, Int64)
        if n > 0
            r = _read_f64_vector(io)
        else
            r = Float64[]
        end
        if n > 1
            vp = _read_f64_vector(io)
        else
            vp = Float64[]
        end
        if n > 2
            vs = _read_f64_vector(io)
        else
            vs = Float64[]
        end
        if n > 3
            dens = _read_f64_vector(io)
        else
            dens = Float64[]
        end
        if n > 4
            Qk = _read_f64_vector(io)
        else
            Qk = Float64[]
        end
        if n > 5
            Qm = _read_f64_vector(io)
        else
            Qm = Float64[]
        end
        close(io)
        return SphereLayer(r, vp, vs, dens, Qk, Qm)
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
    vp::Array{Float64}
    vs::Array{Float64}
    dens::Array{Float64}
    Qk::Array{Float64}
    Qm::Array{Float64}
end

function PsudoRegular3D(x::Vector{<:Real} = Float64[],
                        y::Vector{<:Real} = Float64[],
                        d::Array{<:Real} = zeros(0, 0, 0),
                        vp::Array{<:Real} = zeros(0, 0, 0),
                        vs::Array{<:Real} = zeros(0, 0, 0),
                        dens::Array{<:Real} = zeros(0, 0, 0),
                        Qk::Array{<:Real} = zeros(0, 0, 0),
                        Qm::Array{<:Real} = zeros(0, 0, 0))
    return PsudoRegular3D(x, y, d, vp, vs, dens, Qk, Qm)
end

function PsudoRegular3D(fp_or_nm::String)
    if isfile(fp_or_nm)
        io = open(fp_or_nm, "r")
        n = read(io, Int64)
        if n > 0
            x = _read_f64_vector(io)
        else
            x = Float64[]
        end
        if n > 1
            y = _read_f64_vector(io)
        else
            y = Float64[]
        end
        if n > 2
            d = _read_f64_array(io)
        else
            d = zeros(0, 0, 0)
        end
        if n > 3
            vp = _read_f64_array(io)
        else
            vp = zeros(0, 0, 0)
        end
        if n > 4
            vs = _read_f64_array(io)
        else
            vs = zeros(0, 0, 0)
        end
        if n > 5
            dens = _read_f64_array(io)
        else
            dens = zeros(0, 0, 0)
        end
        if n > 6
            Qk = _read_f64_array(io)
        else
            Qk = zeros(0, 0, 0)
        end
        if n > 7
            Qm = _read_f64_array(io)
        else
            Qm = zeros(0, 0, 0)
        end
        close(io)
        return PsudoRegular3D(x, y, d, vp, vs, dens, Qk, Qm)
    end
    available_velocity_models = String[]
    for md in ("CRUST1.0",)
        if isfile(joinpath(@__DIR__, "..", "external", "VelocityModel", md, "psr3d.bin"))
            push!(available_velocity_models, md)
        end
    end
    if fp_or_nm ∈ available_velocity_models
        return PsudoRegular3D(joinpath(@__DIR__, "..", "external", "VelocityModel", fp_or_nm, "psr3d.bin"))
    end
    error("Input is neither a file path nor a available model name")
    return nothing
end

function _locate_point_in_vector(x::Real, v::Vector{<:Real})
    ix = findlast(x .> v)
    if isnothing(ix)
        ix = 1
        wx = 0.0
    elseif ix == length(v)
        ix = length(v) - 1
        wx = 1.0
    else
        wx = (x - v[ix]) / (v[ix+1] - v[ix])
    end
    return (ix, wx)
end

function _bilinear_interp(m::Array, iy::Integer, wy::Real, ix::Integer, wx::Real)
    if isempty(m)
        return eltype(m)[]
    end
    w11 = m[:, iy, ix]
    w12 = m[:, iy, ix+1]
    w21 = m[:, iy+1, ix]
    w22 = m[:, iy+1, ix+1]
    return w11 .* ((1.0 - wx) * (1.0 - wy)) +
           w12 .* (wx * (1.0 - wy)) +
           w21 .* ((1.0 - wx) * wy) +
           w22 .* (wx * wy)
end

function FlatLayer(prm::PsudoRegular3D, x::Real, y::Real)
    (ix, wx) = _locate_point_in_vector(x, prm.x)
    (iy, wy) = _locate_point_in_vector(y, prm.y)
    d = _bilinear_interp(prm.d, iy, wy, ix, wx)
    vp = _bilinear_interp(prm.vp, iy, wy, ix, wx)
    vs = _bilinear_interp(prm.vs, iy, wy, ix, wx)
    dens = _bilinear_interp(prm.dens, iy, wy, ix, wx)
    Qk = _bilinear_interp(prm.Qk, iy, wy, ix, wx)
    Qm = _bilinear_interp(prm.Qm, iy, wy, ix, wx)
    return FlatLayer(d, vp, vs, dens, Qk, Qm)
end

struct PsudoIrregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    d::Matrix{Float64}
    vp::Matrix{Float64}
    vs::Matrix{Float64}
    dens::Matrix{Float64}
end

struct Regular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    vp::Array{Float64}
    vs::Array{Float64}
    dens::Array{Float64}
end

struct Irregular3D <: AbstractVelocityModel
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    vp::Matrix{Float64}
    vs::Matrix{Float64}
    dens::Matrix{Float64}
end

end