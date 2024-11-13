module Source

using LinearAlgebra, Statistics
import Base: show, isequal, +, -, Matrix, isapprox

include("basic.jl")

export MomentTensor, show, isequal, decompose, beachball_bitmap, beachball_sdrline, kagan, M0

"""
    ```
    struct MomentTensor
        values::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
    end
    ```

    the values stored (m11, m22, m33, m12, m13, m23)
"""
struct MomentTensor
    values::Tuple{Float64,Float64,Float64,Float64,Float64,Float64}
end

"""
```
MomentTensor(m11::Real, m22::Real, m33::Real, m12::Real, m13::Real, m23::Real) -> MomentTensor
```
"""
function MomentTensor(m11::Real, m22::Real, m33::Real, m12::Real, m13::Real, m23::Real)
    return MomentTensor((Float64(m11), Float64(m22), Float64(m33), Float64(m12), Float64(m13), Float64(m23)))
end

"""
```
MomentTensor(x::Union{AbstractVector{<:Real}, Tuple{<:Real,<:Real,<:Real,<:Real,<:Real,<:Real}}) -> MomentTensor
```
"""
function MomentTensor(x::Union{AbstractVector{<:Real},Tuple{<:Real,<:Real,<:Real,<:Real,<:Real,<:Real}})
    @must length(x) == 6 "length of input must be 6"
    return MomentTensor(Tuple(Float64.(x)))
end

"""
```
MomentTensor(x::AbstractMatrix{<:Real}) -> MomentTensor
```
"""
function MomentTensor(x::AbstractMatrix{<:Real})
    @must issymmetric(x) "Input must be symmetric matrix"
    @must size(x) == (3, 3) "Size of input must be (3, 3)"
    return MomentTensor(x[1, 1], x[2, 2], x[3, 3], x[1, 2], x[1, 3], x[2, 3])
end

"""
```
MomentTensor(strike::Real, dip::Real, rake::Real; m0::Real=1.0) -> MomentTensor
```

using focal mechanism format to create MomentTensor
"""
function MomentTensor(strike::Real, dip::Real, rake::Real; m0::Real = 1.0)
    m = zeros(6)
    m[1] = -m0 * (sind(2 * strike) * sind(dip) * cosd(rake) + (sind(strike))^2 * sind(2 * dip) * sind(rake))
    m[2] = m0 * (sind(2 * strike) * sind(dip) * cosd(rake) - (cosd(strike))^2 * sind(2 * dip) * sind(rake))
    m[3] = m0 * sind(2 * dip) * sind(rake)
    m[4] = m0 * (cosd(2 * strike) * sind(dip) * cosd(rake) + 0.5 * sind(2 * strike) * sind(2 * dip) * sind(rake))
    m[5] = -m0 * (cosd(strike) * cosd(dip) * cosd(rake) + sind(strike) * cosd(2 * dip) * sind(rake))
    m[6] = -m0 * (sind(strike) * cosd(dip) * cosd(rake) - cosd(strike) * cosd(2 * dip) * sind(rake))
    return MomentTensor(m)
    # (n1, n2) = _sdr2normvec(strike, dip, rake)
    # M = n2 * permutedims(n1) + n1 * permutedims(n2)
    # return MomentTensor(M)
end

function Matrix(m::MomentTensor)
    return [m.values[1] m.values[4] m.values[5];
         m.values[4] m.values[2] m.values[6];
         m.values[5] m.values[6] m.values[3]]
end

function show(io::IO, m::MomentTensor)
    show(io, Matrix(m))
end

function isequal(m1::MomentTensor, m2::MomentTensor)
    return all(map(isequal, m1.values, m2.values))
end

for sym in (:(+), :(-))
    @eval begin
        function $(sym)(m1::MomentTensor, m2::MomentTensor)
            return MomentTensor(m1.values .+ m2.values)
        end
    end
end

@doc raw"""
```
Fnorm(m::MomentTensor)
```

return Frobenius norm of `MomentTensor`
```math
F(M) = \sqrt{\sum_{i=1,3}\sum_{j=1,3}m_{ij}^2}
```
"""
function Fnorm(m::MomentTensor)
    return sqrt(sum(abs2, m.values) + sum(abs2, m.values[4:end]))
end

@doc raw"""
```
M0(m::MomentTensor)
```

return `M0` of `MomentTensor`
```math
F(M) = \frac{1}{\sqrt{2}}\sqrt{\sum_{i=1,3}\sum_{j=1,3}m_{ij}^2}
```
"""
M0(m::MomentTensor) = Fnorm(m)/sqrt(2)

"""
```
decompose_eigen(v, method) -> (iso, dc, clvd)
```
method:
    - `:DC_DC` major/minor double couple
    - `:DC_CLVD_1` double couple and clvd with sum of norm 1 == 1
    - `:DC_CLVD_2` double couple and clvd with sum of norm 2 == 1
"""
function decompose_eigen(v::Vector{<:Real}; method::Symbol=:DC_CLVD_2)
    @must issorted(v)
    (m1, m2, m3) = v
    iso = (m1 + m2 + m3) / 3.0
    d1 = (2*m1 - m2 - m3) / 3.0
    d2 = (2*m2 - m3 - m1) / 3.0
    d3 = (2*m3 - m1 - m2) / 3.0
    if method == :DC_DC
        if d2 < 0.0
            return (iso=[iso,iso,iso],dc1=[-d3,0.0,d3],dc2=[-d2,d2,0.0])
        else
            return (iso=[iso,iso,iso],dc1=[d1,0.0,-d1],dc2=[0.0,d2,-d2])
        end
    elseif method == :DC_CLVD_1
        if d2 < 0.0
            return (iso=[iso,iso,iso],dc=[d1-d2,0.0,d3+2*d2],clvd=[d2,d2,-2*d2])
        else
            return (iso=[iso,iso,iso],dc=[d1+2*d2,0.0,d3-d2],clvd=[-2*d2,d2,d2])
        end
    elseif method == :DC_CLVD_2
        return (iso=[iso,iso,iso],dc=[d1+d2/2,0.0,d3+d2/2],clvd=[-d2/2,d2,-d2/2])
    else
        error("illegal method $method")
    end
end

@doc raw"""
```
decompose(m; method) -> (iso, dc, clvd)/(iso, dc1, dc2)
```

Decompose `MomentTensor` `m` into ``M_{ISO}``, ``M_{DC}``(double couple) and
``M_{CLVD}`` with same coordinate system.
See function `decompose_eigen` for more detail

"""
function decompose(m::MomentTensor; method::Symbol=:DC_CLVD_2)
    (v, P) = eigen(Matrix(m))
    PT = permutedims(P)
    decomposed = decompose_eigen(v; method=method)
    Miso = P * diagm(decomposed.iso) * PT
    if method == :DC_DC
        Mdc1 = P * diagm(decomposed.dc1) * PT
        Mdc2 = P * diagm(decomposed.dc2) * PT
        return (iso = MomentTensor((Miso + permutedims(Miso)) ./ 2),
                dc1 = MomentTensor((Mdc1 + permutedims(Mdc1)) ./ 2),
                dc2 = MomentTensor((Mdc2 + permutedims(Mdc2)) ./ 2))
    else
        Mdc = P * diagm(decomposed.dc) * PT
        Mclvd = P * diagm(decomposed.clvd) * PT
        return (iso = MomentTensor((Miso + permutedims(Miso)) ./ 2),
                dc = MomentTensor((Mdc + permutedims(Mdc)) ./ 2),
                clvd = MomentTensor((Mclvd + permutedims(Mclvd)) ./ 2))
    end
end

function _get_eigen_angle(u1, u3, v1, v3)
    u2 = cross(u3, u1)
    v2 = cross(v3, v1)
    c = (tr([u1 u2 u3]*[v1 v2 v3]') - 1.0)*0.5
    return acosd(max(-1.0, min(1.0, c)))
end

"""
```
function kagan(mt1::MomentTensor, mt2::MomentTensor)
```

kagan angle of the principle axis between mt1 and mt2.

see Kagan, Y.Y. 1991. 3-D rotation of double-couple earthquake sources 106, 709–716.
DOI: 10.1111/j.1365-246X.1991.tb06343.x
"""
function kagan(mt1::MomentTensor, mt2::MomentTensor)
    m1 = Matrix(mt1)
    m2 = Matrix(mt2)
    E1 = eigen(m1)
    E2 = eigen(m2)

    if norm(E1.values-E2.values) > 1e-5
        error("vector position is not correct")
    end

    return min(
        _get_eigen_angle(E1.vectors[:,1], E1.vectors[:,3], E2.vectors[:,1], E2.vectors[:,3]),
        _get_eigen_angle(E1.vectors[:,1], E1.vectors[:,3], -E2.vectors[:,1], E2.vectors[:,3]),
        _get_eigen_angle(E1.vectors[:,1], E1.vectors[:,3], E2.vectors[:,1], -E2.vectors[:,3]),
        _get_eigen_angle(E1.vectors[:,1], E1.vectors[:,3], -E2.vectors[:,1], -E2.vectors[:,3])
    )
end

#=
@doc raw"""
```
focalmechanism(m::MomentTensor; digits::Integer = 0) -> (plane1=(strike1, dip1, rake1), plane2=(strike2, dip2, rake2))
```

calculate focal mechanism expression of `MomentTensor` m's double-couple component
"""
function focalmechanism(m::MomentTensor; digits::Integer = 0)
    M = [m.values[1] m.values[4] m.values[5];
         m.values[4] m.values[2] m.values[6];
         m.values[5] m.values[6] m.values[3]]
    (v, P) = eigen(M)
    Paxis = P[:, 1]
    # Baxis = P[:, 2]
    Taxis = P[:, 3]
    n1 = normalize(Paxis + Taxis)
    n1 .*= n1[3] > 0.0 ? -1.0 : 1.0
    n2 = normalize(M * n1)
    return (plane1 = round.(_normvec2sdr(n1, n2), digits = digits),
            plane2 = round.(_normvec2sdr(n2, n1), digits = digits))
end
=#

function _linetrace(n::Vector{Float64}, theta::AbstractVector{Float64})
    xv = cross(n, [0.0, 0.0, 1.0]) |> normalize
    yv = cross(n, xv) |> normalize
    if yv[3] < 0.0
        yv .*= -1
    end
    return [xv yv] * permutedims([cosd.(theta) sind.(theta)])
end

function _projectcoor(c::Matrix{<:Real})
    xy = Tuple{Float64,Float64}[]
    for i in axes(c, 2)
        push!(xy, (c[1, i] / (1.0 + c[3, i]), c[2, i] / (1.0 + c[3, i])))
    end
    return xy
end

"""
```
beachball_sdrline(m::MomentTensor, dtheta::Real=1.0; innerdecompose::Bool=true) -> (l1=xy1, l2=xy2, edge=xy3)
```

xy? is `Vector{Tuple{Float64,Float64}}` like `[(1.0, 2.0), (2.0, 3.0)]`
"""
function beachball_sdrline(m::MomentTensor, dtheta::Real = 1.0; innerdecompose::Bool = true)
    dm = innerdecompose ? decompose(m).dc : m
    M = [dm.values[1] dm.values[4] dm.values[5];
         dm.values[4] dm.values[2] dm.values[6];
         dm.values[5] dm.values[6] dm.values[3]]
    (_, V) = eigen(M)
    P = V[:, 1]
    T = V[:, 3]
    n1 = normalize(P + T)
    n2 = normalize(P - T)
    theta = range(; start = 0.0, stop = 180.0, step = dtheta)
    trace1 = _linetrace(n1, theta)
    trace2 = _linetrace(n2, theta)
    l1 = _projectcoor(trace1)
    l2 = _projectcoor(trace2)
    return (l1 = l1, l2 = l2,
            edge = map(v -> (cosd(v), sind(v)), range(; start = 0.0, stop = 360.0, step = dtheta)))
end

"""
```
function beachball_bitmap(m::MomentTensor; resolution=(201,201)) -> Matrix
```

get a map of values to plot `MomentTensor`. the first dimension of `Matrix` is north, and the second is east
"""
function beachball_bitmap(m::MomentTensor; resolution::Tuple{<:Integer,<:Integer} = (201, 201))
    M = [m.values[1] m.values[4] m.values[5];
         m.values[4] m.values[2] m.values[6];
         m.values[5] m.values[6] m.values[3]]
    vmap = zeros(resolution)
    c = zeros(3)
    for j in axes(vmap, 2), i in axes(vmap, 1)
        x = 2.0 * (i - 1) / (resolution[1] - 1) - 1.0
        y = 2.0 * (j - 1) / (resolution[2] - 1) - 1.0
        nr = x^2 + y^2
        if nr > 1.0
            vmap[i, j] = NaN
            continue
        end
        c[3] = (1.0 - nr) / (1.0 + nr)
        r = 2.0 / (nr + 1.0)
        c[1] = r * x
        c[2] = r * y
        vmap[i, j] = c' * M * c
    end
    return vmap
end

export SDR

"""
```
struct SDR
    strike1::Float64
    dip1::Float64
    rake1::Float64
    strike2::Float64
    dip2::Float64
    rake2::Float64
    m0::Float64
end
```
"""
struct SDR
    strike1::Float64
    dip1::Float64
    rake1::Float64
    strike2::Float64
    dip2::Float64
    rake2::Float64
    m0::Float64
end

function _normvec2sdr(planenorm::AbstractVector{<:Real}, slipdirec::AbstractVector{<:Real})
    factor = planenorm[3] > 0.0 ? -1.0 : 1.0
    n1 = planenorm .* factor
    n2 = slipdirec .* factor
    strike = mod(atand(n1[2], n1[1]) - 90.0, 360.0)
    dip = 180.0 - acosd(n1[3])
    refA = [cosd(strike), sind(strike), 0.0]
    refB = cross(n1, refA)
    rake = atand(dot(n2, refB), dot(n2, refA))
    return (strike, dip, rake)
end

function _sdr2normvec(strike::Real, dip::Real, rake::Real)
    n1 = [-sind(strike)*sind(dip), cosd(strike) * sind(dip), -cosd(dip)]
    refA = [cosd(strike), sind(strike), 0.0]
    refB = normalize(cross(n1, refA))
    n2 = normalize(refA .* cosd(rake) + refB .* sind(rake))
    return (n1, n2)
end

@doc raw"""
```
function SDR(strike::Real, dip::Real, rake::Real, m0::Real=1) -> SDR
```
"""
function SDR(strike::Real, dip::Real, rake::Real, m0::Real=1)
    (n1, n2) = _sdr2normvec(strike, dip, rake)
    sdr1 = _normvec2sdr(n1, n2)
    sdr2 = _normvec2sdr(n2, n1)
    return SDR(Float64(sdr1[1]), Float64(sdr1[2]), Float64(sdr1[3]),
        Float64(sdr2[1]), Float64(sdr2[2]), Float64(sdr2[3]), Float64(m0))
end

@doc raw"""
```
SDR(m::MomentTensor) -> SDR
```

calculate focal mechanism expression of `MomentTensor` m's double-couple component
"""
function SDR(m::MomentTensor)
    dcmp = decompose(m)
    M = Matrix(dcmp.dc)
    (v, P) = eigen(M)
    Paxis = P[:, 1]
    # Baxis = P[:, 2]
    Taxis = P[:, 3]
    n1 = normalize(Paxis + Taxis)
    n1 .*= n1[3] > 0.0 ? -1.0 : 1.0
    n2 = normalize(M * n1)
    return SDR(_normvec2sdr(n1, n2)..., M0(dcmp.dc))
end

"""
```
MomentTensor(sdr::SDR) -> MomentTensor
```
"""
MomentTensor(sdr::SDR) = MomentTensor(sdr.strike1, sdr.dip1, sdr.rake1, m0=sdr.m0)

"""
```
kagan(sdrA::SDR, sdrB::SDR)
```
"""
function kagan(sdrA::SDR, sdrB::SDR)
    (nA1, nA2) = _sdr2normvec(sdrA.strike1, sdrA.dip1, sdrA.rake1)
    (nB1, nB2) = _sdr2normvec(sdrB.strike1, sdrB.dip1, sdrB.rake1)

    return min(
        _get_eigen_angle(nA1, nA2, nB1, nB2),
        _get_eigen_angle(nA1, nA2, -nB1, nB2),
        _get_eigen_angle(nA1, nA2, nB1, -nB2),
        _get_eigen_angle(nA1, nA2, -nB1, -nB2)
    )
end

function isapprox(sdrA::SDR, sdrB::SDR)
    return Base.isapprox(kagan(sdrA,sdrB), 0.0) && Base.isapprox(sdrA.m0, sdrB.m0)
end
"""
```
beachball_sdrline(m::SDR, dtheta::Real=1.0) -> (l1=xy1, l2=xy2, edge=xy3)
```
"""
beachball_sdrline(sdr::SDR, dtheta::Real = 1.0) =
    beachball_sdrline(MomentTensor(sdr), dtheta; innerdecompose=false)

"""
```
function beachball_bitmap(m::SDR; resolution=(201,201)) -> Matrix
```

get a map of values to plot `SDR`. the first dimension of `Matrix` is north, and the second is east
"""
beachball_bitmap(sdr::SDR, resolution::Tuple{<:Integer,<:Integer} = (201, 201)) =
    beachball_bitmap(MomentTensor(sdr); resolution=resolution)


end
