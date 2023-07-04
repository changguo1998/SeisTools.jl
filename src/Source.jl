module Source

using LinearAlgebra, Statistics
import Base: show, isequal, +, -, Matrix

include("macros.jl")

export MomentTensor, show, isequal, decompose, focalmechanism, beachball_bitmap, beachball_sdrline

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

"""
```
MomentTensor(strike::Real, dip::Real, rake::Real; scale::Real=1.0) -> MomentTensor
```

using focal mechanism format to create MomentTensor
"""
function MomentTensor(strike::Real, dip::Real, rake::Real; scale::Real = 1.0)
    m = zeros(6)
    m[1] = -scale * (sind(2 * strike) * sind(dip) * cosd(rake) + (sind(strike))^2 * sind(2 * dip) * sind(rake))
    m[2] = scale * (sind(2 * strike) * sind(dip) * cosd(rake) - (cosd(strike))^2 * sind(2 * dip) * sind(rake))
    m[3] = scale * sind(2 * dip) * sind(rake)
    m[4] = scale * (cosd(2 * strike) * sind(dip) * cosd(rake) + 0.5 * sind(2 * strike) * sind(2 * dip) * sind(rake))
    m[5] = -scale * (cosd(strike) * cosd(dip) * cosd(rake) + sind(strike) * cosd(2 * dip) * sind(rake))
    m[6] = -scale * (sind(strike) * cosd(dip) * cosd(rake) - cosd(strike) * cosd(2 * dip) * sind(rake))
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

@doc raw"""
```
decompose(m::MomentTensor) -> (iso=MomentTensor, dc=MomentTensor, clvd=MomentTensor)
```

Decompose `MomentTensor` `m` into iso, double couple and CLVD with same coordinate system.
The function run as steps below:

1. get eigen value of m as ``e_1 ≤ e_2 ≤ e_3`` and eigen vector matrix P
2. iso = (e1 + e2 + e3) / 3.0
3.
```math
clvd = \left(\begin{pmatrix}
(e1+e3-2e2)/6 & & \\
& (2e2-e1-e3)/3 & \\
&  & (e1+e3-2e2)/6
\end{pmatrix}\right)
```
4. dc  =
```math
dc = \left(\begin{pmatrix}
(e1-e3)/2 & & \\
& 0 & \\
&  & (e3-e1)/2
\end{pmatrix}\right)
```
5. Miso = P[iso, iso, iso]P'
6. Mclvd = PclvdP'
7. Mdc = PdcP'

This decomposation keep that ``F(M)^2 = F(Miso)^2 + F(Mclvd)^2 + F(Mdc)^2``
"""
function decompose(m::MomentTensor)
    (v, P) = eigen(Matrix(m))
    PT = permutedims(P)
    iso = mean(v)
    c2 = (2 * v[2] - v[1] - v[3]) / 3.0
    c1 = -c2 / 2.0
    c3 = -c2 / 2.0
    s1 = (v[1] - v[3]) / 2.0
    s2 = 0.0
    s3 = (v[3] - v[1]) / 2.0
    Miso = P * diagm([iso, iso, iso]) * PT
    Mdc = P * diagm([s1, s2, s3]) * PT
    Mclvd = P * diagm([c1, c2, c3]) * PT
    return (iso = MomentTensor((Miso + permutedims(Miso)) ./ 2),
            dc = MomentTensor((Mdc + permutedims(Mdc)) ./ 2),
            clvd = MomentTensor((Mclvd + permutedims(Mclvd)) ./ 2))
end

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

end
