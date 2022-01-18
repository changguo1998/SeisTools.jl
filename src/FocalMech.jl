module FocalMech

function momenttensor_r(s::Real, d::Real, r::Real)
    m = zeros(6)
    m[1] = -1 * (sin(2 * s) * sin(d) * cos(r) + (sin(s))^2 * sin(2 * d) * sin(r))
    m[2] = sin(2 * s) * sin(d) * cos(r) - (cos(s))^2 * sin(2 * d) * sin(r)
    m[3] = sin(2 * d) * sin(r)
    m[4] = cos(2 * s) * sin(d) * cos(r) + 0.5 * sin(2 * s) * sin(2 * d) * sin(r)
    m[5] = -1 * (cos(s) * cos(d) * cos(r) + sin(s) * cos(2 * d) * sin(r))
    m[6] = -1 * (sin(s) * cos(d) * cos(r) - cos(s) * cos(2 * d) * sin(r))
    return m
end

momenttensor_r(x::Tuple{R,S,T}) where {R<:Real,S<:Real,T<:Real} = momenttensor_r(x...)
momenttensor_r(x::Vector{T}) where {T<:Real} = momenttensor_r(x[1], x[2], x[3])
momenttensor_r(; strike::Real = 0.0, dip::Real = 0.0, rake::Real = 0.0) = momenttensor_r(strike, dip, rake)
momenttensor(x...) = momenttensor_r(deg2rad.(x)...)
momenttensor(; strike::Real = 0.0, dip::Real = 0.0, rake::Real = 0.0) = momenttensor_r(deg2rad(strike), deg2rad(dip),
                                                                                       deg2rad(rake))

function beachball_line() end
function beachball_map(s::Real, d::Real, r::Real; xgrid::Int = 10, ygrid::Int = 10)
    x = range(-1.0, 1.0; length = xgrid)
    y = range(-1.0, 1.0; length = ygrid)
    c = zeros(xgrid, ygrid)
    t = momenttensor(s, d, r)
end
end
