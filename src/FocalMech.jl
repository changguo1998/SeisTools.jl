module FocalMech

function momenttensor(s::T, d::T, r::T) where {T<:Real}
    m = zeros(6)
    m[1] = -1 * (sin(2 * s) * sin(d) * cos(r) + (sin(s))^2 * sin(2 * d) * sin(r))
    m[2] = sin(2 * s) * sin(d) * cos(r) - (cos(s))^2 * sin(2 * d) * sin(r)
    m[3] = sin(2 * d) * sin(r)
    m[4] = cos(2 * s) * sin(d) * cos(r) + 0.5 * sin(2 * s) * sin(2 * d) * sin(r)
    m[5] = -1 * (cos(s) * cos(d) * cos(r) + sin(s) * cos(2 * d) * sin(r))
    m[6] = -1 * (sin(s) * cos(d) * cos(r) - cos(s) * cos(2 * d) * sin(r))
    return m
end

momenttensor(x::Tuple{R,S,T}) where {R<:Real,S<:Real,T<:Real} = momenttensor(x...)
momenttensor(x::Vector{T}) where {T<:Real} = momenttensor(x[1], x[2], x[3])
momenttensor(; strike::Real = 0.0, dip::Real = 0.0, rake::Real = 0.0) = momenttensor(strike, dip, rake)
momenttensor_d(x...) = momenttensor(deg2rad.(x)...)
momenttensor_d(; strike::Real = 0.0, dip::Real = 0.0, rake::Real = 0.0) = momenttensor(deg2rad(strike), deg2rad(dip),
                                                                                       deg2rad(rake))
end
