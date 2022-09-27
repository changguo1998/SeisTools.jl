module Geodesy

using LinearAlgebra, Statistics

include("macros.jl")

abstract type Point <: Any end

struct LatLon <: Point
    lat::Float64
    lon::Float64
    el::Float64
end

struct UTM <: Point
    x::Float64
    y::Float64
    r::String
    el::Float64
end

"""
    distance(lat1, lon1, lat2, lon2) -> (dist, az, gcarc)

compute distance and azimuth between two points on the Earth according to the reference Sphere.
distance is in km, az in degree, centered at point 1, and garc in radius degree
"""
function distance(lat1, lon1, lat2, lon2)
    R = 6371.0
    θ1 = deg2rad(90.0 - lat1)
    θ2 = deg2rad(90.0 - lat2)
    φ1 = deg2rad(lon1)
    φ2 = deg2rad(lon2)
    n1 = [sin(θ1) * cos(φ1), sin(θ1) * sin(φ1), cos(θ1)]
    n2 = [sin(θ2) * cos(φ2), sin(θ2) * sin(φ2), cos(θ2)]
    gcarc = acos(n1' * n2)
    dist = R * gcarc
    t12 = normalize(n2 .- (n1' * n2) .* n1)
    tnorth = [sin(θ1 - pi / 2) * cos(φ1), sin(θ1 - pi / 2) * sin(φ1), cos(θ1 - pi / 2)]
    teast = [cos(φ1 + pi / 2), sin(φ1 + pi / 2), 0.0]
    az = atand(t12' * teast, t12' * tnorth) |> x -> mod(x, 360.0)
    return (dist, az, gcarc)
end

distance(p1::LatLon, p2::LatLon) = distance(p1.lat, p1.lon, p2.lat, p2.lon)

end
