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


abstract type ReferenceModel <: Any end

struct ReferenceSphere <: ReferenceModel
    unit::String
    r::Float64
end

EarthSphere = ReferenceSphere("m", 6371000.0)

function distance(lat1::AbstractFloat, lon1::AbstractFloat,
    lat2::AbstractFloat, lon2::AbstractFloat, ref::ReferenceSphere=EarthSphere)
    dc = max(0.0, (Float64(1.0)-cosd(lon1-lon2))*cosd(lat1)*cosd(lat2))
    cosδ = cosd(lat1-lat2) - dc
    return ref.r * acos(cosδ)
end

distance(p1::LatLon, p2::LatLon, ref::ReferenceSphere=EarthSphere) =
    distance(p1.lat, p1.lon, p2.lat, p2.lon, ref)

function azimuth(lat1::AbstractFloat, lon1::AbstractFloat,
    lat2::AbstractFloat, lon2::AbstractFloat, ref::ReferenceSphere=EarthSphere)
    n1 = [cosd(lon1)*cosd(lat1), sind(lon1)*cosd(lat1), sind(lat1)]
    n2 = [cosd(lon2)*cosd(lat2), sind(lon2)*cosd(lat2), sind(lat2)]
    xn = [-cosd(lon1)*sind(lat1), -sind(lon1)*sind(lat1), cosd(lat1)]
    xe = [-sind(lon1), cosd(lon1), 0.0]
    return mod(atand(dot(xe, n2-n1), dot(xn, n2-n1)), 360.0)
end


struct ReferenceEllipsoid <: ReferenceModel
    unit::String
    a::Float64
    b::Float64
end

WGS84 = ReferenceEllipsoid("m", 6378137.0, 6356752.0)

#=
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
    gcarc = acos(dot(n1, n2))
    dist = R * gcarc
    t12 = normalize(n2 .- dot(n1, n2) .* n1)
    tnorth = [sin(θ1 - pi / 2) * cos(φ1), sin(θ1 - pi / 2) * sin(φ1), cos(θ1 - pi / 2)]
    teast = [cos(φ1 + pi / 2), sin(φ1 + pi / 2), 0.0]
    az = atand(t12' * teast, t12' * tnorth) |> x -> mod(x, 360.0)
    return (dist, az, gcarc)
end

distance(p1::LatLon, p2::LatLon) = distance(p1.lat, p1.lon, p2.lat, p2.lon)
=#
end
