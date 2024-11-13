module Geodesy

using LinearAlgebra, Statistics

include("basic.jl")

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

function distance(lat1::Real, lon1::Real, lat2::Real, lon2::Real, ref::ReferenceSphere=EarthSphere)
    dc = max(0.0, (Float64(1.0)-cosd(lon1-lon2))*cosd(lat1)*cosd(lat2))
    cosδ = cosd(lat1-lat2) - dc
    return ref.r * acos(cosδ)
end

distance(p1::LatLon, p2::LatLon, ref::ReferenceSphere=EarthSphere) =
    distance(p1.lat, p1.lon, p2.lat, p2.lon, ref)

function azimuth(lat1::Real, lon1::Real, lat2::Real, lon2::Real, ref::ReferenceSphere=EarthSphere)
    if ref != EarthSphere
        error("only for sphere model")
        return nothing
    end
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
    fi::Float64
end

function ReferenceEllipsoid(; unit::String="m", a::Real=0, b::Real=0, fi::Real=-1)
    if iszero(b)
        if fi >= 0
            return ReferenceEllipsoid(unit, a, (1.0-1.0/fi)*a, fi)
        else
            error("invalid input, a: $a, b: $b, fi: $fi")
        end
    end
    if iszero(a)
        error("invalid input, a: $a, b: $b, fi: $fi")
    end
    return ReferenceEllipsoid(unit, a, b, a/(a-b))
end

WGS84 = ReferenceEllipsoid(; a = 6378137.0, fi = 298.257223563)
NAD83 = ReferenceEllipsoid(; a = 6378137.0, fi = 298.257222101)
GRS80 = ReferenceEllipsoid(; a = 6378137.0, fi = 298.257222101)
NAD27 = ReferenceEllipsoid(; a = 6378206.4, fi = 294.978698214)
INT24 = ReferenceEllipsoid(; a = 6378388.0, fi = 297.000000000)
CLK66 = ReferenceEllipsoid(; a = 6378206.4, fi = 294.978698214)

function _polyval(p::Vector{Float64}, x::Float64)
    v = 0.0
    for p_element = p
        v *= x
        v += p_element
    end
    return v
end

function _utm_to_ll_coef(e::Float64, m::Int=0)
    if m == 0
        c0 = [-175/16384  0    -5/256 0  -3/64 0 -1/4 0 1;
               -105/4096  0  -45/1024 0  -3/32 0 -3/8 0 0;
               525/16384  0   45/1024 0 15/256 0    0 0 0;
              -175/12288  0  -35/3072 0      0 0    0 0 0;
              315/131072  0         0 0      0 0    0 0 0];
    elseif m == 1
        c0 = [-175/16384 0   -5/256 0  -3/64 0 -1/4 0 1;
                 1/61440 0   7/2048 0   1/48 0  1/8 0 0;
              559/368640 0   3/1280 0  1/768 0    0 0 0;
              283/430080 0 17/30720 0      0 0    0 0 0;
           4397/41287680 0        0 0      0 0    0 0 0];
    elseif m == 2
        c0 = [-175/16384 0   -5/256 0  -3/64 0 -1/4 0 1;
             -901/184320 0  -9/1024 0  -1/96 0  1/8 0 0;
             -311/737280 0  17/5120 0 13/768 0    0 0 0;
              899/430080 0 61/15360 0      0 0    0 0 0;
          49561/41287680 0        0 0      0 0    0 0 0];
    end
    return map(axes(c0, 1)) do irow
        _polyval(c0[irow,:],e)
    end
end

"""
utm2ll(x,y,f,ref=WGS84) -> (lat, lon)
"""
function utm2ll(x::Real, y::Real, f::Integer, ref::ReferenceEllipsoid=WGS84)

    D0 = 180.0/pi
    maxiter = 100;
    EPS = 1e-11

    A1 = ref.a
    F1 = ref.fi
    K0 = 0.9996
    X0 = 500000.0
    Y0 = 1e7*((f<0) ? 1 : 0);
    P0 = 0
    L0 = (6 * abs(f) - 183) / D0
    E1 = sqrt((A1^2-(A1*(1-1/F1))^2)/A1^2)
    N = K0 * A1
    # computing parameters for Mercator Transverse projection
    C = _utm_to_ll_coef(E1, 0)
    YS = Y0 - N * ( C[1] * P0 +
                    C[2] * sin(2 * P0) +
                    C[3] * sin(4 * P0) +
                    C[4] * sin(6 * P0) +
                    C[5] * sin(8 * P0)
                  )
    C = _utm_to_ll_coef(E1, 1)
    zt = (y-YS)/N/C[1] + (x-X0)/N/C[1]*im
    z = zt - C[2] * sin(2 * zt) -
             C[3] * sin(4 * zt) -
             C[4] * sin(6 * zt) -
             C[5] * sin(8 * zt)
    L = real(z)
    LS = imag(z)

    l = L0 + atan(sinh(LS)/cos(L))
    # @info "lon: $(l*D0)"
    p = asin(sin(L)/cosh(LS))

    L = log(tan(pi/4 + p/2))

    # calculates latitude from the isometric latitude
    p = 2 * atan(exp(L)) - pi/2
    p0 = NaN
    n = 0
    # @info "p0 before: $p0, p before: $p"
    while (isnan(p0) || abs(p-p0) > EPS) && (n < maxiter)
        p0 = p
        es = E1 * sin(p0)
        p = 2 * atan(((1+es)/(1-es))^(E1/2)*exp(L)) - pi/2
        # @info "itr $n lat: $(p*D0)"
        n += 1
    end
    # @info "p0 after: $p0, p after: $p"
    lat = p * D0
    lon = l * D0
    # @info join([lat, " ", lon])
    return (lat, lon)
end

function LatLon(c::UTM, ref::ReferenceEllipsoid=WGS84)
    (lat, lon) = utm2ll(c.x, c.y, c.r, ref)
    return LatLon(lat, lon, c.el)
end

"""
ll2utm(lat,lon; ref=WGS84) -> (x, y, f)
"""
function ll2utm(lat::Real, lon::Real, ref::ReferenceEllipsoid=WGS84)
    D0 = 180.0/pi
    K0 = 0.9996
    X0 = 500000.0
    A1 = ref.a
    F1 = ref.fi

    p1 = lat/D0
    l1 = lon/D0
    F0 = round((l1*D0+183)/6)
    B1 = A1 * (1 - 1/F1)
    E1 = sqrt((A1*A1-B1*B1)/(A1*A1))
    P0 = 0/D0;
    L0 = (6*F0-183)/D0
    Y0 = 1e7*(p1 < 0 ? 1 : 0)
    N = K0 * A1
    C = _utm_to_ll_coef(E1, 0)
    B = C[1] * P0 + sum(C[2:end] .* map(_t->sin(_t*P0), [2,4,6,8]))
    YS = Y0 - N*B
    C = _utm_to_ll_coef(E1, 2)
    L = log(tan(pi/4 + p1/2)*(
        ((1-E1*sin(p1))/(1+E1*sin(p1)))^(E1/2)
        ))
    z = Complex(atan(sinh(L)/cos(l1-L0)), log(tan(pi/4 + asin(sin(l1-L0)/cosh(L)) /2)))
    Z = N * C[1] * z + N * sum(C[2:end] .* map(_t->sin(_t*z), [2, 4, 6, 8]))
    x = imag(Z) + X0
    y = real(Z) + YS
    return (x, y, F0)
end

function UTM(c::LatLon, ref::ReferenceEllipsoid=WGS84)
    (x, y, r) = ll2utm(c.lat, c.lon, ref)
    return UTM(x, y, r, c.el)
end

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
