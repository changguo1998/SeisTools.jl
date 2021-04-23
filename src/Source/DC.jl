function tpd2xyz(t, p)
    return [sind.(t) .* cosd.(p), sind.(t) .* sind.(p), cosd.(t)]
end

function xyz2tpd(x, y, z)
    t = acosd(z)
    p = atand(y, x)
    return [t, p]
end

function normalvector(m::DoubleCouple)
    v1 = tpd2xyz(180 - m.dip, 90 + m.strike)
    refA = tpd2xyz(90, m.strike)
    refB = cross(v1, refA[:])
    v2 = refA .* cosd(m.rake) .+ refB .* sind(m.rake)
    return (v1, v2)
end

function normal2DC(v1::AbstractArray, v2::AbstractArray)
    if v1[3] > 0
        v1 .= -v1
        v2 .= -v2
    end
    tp1 = xyz2tpd(v1[1], v1[2], v1[3])
    s = mod(tp1[2] - 90, 360)
    d = 180 - tp1[1]
    refv1 = tpd2xyz(90, s)
    refv2 = cross(v1, refv1)
    r = atand(v2' * refv2, v2' * refv1)
    return DoubleCouple(s, d, r)
end