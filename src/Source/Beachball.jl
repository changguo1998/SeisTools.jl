function fillcolor(m::DoubleCouple, X, Y, Z)
    (v1, v2) = normalvector(m)
    t = ones(size(X))
    V1 = cat(t .* v1[1], t .* v1[2], t .* v1[3]; dims = 3)
    V2 = cat(t .* v2[1], t .* v2[2], t .* v2[3]; dims = 3)
    co = cat(X, Y, Z; dims = 3)
    C = sign.(sum(V1 .* co; dims = 3) .* sum(V2 .* co; dims = 3))
    return C[:, :, 1]
end

function fillcolor(m::MomentTensor, X, Y, Z)
    if m.coor == :NED
        C = -ones(size(X))
        for i = 1:size(X, 1), j = 1:size(X, 2)
            v = [X[i, j]; Y[i, j]; Z[i, j]]
            C[i, j] = v' * m.value * v
        end
        return C
    else
        error("The moment tensor should be expressed in NED coordinate system.")
    end
end

"""
beachball3d(mechanism)

    plot 3d beachball
    ========================================================
    beachball3d(mechanism::Mechanism; grid::AbstractFloat=1.0, part::Symbol=:lower)

    Arguments:
        **mechanism**: mechanism of source.
        **grid**: meshgrid, unit is degree
    Output:
        (x, y, c)
        x, y is the grid of unit circle
        c shows colors
"""
function beachball3d(mechanism::Mechanism; ngrid::Int = 100, part::Symbol = :lower)
    if part == :upper
        t = [range(90.0, 180.0; length = ngrid)...]
        p = [range(0.0, 360.0; length = ngrid)...]
    elseif part == :north
        t = [range(0.0, 180.0; length = ngrid)...]
        p = [range(-90.0, 90.0; length = ngrid)...]
    elseif part == :south
        t = [range(0.0, 180.0; length = ngrid)...]
        p = [range(90.0, 270.0; length = ngrid)...]
    elseif part == :east
        t = [range(0.0, 180.0; length = ngrid)...]
        p = [range(0.0, 180.0; length = ngrid)...]
    elseif part == :west
        t = [range(0.0, 180.0; length = ngrid)...]
        p = [range(180.0, 360.0; length = ngrid)...]
    else
        t = [range(0.0, 90.0; length = ngrid)...]
        p = [range(0.0, 360.0; length = ngrid)...]
    end
    T = repeat(t, 1, length(p))
    P = repeat(p', length(t), 1)
    (X, Y, Z) = tpd2xyz(T, P)
    C = fillcolor(mechanism, X, Y, Z)
    return (X, Y, Z, C)
end

"""
beachball2d(mechanism)

    plot 2d beachball
    ========================================================
    beachball2d(mechanism; grid=1)

    Arguments:
        **mechanism**: mechanism of source.
        **grid**: meshgrid, unit is degree
    Output:
        (x, y, c)
        x, y is the grid of unit circle
        c shows colors
"""
# function beachball(mechanism::Mechanism; fillcolor=:black, grid=1, show=true)
function beachball2d(mechanism::Mechanism; ngrid::Int = 100)
    t = [range(-1.0, 1.0; length = ngrid)...]
    x = zeros(ngrid, ngrid)
    y = zeros(size(x))
    z = zeros(size(x))
    for i = 1:ngrid, j = 1:ngrid
        n2 = t[i]^2 + t[j]^2
        tz = (1.0 - n2) / (1.0 + n2)
        r = 2.0 / (n2 + 1.0)
        x[i, j] = (tz >= 0) ? t[i] * r : NaN
        y[i, j] = (tz >= 0) ? t[j] * r : NaN
        z[i, j] = (tz >= 0) ? tz : NaN
    end
    c = fillcolor(mechanism, x, y, z)
    c = sign.(c)
    return (x = x, y = y, c = c)
end

function beachball(mechanism::Mechanism; ngrid::Int = 200)
    t = beachball2d(mechanism; ngrid = ngrid)
    return (x = t.x, y = t.y, c = t.c)
end
