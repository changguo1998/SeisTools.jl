function fillcolor(m::DoubleCouple, X, Y, Z)
    (v1, v2) = normalvector(m)
    t = ones(size(X))
    V1 = cat(t .* v1[1], t .* v1[2], t .* v1[3], dims=3)
    V2 = cat(t .* v2[1], t .* v2[2], t .* v2[3], dims=3)
    co = cat(X, Y, Z, dims=3)
    C = sign.(sum(V1 .* co, dims=3) .* sum(V2 .* co, dims=3))
    return C[:, :, 1]
end

function fillcolor(m::MomentTensor, X, Y, Z)
    if m.coor == :NED
        C = -ones(size(X))
        for i in 1:size(X, 1), j in 1:size(X, 2)
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
function beachball3d(mechanism::Mechanism; grid::AbstractFloat=1.0, part::Symbol=:lower)
    if part == :upper
        t = 90.0:grid:180.0
        p = 0.0:grid:360.0
    elseif part == :north
        t = 0.0:grid:180.0
        p = -90.0:grid:90.0
    elseif part == :south
        t = 0.0:grid:180.0
        p = 90.0:grid:270.0
    elseif part == :east
        t = 0.0:grid:180.0
        p = 0.0:grid:180.0
    elseif part == :west
        t = 0.0:grid:180.0
        p = 180.0:grid:360.0
    else
        t = 0.0:grid:90.0
        p = 0.0:grid:360.0
    end
    T = repeat(t,  1, length(p))
    P = repeat(p', length(t), 1)
    (X, Y, Z) = tpd2xyz(T, P)
    C = fillcolor(mechanism, X, Y, Z)
    return(X, Y, Z, C)
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
function beachball2d(mechanism::Mechanism; grid::AbstractFloat=1.0)
    (X, Y, Z, C) = beachball3d(mechanism; grid=grid, part=:lower)
    k = 1.0 ./ (1.0 .+ Z)
    x = Y .* k
    y = X .* k
    c = sign.(C)
    # if show
    #     contour(x, y, c, fill=true, levels=1, colorbar=:none, seriescolor=[:white, fillcolor], 
    #     showaxis=false, grid=:hide, aspect_ratio=:equal)
    #     t = 0.0:grid:360.0
    #     plot!(cosd.(t), sind.(t), linecolor=:black, linewidth=1, label="")
    # end
    return (x = x, y = y, c = c)
end

function beachball(mechanism::Mechanism; grid::AbstractFloat=1.0, fillcolor::Symbol=:blue)
    t = beachball2d(mechanism, grid=grid)
    contour(t.x, t.y, t.c, fill=true, levels=1, colorbar=:none, seriescolor=[:white, fillcolor], 
        showaxis=false, grid=:hide, aspect_ratio=:equal)
    t = 0.0:grid:360.0
    plot!(cosd.(t), sind.(t), linecolor=:black, linewidth=1, label="")
end