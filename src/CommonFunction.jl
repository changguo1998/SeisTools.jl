module CommonFunction

export ricker

function _linearscale!(x::AbstractVector{<:AbstractFloat}, vmin::Real = 0.0, vmax::Real = 1.0)
    v0 = minimum(x)
    k0 = (vmax - vmin) / (maximum(x) - v0)
    for i = 1:length(x)
        x[i] = k0 * (x[i] - v0) + vmin
    end
    return nothing
end

"""
    ricker(t::Real, f0::Real)

time domain 0-shift ricker wavelet
"""
function ricker(t::Real, f0::Real)
    return (1.0 - 2.0 * (π * f0 * t)^2) * exp(-(π * f0 * t)^2)
end

"""
    ricker(f::Complex, f0::Real)

frequency domain 0-shift ricker wavelet
"""
function ricker(f::Complex, f0::Real)
    return f^2 * exp(-(f / f0)^2) / (π^1.5 * f0)
end

"""
    smoothramp(t::Real, thalf::Real, lim::Real = 0.95)

time domain 0-shift smoothramp wavelet

This function using `tanh` function to generate wavelet. the lim parameter is the value when t == thalf.
"""
function smoothramp(t::Real, thalf::Real, lim::Real = 0.95)
    if lim >= 1.0
        error("lim must be smaller than 1")
    end
    return tanh(log((1 + lim) / (1 - lim)) * t / thalf / 2.0) * 0.5 + 0.5
end

"""
"""
function smoothramp(f::Complex, thalf::Real) end

function cosinewindow!(x::AbstractVector)
    lx = length(x)
    for i = 1:lx
        x[i] = 0.5 * (1.0 - cos(π * (i - 1) / (lx - 1)))
    end
    return nothing
end

function cosinewindow(n::Integer)
    x = zeros(n)
    cosinewindow!(x)
    return x
end

end
