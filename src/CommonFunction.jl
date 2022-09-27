module CommonFunction

using SpecialFunctions, SpecialPolynomials

include("macros.jl")

export ricker

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

function cosinewindow(x::Real)
    return (1.0 - cospi(x))/2.0
end

function hanning(x::Real)
    return (1.0 + cospi(x))/2.0
end

end
