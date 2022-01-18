module SSP
using Statistics, FFTW

function detrend!(x::AbstractVector; type::AbstractString = "LeastSquare")
    N = length(x)
    ZERO = convert(eltype(x), 0.0)
    if type == "Mean"
        k = ZERO
        b = mean(x)
    elseif type == "SimpleLinear"
        k = (x[end] - x[1]) / (N - 1)
        b = x[1]
    elseif type == "LeastSquare"
        xm = (N + 1.0) / 2.0
        x2m = (N + 1.0) * (2.0 * N + 1.0) / 6.0
        ym = mean(x)
        xym = ZERO
        for i = 1:N
            xym += x[i] * i
        end
        xym /= N
        k = (xm * ym - xym) / (xm^2 - x2m)
        b = (xm * xym - ym * x2m) / (xm^2 - x2m)
    else
        k = ZERO
        b = ZERO
    end
    for i = 1:N
        x[i] -= k * i + b
    end
    return nothing
end

function detrend(x::AbstractVector; type::AbstractString = "LeastSquare")
    y = deepcopy(x)
    detrend!(y; type = type)
    return y
end

function taper!(f::Function, x::AbstractVector; ratio::Real = 0.05)
    @assert ((ratio >= 0.0) && (ratio <= 0.5)) "ratio should between 0 and 0.5"
    @assert ((f(0.0) == 0.0) && (f(1.0) == 1.0)) "weight function w(x) should satisfy: w(0)==0, w(1)==1"
    N = length(x)
    M = round(Int, N * ratio)
    for i = 1:M
        x[i] *= f((i - 1) / (M - 1))
        x[N-i+1] *= f((i - 1) / (M - 1))
    end
    return nothing
end

function taper(f::Function, x::AbstractVector; ratio::Real = 0.05)
    y = deepcopy(x)
    taper!(f, x; ratio = ratio)
    return y
end

function taper!(x::AbstractVector; ratio::Real = 0.05)
    taper!(identity, x; ratio = ratio)
    return nothing
end

function taper(x::AbstractVector; ratio::Real = 0.05)
    return taper(identity, x; ratio = ratio)
end

function bandpass(x::AbstractVector, w1::AbstractFloat, w2::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    @assert fs > 0.0
    ftr = digitalfilter(Bandpass(w1, w2; fs = fs), Butterworth(n))
    return filtfilt(ftr, x)
end

function bandpass!(x::AbstractVector, w1::AbstractFloat, w2::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    y = bandpass(x, w1, w2, fs; n = n)
    x .= y
    return nothing
end

function lowpass(x::AbstractVector, w::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    @assert fs > 0.0
    ftr = digitalfilter(Lowpass(w; fs = fs), Butterworth(n))
    return filtfilt(ftr, x)
end

function lowpass!(x::AbstractVector, w::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    y = lowpass(x, w, fs; n = n)
    x .= y
    return nothing
end

function highpass(x::AbstractVector, w::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    @assert fs > 0.0
    ftr = digitalfilter(Highpass(w; fs = fs), Butterworth(n))
    return filtfilt(ftr, x)
end

function highpass!(x::AbstractVector, w::AbstractFloat, fs::Real = 0.0; n::Int = 4)
    y = highpass(x, w, fs; n = n)
    x .= y
    return nothing
end

function trans(x::AbstractVector, from::NamedTuple{(:z, :p, :k),Tuple{Vector{T},Vector{T},S}},
               to::NamedTuple{(:z, :p, :k),Tuple{Vector{T},Vector{T},S}}, fs::Real) where {T<:Complex,S<:Real}
    n = length(x)
    m = floor(Int, n / 2)
    X = fft(x)
    Y = zeros(ComplexF64, n)
    rmax = -Inf
    for i = 1:m
        ω = 2 * π * 1im * i * fs / n
        r = to.k / from.k
        for p in from.p
            r *= ω - p
        end
        for z in to.z
            r *= ω - z
        end
        for p in to.p
            r /= ω - p
        end
        for z in from.z
            r /= ω - z
        end
        rmax = max(rmax, abs(r))
        Y[i] = X[i]*r
        Y[n-i+1] = conj(Y[i])
    end
    if rmax > 1e5
        @warn "Waveform in some frequency is scaled larger than 1e5"
    end
    ifft!(Y)
    y = real.(Y)
    return y
end
end
