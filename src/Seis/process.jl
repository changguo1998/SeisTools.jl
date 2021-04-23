"""
detrend

    return waves that are removed trend
    =================================
    detrend(x::AbstractArray; type="lsq")
"""
function detrend(x::AbstractArray; type="lsq")
    N = length(x)
    if type == "simple"
        d = (x[N] - x[1]) .* ( (1:N) .- 1) ./ (N - 1) .+ x[1]
    elseif type == "lsq"
        xm = (N + 1.0) / 2.0
        x2m = (N + 1.0) * (2.0 * N + 1.0) / 6.0
        ym = mean(x)
        xym = mean(x .* (1:N))
        d = ((xm * ym - xym) .* (1:N) .+ xm * xym .- ym * x2m) ./ (xm * xm - x2m)
    end
    return x .- d
end

"""
taper!

    taper!(out::AbstractArray, f::function, w::AbstractArray, ratio::Real)
    similar with taper
"""
function taper!(out::AbstractArray, f::Function, w::AbstractArray, ratio::Real)
    # check parameter
    if (ratio <= 0) || (ratio > 0.5)
        error("ratio should valued in [0, 0.5].");
    end
    if (f(0.0) != 0.0) && (f(1.0) != 1.0)
        error("weight function w should satisfy: w(0)==0, w(1)==1.")
    end
    N = floor(Int, length(w) * ratio)
    weight = [f.((0:(N - 1)) ./ (N - 1));ones(length(w) - 2 * N);f.(((N - 1):-1:0) ./ (N - 1))]
    out[:] = w .* weight
    return nothing
end

"""
taper

    taper waves
    ==============================
    taper(f::Function, w::AbstractArray, ratio::Real=0.05)
    f should satisfy: f(0)=0 and f(1) = 1. for example, sin(x*pi/2)
    ratio should be between 0 and 0.5
"""
function taper(f::Function, w::AbstractArray, ratio::Real=0.05)
    v = zeros(size(w))
    taper!(v, f, w, ratio)
    return v
end

"""
filt!
"""
function filt!(fw::AbstractArray, w::AbstractArray, response::DSP.FilterType=Bandpass(0.01, 0.2),
    method::DSP.FilterCoefficients=Butterworth(3))
    DSP.filt!(fw, digitalfilter(response, method), w)
    return nothing
end

function filt!(fw::AbstractArray, w::AbstractArray, response::Union{AbstractString,Symbol}=:bandpass,
    design::Union{AbstractString,Symbol}=:butter; fs::Real=100, band::Tuple=(0.1, 40), order::Int=3)
    if (response in ("bandpass", "Bandpass", "Band", "band", "BandPass")) || (response in (:bandpass, :Bandpass, :Band, :band, :BandPass))
        responsetype = DSP.Filters.Bandpass(band[1], band[2], fs=fs)
    elseif (response in ("lowpass", "Lowpass", "Low", "low", "LowPass")) || (response in (:lowpass, :Lowpass, :Low, :low, :LowPass))
        responsetype = DSP.Filters.Lowpass(band[1], fs=fs)
    elseif (response in ("highpass", "Highpass", "High", "high", "HighPass")) || (response in (:highpass, :Highpass, :High, :high, :HighPass))
        responsetype = DSP.Filters.Highpass(band[1], fs=fs)
    else
        responsetype = DSP.Filters.Bandstop(band[1], band[2], fs=fs)
    end
    if (design in ("butter", "Butter", "butterworth", "Butterworth")) || (design in (:butter, :Butter, :butterworth, :Butterworth))
        designmethod = DSP.Filters.Butterworth(order)
    else
        designmethod = DSP.Filters.Butterworth(order)
    end
    filt!(fw, w, responsetype, designmethod)
    return nothing
end

"""
filt

    filt wave using DSP.filt function

    ==================================

    filt(w::AbstractArray, response::Union{AbstractString,Symbol}=:bandpass, design::Union{AbstractString,Symbol}=:butter;
         fs::Real=100, band::Tuple=(0.1, 40), order::Int=3)

    response can be choosen from banpass, bandstop, lowpass and highpass
    design can only be set as butter, for other filters like Chebyshev I/II, Ellipsoid will be add later
    fs sampling rate, default is 100Hz
    band must be Tuple, if you are going to use lowpass or hipass, you can set it like (1.0,)
    order is the order of the filter

    ==================================

    filt!(w::AbstractArray, response::FilterType=Bandpass(0.01, 0.2), method::FilterCoefficients=Butterworth(3))

    design the filter yourself and use `filt` to filt the wave
"""
function filt(w::AbstractArray, response::Union{AbstractString,Symbol}=:bandpass,
    design::Union{AbstractString,Symbol}=:butter; fs::Real=100, band::Tuple=(0.1, 40), order::Int=3)
    fw = zeros(size(w))
    filt!(fw, w, response, design; fs=fs,band=band, order=order)
    return fw
end

function filt(w::AbstractArray, response::FilterType=Bandpass(0.01, 0.2),
    method::FilterCoefficients=Butterworth(3))
    fw = zeros(size(w))
    DSP.filt!(fw, digitalfilter(response, method), w)
    return fw
end
