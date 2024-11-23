module DataProcess

using Statistics, FFTW, Dates, LinearAlgebra, DSP

TIME_PRECISION = Microsecond(1)
TIME_PRECISION_SECOND_RATIO = Second(1)/TIME_PRECISION

_Second(t::Real) = round(Int, t*TIME_PRECISION_SECOND_RATIO)*TIME_PRECISION

export @linearscale

include("basic.jl")

export detrend, detrend!

"""
```
detrend!(x::AbstractVecOrMat; type=:LeastSquare)
```

Remove the linear content of x. See `detrend` for more information
"""
function detrend!(x::AbstractVecOrMat; type::Symbol = :LeastSquare)
    @must (type in (:LeastSquare, :Mean, :SimpleLinear)) "type must be one of :LeastSquare, Mean or SimpleLinear"
    N = size(x, 1)
    ZERO = convert(eltype(x), 0.0)
    for col in eachcol(x)
        if type == :Mean
            k = ZERO
            b = mean(col)
        elseif type == :SimpleLinear
            k = (col[end] - col[1]) / (N - 1)
            b = col[1]
        elseif type == :LeastSquare
            xm = (N + 1.0) / 2.0
            x2m = (N + 1.0) * (2.0 * N + 1.0) / 6.0
            ym = mean(col)
            xym = ZERO
            for i = 1:N
                xym += col[i] * i
            end
            xym /= N
            k = (xm * ym - xym) / (xm^2 - x2m)
            b = (xm * xym - ym * x2m) / (xm^2 - x2m)
        else
            k = ZERO
            b = ZERO
        end
        for i = 1:N
            col[i] -= k * i + b
        end
    end
    return nothing
end

"""
```
detrend(x::AbstractVecOrMat; type=:LeastSquare)
```

Remove the linear content of x. The linear type can be

  - `:LeastSquare`(default) using Least Square method to get the linear content
  - `:Mean` using mean of x
  - `:SimpleLinear` using the first and last sample to get linear content
"""
function detrend(x::AbstractVecOrMat; type::Symbol = :LeastSquare)
    y = deepcopy(x)
    detrend!(y; type = type)
    return y
end

export taper, taper!

"""
```
taper!(f::Function, x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
```

see `taper`
"""
function taper!(f::Function, x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
    @must ((ratio >= 0.0) && (ratio <= 0.5)) "ratio should between 0 and 0.5"
    @must ((f(0.0) == 0.0) && (f(1.0) == 1.0)) "weight function w(x) should satisfy: w(0)==0, w(1)==1"
    @must (side in (:Head, :Tail, :Both)) "specify which side to be tapered"
    N = size(x, 1)
    M = max(1, round(Int, N * ratio))
    if M == 1
        for xc in eachcol(x)
            if (side == :Head) || (side == :Both)
                xc[1] *= 0.0
            end
            if (side == :Tail) || (side == :Both)
                xc[end] *= 0.0
            end
        end
        return nothing
    end
    for xc in eachcol(x)
        if (side == :Head) || (side == :Both)
            for i = 1:M
                xc[i] *= f((i - 1) / (M - 1))
            end
        end
        if (side == :Tail) || (side == :Both)
            for i = 1:M
                xc[N-i+1] *= f((i - 1) / (M - 1))
            end
        end
    end
    return nothing
end

"""
```
taper!(x::AbstractVecOrMat; ratio::Real=0.05, side::Symbol=:Both)
```
"""
function taper!(x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
    taper!(identity, x; ratio = ratio, side = side)
    return nothing
end

"""
```
taper(f::Function, x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
```

Add taper to waveform.

  - `f` is window function. default is f(x) = x.
  - `ratio` is window length, between 0 - 0.5
  - `side` should be one of `:Both`, `:Head` or `Tail`
"""
function taper(f::Function, x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
    y = deepcopy(x)
    taper!(f, y; ratio = ratio, side = side)
    return y
end

function taper(x::AbstractVecOrMat; ratio::Real = 0.05, side::Symbol = :Both)
    y = deepcopy(x)
    taper!(identity, y; ratio = ratio, side = side)
    return y
end

export bandpass, highpass, lowpass

"""
```
bandpass(x::AbstractVector, w1::Real, w2::Real, fsample::Real = 0.0; n::Int = 4)
```
"""
function bandpass(x::AbstractVecOrMat, w1::Real, w2::Real, fsample::Real = 0.0; n::Int = 4)
    @must fsample > 0.0
    ftr = digitalfilter(Bandpass(w1, w2; fs = fsample), Butterworth(n))
    return filtfilt(ftr, x)
end

# function bandpass(x::AbstractVecOrMat, w1::Real, w2::Real, fs::Real = 1.0)
#     y = deepcopy(x)
#     L = size(x, 1)
#     df = fs / L
#     for ix in axes(x, 2)
#         X = fft(x[:, ix])
#         Y = zeros(eltype(X), length(X))
#         for i = 1:round(Int, L / 2)
#             f = (i - 1) * df
#             if (f >= w1) && (f <= w2)
#                 Y[i] = X[i]
#                 Y[L-i+1] = conj(X[i])
#             end
#         end
#         ifft!(Y)
#         y[:, ix] .= real.(Y)
#     end
#     return y
# end

"""
```
lowpass(x::AbstractVector, w::Real, fsample::Real = 0.0; n::Int = 4)
```
"""
function lowpass(x::AbstractVector, w::Real, fsample::Real = 0.0; n::Int = 4)
    @must fsample > 0.0
    ftr = digitalfilter(Lowpass(w; fs = fsample), Butterworth(n))
    return filtfilt(ftr, x)
end

"""
```
highpass(x::AbstractVector, w::Real, fsample::Real = 0.0; n::Int = 4)
```
"""
function highpass(x::AbstractVector, w::Real, fsample::Real = 0.0; n::Int = 4)
    @must fsample > 0.0
    ftr = digitalfilter(Highpass(w; fs = fsample), Butterworth(n))
    return filtfilt(ftr, x)
end

export ZPK, trans, DISPLACEMENT, VELOCITY, ACCELERATION
struct ZPK
    z::Vector{ComplexF64}
    p::Vector{ComplexF64}
    k::Float64
end

function ZPK(z::Vector = ComplexF64[], p::Vector = ComplexF64[], k::Real = 1.0)
    return ZPK(ComplexF64.(z), ComplexF64.(p), Float64(k))
end

const DISPLACEMENT = ZPK()
const VELOCITY = ZPK([0.0])
const ACCELERATION = ZPK([0.0, 0.0])

"""
```
trans(x::AbstractVector, from::ZPK, to::ZPK, fs::Real)
```

add or remove response. the response is applied in zpk form
"""
function trans(x::AbstractVector, from::ZPK, to::ZPK, fs::Real)
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
        Y[i] = X[i] * r
        if i > 1
            Y[n-i+2] = conj(Y[i])
        end
    end
    Y[m+1] = 0.0
    # if rmax > 1e5
    #     @warn "Waveform in some frequency is scaled larger than 1e5"
    # end
    ifft!(Y)
    y = real.(Y)
    return y
end

"""
```
conv_f!(z::AbstractArray{<:Real}, x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real})
```

Convolve `x` and `y` at frequency domain, and save result into `z`.

The first dimension will be consider as time, and other dimension will be iterated.
"""
function conv_f!(z::AbstractArray{<:Real}, x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real})
    L = size(x, 1) + size(y, 1) - 1
    Lz = min(size(z, 1), L)
    z .= convert(eltype(z), 0)
    Threads.@threads for (xcol, ycol) in Tuple.(CartesianIndices((size(x, 2), size(y, 2))))
        X = zeros(ComplexF64, L)
        Y = zeros(ComplexF64, L)
        for xr in axes(x, 1)
            X[xr] = x[xr, xcol]
        end
        for yr in axes(y, 1)
            Y[yr] = y[yr, ycol]
        end
        fft!(X)
        fft!(Y)
        X .*= Y
        ifft!(X)
        for zr = 1:Lz
            z[zr, xcol, ycol] = real(X[zr])
        end
    end
    return nothing
end

"""
```
conv_f(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real}) -> Array
```

Convolve each vector from the first dimension of `x` and `y` at frequency domain,
save result into `z`.

The first dimension will be consider as time, and other dimension will be iterated.
Size of `z` will be `(size(x, 1)+size(y, 1)-1, size(x, 2), size(y, 2))`.
The result of `conv(x[:, ix], y[:, iy])` will be stored at `z[:, ix, iy]`
"""
function conv_f(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real})
    z = zeros(promote_type(eltype(x), eltype(y)), size(x, 1) + size(y, 1) - 1, size(x, 2), size(y, 2))
    conv_f!(z, x, y)
    return z
end

"""
```
conv_t!(z::AbstractArray{<:Real}, x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real})
```

Convolve `x` and `y` at time domain and save result into `z`.

The first dimension will be consider as time, and other dimension will be iterated.
The first `size(z, 1)` element will be stored in `z`
"""
function conv_t!(z::AbstractArray{<:Real}, x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real})
    Lx = size(x, 1)
    Ly = size(y, 1)
    ZERO = convert(eltype(z), 0)
    Threads.@threads for (xcol, ycol) in Tuple.(CartesianIndices((size(x, 2), size(y, 2))))
        for zrow in axes(z, 1)
            z[zrow, xcol, ycol] = ZERO
            for xrow = 1:min(Lx, zrow)
                yrow = zrow - xrow + 1
                if yrow > Ly
                    continue
                end
                z[zrow, xcol, ycol] += x[xrow, xcol] * y[yrow, ycol]
            end
        end
    end
    return nothing
end

"""
```
conv_t(x, y, len::Integer=size(x, 1) + size(y, 1) - 1) -> Array
```

convolve each vector from the first dimension of `x` and `y` at time domain,
save the first `len` result into `z`. the first dimension will be consider as time,
and other dimension will be iterated. Size of `z` will be `(len, size(x, 2), size(y, 2))`.
The result of `conv(x[:, ix], y[:, iy])` will be stored at `z[:, ix, iy]`
"""
function conv_t(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real}, len::Integer = size(x, 1) + size(y, 1) - 1)
    z = zeros(promote_type(eltype(x), eltype(y)), len, size(x, 2), size(y, 2))
    conv_t!(z, x, y)
    return z
end

"""
```
xcorr_t!(r, x, y, yshiftstart::Integer)
```

calculate time domain cross correlation of `x` and `y` (shift `y` relative to `x`), store the result in `r`.
the first element of `r` is related to `yshiftstart`
"""
function xcorr_t!(r::AbstractArray{<:Real}, x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real},
                  yshiftstart::Integer)
    Lx = size(x, 1)
    Ly = size(y, 1)
    ZERO = convert(eltype(r), 0)
    Threads.@threads for (xcol, ycol) in Tuple.(CartesianIndices((size(x, 2), size(y, 2))))
        for rrow in axes(r, 1)
            r[rrow, xcol, ycol] = ZERO
            shift = yshiftstart + rrow - 1
            for xrow = max(1, shift + 1):Lx
                yrow = xrow - shift
                if yrow > Ly
                    break
                end
                r[rrow, xcol, ycol] += x[xrow, xcol] * y[yrow, ycol]
            end
        end
    end
    return nothing
end

"""
```
xcorr_t(x, y, shiftstart::Integer=1-size(y, 1), shiftend::Integer=size(x, 1)-1) -> (lag, corr)
```

calculate time domain cross correlation of `x` and `y` (shift `y` relative to `x`).
the shift stored in lag, showing how much the `y` shift relative to `x`
"""
function xcorr_t(x::AbstractVecOrMat{<:Real}, y::AbstractVecOrMat{<:Real}, shiftstart::Integer = 1 - size(y, 1),
                 shiftend::Integer = size(x, 1) - 1)
    @must shiftend > shiftstart "shiftend must be larger than shiftstart"
    z = zeros(promote_type(eltype(x), eltype(y)), shiftend - shiftstart + 1, size(x, 2), size(y, 2))
    xcorr_t!(z, x, y, shiftstart)
    return (lag = range(shiftstart; step = 1, length = shiftend - shiftstart + 1), c = z)
end

"""
```
fap(x::AbstractVecOrMat{<:Real}, dt::Real) -> (f, Amplitude, Phase)
```
calculate amplitude and phase vs frequency
"""
function fap(x::AbstractVecOrMat{<:Real}, dt::Real)
    f = range(0.0, 1.0, length=size(x, 1))./dt
    X = fft(x)
    A = abs.(X)
    p = angle.(X)
    return (f, A, p)
end

export fap, resample, resample!

"""
```
resample!(y::AbstractVecOrMat{<:Real}, x::AbstractVecOrMat{<:Real})
```
"""
function resample!(y::AbstractVecOrMat{<:Real}, x::AbstractVecOrMat{<:Real})
    @must size(y, 2) == size(x, 2) "column of x and y must be equal"
    X = zeros(Complex{eltype(x)}, size(x))
    Y = zeros(Complex{eltype(x)}, size(y))
    Lx = size(X, 1)
    Ly = size(Y, 1)
    X .= x
    for col in eachcol(X)
        fft!(col)
    end
    Y[1, :] .= X[1, :]
    for col in axes(y, 2), row = 2:floor(Int, min(Lx, Ly) / 2)
        Y[row, col] = X[row, col]
        Y[Ly-row+2, col] = conj(X[row, col])
    end
    for col in eachcol(Y)
        ifft!(col)
    end
    y .= real.(Y) .* (Ly / Lx)
    return nothing
end

"""
```
resample(x::AbstractVecOrMat{<:Real}, N::Integer) -> VecOrMat
```
"""
function resample(x::AbstractVecOrMat{<:Real}, N::Integer)
    s = collect(size(x))
    s[1] = N
    y = zeros(eltype(x), Tuple(s))
    resample!(y, x)
    return y
end

"""
```
resample(x::AbstractVecOrMat{<:Real}, ratio::Real) -> VecOrMat
```
"""
resample(x::AbstractVecOrMat{<:Real}, ratio::Real) = resample(x, round(Int, size(x, 1) * ratio))

"""
```
resample(x::AbstractVecOrMat{<:Real}, xdt::Real, ydt::Real) -> VecOrMat
```
"""
resample(x::AbstractVecOrMat{<:Real}, xdt::Real, ydt::Real) = resample(x, round(Int, size(x, 1) * xdt / ydt))

export cut, cut!

"""
```
cut!(y::AbstractVecOrMat{<:Real}, x::AbstractVecOrMat{<:Real}, startx::Integer, starty::Integer, len::Integer)
```

copy rows start at `startx` to y. if x is too short, stop at the end of x, else stop when copied `len` row
"""
function cut!(y::AbstractVecOrMat{<:Real}, x::AbstractVecOrMat{<:Real}, startx::Integer, starty::Integer, len::Integer)
    Lx = size(x, 1)
    Ly = size(y, 1)
    L = min(Lx - startx + 1, len, Ly - starty + 1)
    for col in axes(y, 2)
        for row = 1:L
            y[row+starty-1, col] = x[row+startx-1, col]
        end
    end
    return nothing
end

"""
```
cut(x::AbstractVecOrMat{<:Real}, startx::Integer, len::Integer; fillval=0.0) -> VecOrMat
```

cut `len` rows from x while start at `startx` row
"""
function cut(x::AbstractVecOrMat{<:Real}, startx::Integer, len::Integer; fillval = 0.0)
    if typeof(x) <: AbstractVector
        y = fill(convert(eltype(x), fillval), len)
    else
        y = fill(convert(eltype(x), fillval), len, size(x, 2))
    end
    if startx < 1
        if len - 1 + startx < 1
            return y
        end
        cut!(y, x, 1, 2 - startx, len - 1 + startx)
    else
        cut!(y, x, startx, 1, len)
    end
    return y
end

"""
```
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, len::Period, dt::TimePeriod; fillval=0.0)
-> (ybegin::DateTime, y::VecOrMat, dt::Millisecond)
```

cut rows from x, start from time `start` with time length `len`
"""
function cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, len::Period, dt::TimePeriod;
             fillval = 0.0)
    npts = round(Int, len / dt)
    xstart = round(Int, (start - xbegin) / dt) + 1
    y = cut(x, xstart, npts; fillval = fillval)
    return (xbegin + (xstart - 1) * dt, y, dt)
end

"""
```
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, len::Real, dt::Real; fillval=0.0)
-> (ybegin::DateTime, y::VecOrMat, dt::Millisecond)
```

cut rows from x, start from time `start` with time length `len`. `len` and `dt` is supposed to use unit `s`
"""
function cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, len::Real, dt::Real;
             fillval = 0.0)
    npts = round(Int, _Second(len) / _Second(dt))
    xstart = round(Int, (start - xbegin) / _Second(dt)) + 1
    y = cut(x, xstart, npts; fillval = fillval)
    return (xbegin + (xstart - 1) * _Second(dt), y, _Second(dt))
end

"""
```
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, stop::DateTime, dt::TimePeriod; fillval=0.0)
-> (ybegin::DateTime, y::VecOrMat, dt::Millisecond)
```

cut rows from x, start from time `start` to time `stop`
"""
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, stop::DateTime, dt::TimePeriod;
fillval = 0.0) = cut(x, xbegin, start, stop - start, dt; fillval = fillval)

"""
```
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, stop::DateTime, dt::Real; fillval=0.0)
-> (ybegin::DateTime, y::VecOrMat, dt::Millisecond)
```

cut rows from x, start from time `start` to time `stop`
"""
cut(x::AbstractVecOrMat{<:Real}, xbegin::DateTime, start::DateTime, stop::DateTime, dt::Real;
fillval = 0.0) = cut(x, xbegin, start, stop - start, _Second(dt); fillval = fillval)

export merge

"""
```
merge(x::Vector{<:AbstractVecOrMat{<:Real}}, xbegins::Vector{DateTime}, dt::TimePeriod; fillval=0.0)
-> (ybegin::DateTime, y::VecOrMat, dt::TimePeriod)
```

merge each element in x into one VecOrMat y, the latter element in x will overwrite data from the former element when
there are overlapes
"""
function merge(x::Vector{<:AbstractVecOrMat{<:Real}}, xbegins::Vector{DateTime}, dt::TimePeriod; fillval = 0.0)
    xstops = map(i -> xbegins[i] + dt * size(x[i], 1), eachindex(x))
    ybegin = minimum(xbegins)
    ystop = maximum(xstops)
    len = round(Int, (ystop - ybegin) / dt)
    if typeof(x[1]) <: AbstractVector
        y = fill(fillval, len)
    else
        y = fill(fillval, len, size(x[1], 2))
    end
    for i in eachindex(x)
        shift = round(Int, (xbegins[i] - ybegin) / dt) + 1
        cut!(y, x[i], 1, shift, size(x[i], 1))
    end
    return (ybegin, y, dt)
end

end # module
