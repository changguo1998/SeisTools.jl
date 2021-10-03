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

function detrend!(f::WaveFrame; type::AbstractString = "LeastSquare")
    detrend!(f.data; type = type)
    return nothing
end

function detrend(x::Union{WaveFrame,AbstractVector}; type::AbstractString = "LeastSquare")
    g = deepcopy(f)
    detrend!(g; type = type)
    return g
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

function taper!(f::Function, frame::WaveFrame; ratio::Real = 0.05)
    taper!(f, frame.data; ratio = ratio)
    return nothing
end

function taper(f::Function, x::Union{WaveFrame,AbstractVector}; ratio::Real = 0.05)
    y = deepcopy(x)
    taper!(f, x; ratio = ratio)
    return y
end

function taper!(x::Union{WaveFrame,AbstractVector}; ratio::Real = 0.05)
    taper!(var -> var, x; ratio = ratio)
    return nothing
end

function taper(x::Union{WaveFrame,AbstractVector}; ratio::Real = 0.05)
    return taper(var -> var, x; ratio = ratio)
end

function bandpass(x::AbstractVector, w1::AbstractFloat, w2::AbstractFloat, n::Int = 4;
                  fs::Real = 0.0)
    if fs > 0.0
        ftr = digitalfilter(Bandpass(w1, w2; fs = fs), Butterworth(n))
    else
        ftr = digitalfilter(Bandpass(w1, w2; fs = fs), Butterworth(n))
    end
    return filtfilt(ftr, x)
end

function bandpass!(x::AbstractVector, w1::AbstractFloat, w2::AbstractFloat, n::Int = 4;
                   fs::Real = 0.0)
    y = bandpass(x, w1, w2, n; fs = fs)
    for i = 1:length(x)
        x[i] = y[i]
    end
    return nothing
end

function bandpass(frame::WaveFrame, w1::AbstractFloat, w2::AbstractFloat, n::Int = 4;
                  fs::Real = 0.0)
    y = deepcopy(frame)
    yw = bandpass(frame.data, w1, w2, n; fs = fs)
    for i = 1:length(yw)
        y.data[i] = yw[i]
    end
    return y
end

function calshift(t::DateTime, ref::DateTime, dt::Period)
    return round(Int, round(t - ref, Nanosecond) / round(dt, Nanosecond))
end

function merge(f::Vector{T}, gv...; fill::Real = 0.0) where {T<:Frame}
    @assert eltype(gv) <: Frame
    flag = (delta = true, id = true)
    for ig = 1:length(gv), k in keys(flag)
        flag[k] &= getproperty(f, k) == getproperty(gv[ig], k)
    end
    for k in keys(flag)
        if !flag[k]
            @error "property: " * string(k) * " is not consistant"
        end
    end
    btimes = DateTime[f.begintime]
    etimes = DateTime[f.begintime+Microsecond(f.delta * (f.npts - 1))]
    newBtime = minimum(btimes)
    newEtime = maximum(etimes)
    newDelta = f.delta
    newNPTS = round(Int, round(newEtime - newBtime, Microsecond) / round(newDelta, Microsecond))
    newdata = fill(convert(eltype(f.data), fill), newNPTS)
    s = calshift(f.begintime, newBtime, newDelta)
    for i = 1:f.npts
        newdata[s+i] = f.data[i]
    end
    for ig = 1:length(gv)
        s = calshift(gv[ig], newBtime, newDelta)
        for i = 1:gv[ig].npts
            newdata[s+i] = gv[ig].data[i]
        end
    end
    return WaveFrame(f.network, f.station, f.device, f.component, newBtime, newDelta, newdata,
                     newNPTS, f.meta)
end

function slice(frame::WaveFrame, starttime::DateTime, endtime::DateTime; fill::Real=0.0)
end
