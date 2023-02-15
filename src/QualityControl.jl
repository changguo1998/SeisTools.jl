module QualityControl

using Dates, Statistics

function maxapproxconst(x::AbstractVector)
    flag = falses(length(x) - 1)
    Threads.@threads for i = 1:length(x)-1
        flag[i] = abs(x[i]-x[i+1]) < max(abs(x[i]), abs(x[i+1])) * 0.05
    end
    nconst = zeros(Int, length(x))
    nconst[1] = 0
    for i = 2:length(x)
        if flag[i-1]
            nconst[i] = nconst[i-1] + 1
        else
            nconst[i] = 0
        end
    end
    return nconst
end

"""
```
constrecord(x::AbstractVector, window::Integer) -> Bool
```

return true if there is constant value in waveform
"""
function hasconstrecord(x::AbstractVector, window::Integer)
    nc = maxapproxconst(x)
    return maximum(nc) >= window - 1
end

"""
```
constrecord(x::AbstractVector, dt::Period, window::Period) -> Bool
```

return true if there is constant value in waveform
"""
hasconstrecord(x::AbstractVector, dt::Period, window::Period) = hasconstrecord(x,
                                                                               round(Int,
                                                                                     Millisecond(window) /
                                                                                     Millisecond(dt)))

"""
```
limitedamplitude(x::AbstractVector, hw::Integer)
```

test if the waveform's amplitude is limited
"""
function limitedamplitude(x::AbstractVector, hw::Integer)
    flag = falses(length(x))
    vmax = maximum(x)
    vmin = minimum(x)
    sigma0 = std(x; corrected=true)
    Threads.@threads for i = eachindex(x)
        sigma = std(x[max(1, i-hw):min(L, i+hw)]; corrected=true)
        flag[i] = (sigma/sigma0) < 0.5
    end
    nconst = zeros(Int, length(x))
    nconst[1] = 0
    for i = 2:length(x)
        if flag[i-1]
            nconst[i] = nconst[i-1] + 1
        else
            nconst[i] = 0
        end
    end
    return maximum(nconst) >= hw - 1
end

function kurtosis(x::AbstractVector)
    mx = mean(x)
    sx = var(x; corrected = false)
    return sum(v -> (v - mx)^4, x) / sx^2 / length(x)
end

function kurtosis_sample(x::AbstractVector)
    l = length(x)
    return kurtosis(x) * (l + 1) / (l - 3) - 3 * (l - 1)^2 / (l - 2) / (l - 3)
end

end
