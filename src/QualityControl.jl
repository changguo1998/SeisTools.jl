module QualityControl

using Dates, Statistics

function maxconst(x::AbstractVector)
    flag = falses(length(x) - 1)
    Threads.@threads for i = 1:length(x)-1
        flag[i] = x[i] == x[i+1]
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
    nc = maxconst(x)
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
