module QualityControl

using Dates, Statistics

"""
```
constrecord(x::AbstractVector, window::Integer) -> Bool
```

return true if there is constant value in waveform
"""
function constrecord(x::AbstractVector, window::Integer)
    flag = falses(length(x)-1)
    for i = 1:length(x)-1
        flag[i] = x[i] == x[i+1]
    end
    fr = false
    for i = 1:length(flag)-window+2
        ft = true
        for j = 1:window-1
            ft &= flag[i+j-1]
        end
        fr |= ft
        if fr
            break
        end
    end
    return fr
end

"""
```
constrecord(x::AbstractVector, dt::Period, window::Period) -> Bool
```

return true if there is constant value in waveform
"""
constrecord(x::AbstractVector, dt::Period, window::Period) = constrecord(x, round(Int, Millisecond(window)/Millisecond(dt)))

function kurtosis(x::AbstractVector)
    mx = mean(x)
    sx = var(x, corrected=false)
    return sum(v->(v-mx)^4, x)/sx^2/length(x)
end

function kurtosis_sample(x::AbstractVector)
    l = length(x)
    return kurtosis(x)*(l+1)/(l-3) - 3*(l-1)^2/(l-2)/(l-3)
end

end