module SeisTools

# * basic macros and functions

"""
    @must(cond, text = "")

like @assert, but is user defined and will always execute.
if `cond` is false, the macro will throw an error with `text`
"""
macro must(cond, text = "")
    return :(if !($(esc(cond)))
                 error($(esc(text)))
             end)
end

"""
    @hadbetter(cond, text = "")

warning when `cond` is false with information `text`
"""
macro hadbetter(cond, text = "")
    return :(if !($(esc(cond)))
                 @warn($text)
             end)
end

# * modules
include("SAC.jl")
include("SEGY.jl")
include("SSP.jl")
include("SACPZ.jl")
include("Geodesy.jl")

end # module
