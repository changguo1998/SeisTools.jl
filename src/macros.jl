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

macro linearscale(x, x1, x2, y1, y2)
    return :(($(esc(x)) - $(esc(x1))) / ($(esc(x2)) - $(esc(x1))) * ($(esc(y2)) - $(esc(y1))) + $(esc(y1)))
end
