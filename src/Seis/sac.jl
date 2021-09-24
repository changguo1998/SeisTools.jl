module SAC
function readhead(io::IO) end
function read(io::IO) end
function read(path::AbstractString)
    return open(readsac, path)
end
function write(io::IO, frame::WaveFrame) end
function write(path::AbstractString, frame::WaveFrame)
    open(path, "w") do io
        return write(io, frame)
    end
    return nothing
end
function emptyheader()
end

function standardname(frame::WaveFrame; standard::AbstractString="iris")
    if standard == ""
        return @sprintf("%04d.%03d")
    else
    end
end
end
