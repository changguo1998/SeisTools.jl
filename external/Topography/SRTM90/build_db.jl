using DelimitedFiles

# include(abspath(@__DIR__, "../../../src/basic.jl"))
include( "../../../src/basic.jl")

function _get_par(io::IO, T::Type)
    t = readline(io)
    l = split(t, ' ', keepempty=false)
    return parse(T, l[2])
end

function _read_asc_header(fn::String)
    io = open(fn)
    ncol = _get_par(io, Int)
    nrow = _get_par(io, Int)
    olon = _get_par(io, Float64)
    olat = _get_par(io, Float64)
    dd = _get_par(io, Float64)
    nanvalue = _get_par(io, Float64)
    close(io)
    return (nrow, ncol, olat, olon, dd, nanvalue)
end

mkpath("buffer")

zip_files = readdir("raw_data")

for zf in zip_files
    local ascf = replace(zf, ".zip"=>".asc")
    if isfile(joinpath("buffer", ascf))
        continue
    end
    cmd = Cmd(["unzip", "-o", joinpath("raw_data", zf), "-d", "buffer"])
    run(cmd)
end

# exit(0)

ascfiles = filter(endswith(".asc"), readdir("buffer"))
for af in ascfiles
    @info "Build block database: $af"
    local splitnm = splitext(af)
    local newname = "block_" * splitnm[1] * ".bin"
    if isfile(newname)
        continue
    end
    local nrow, ncol, olat, olon, dd, nanvalue
    (nrow, ncol, olat, olon, dd, nanvalue) = _read_asc_header(joinpath("buffer", af))
    local lon = range(olon, olon+5.0, ncol)
    local lat = range(olat, olat+5.0, nrow)
    local m = readdlm(joinpath("buffer", af); skipstart=6)
    reverse!(m, dims=1)
    m = permutedims(m)
    m[m.==nanvalue] .= NaN
    local io = open(newname, "w")
    _write_f64_vector!(io, collect(lat))
    _write_f64_vector!(io, collect(lon))
    _write_f64_array(io, m)
    close(io)
end
