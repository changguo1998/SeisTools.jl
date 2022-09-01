module SACPZ
using Dates, Printf
using Base: parse as baseparse

HEAD = ("NETWORK", "STATION", "LOCATION", "CHANNEL", "CREATED", "START", "END", "DESCRIPTION", "LATITUDE", "LONGITUDE",
        "ELEVATION", "DEPTH", "DIP", "AZIMUTH", "SAMPLE_RATE", "INPUT_UNIT", "OUTPUT_UNIT", "INSTTYPE", "INSTGAIN",
        "COMMENT", "SENSITIVITY", "A0")
FMT = ("* NETWORK   (KNETWK): %s\n", "* STATION    (KSTNM): %s\n", "* LOCATION   (KHOLE): %s\n",
       "* CHANNEL   (KCMPNM): %s\n", "* CREATED           : %04d-%02d-%02dT%02d:%02d:%02d\n",
       "* START             : %04d-%02d-%02dT%02d:%02d:%02d\n", "* END               : %04d-%02d-%02dT%02d:%02d:%02d\n",
       "* DESCRIPTION       :%s\n", "* LATITUDE          : %.6f\n", "* LONGITUDE         : %.6f\n",
       "* ELEVATION         : %.1f\n", "* DEPTH             : %.1f\n", "* DIP               : %.1f\n",
       "* AZIMUTH           : %.1f\n", "* SAMPLE RATE       : %.1f\n", "* INPUT UNIT        : %s\n",
       "* OUTPUT UNIT       : %s\n", "* INSTTYPE          : %s\n", "* INSTGAIN          : %.6e %s\n",
       "* COMMENT           : %s\n", "* SENSITIVITY       : %.6e %s\n", "* A0                : %.6e\n")

function parse_one_station(lines::Vector{T}) where {T<:AbstractString}
    d = Dict{String,Any}()
    # header
    for h in HEAD
        for l in lines
            if contains(l, replace(h, "_" => " "))
                idx = findfirst(x -> x == ':', l)
                if h in
                   ("NETWORK", "STATION", "LOCATION", "CHANNEL", "DESCRIPTION", "INPUT_UNIT", "OUTPUT_UNIT", "INSTTYPE",
                    "COMMENT") # string
                    d[h] = String(strip(l[idx+1:end]))
                elseif h in ("CREATED", "START", "END") # datetime
                    d[h] = DateTime(strip(l[idx+1:end]))
                elseif h in ("LATITUDE", "LONGITUDE", "ELEVATION", "DEPTH", "DIP", "AZIMUTH", "SAMPLE_RATE", "A0") # float
                    d[h] = baseparse(Float64, strip(l[idx+1:end]))
                elseif h in ("INSTGAIN", "SENSITIVITY")
                    tl = split(strip(l); keepempty = false)
                    d[h] = [baseparse(Float64, tl[end-1]), tl[end]]
                end
            end
        end
        if !(h in keys(d))
            # @warn "key $h not found. (line $(findall(x->x==l, lines)))"
            d[h] = "NOT FOUND"
        end
        # println("- $h => $(d[h])")
    end
    # zeros
    idx = findfirst(x -> contains(x, "ZEROS"), lines)
    nz = baseparse(Int, strip(lines[idx][6:end]))
    ZS = zeros(ComplexF64, nz)
    for i = 1:nz
        t = split(strip(lines[idx+i]); keepempty = false)
        ZS[i] = ComplexF64(baseparse(Float64, t[1]), baseparse(Float64, t[2]))
    end
    d["ZEROS"] = ZS
    # poles
    idx = findfirst(x -> contains(x, "POLES"), lines)
    np = baseparse(Int, strip(lines[idx][6:end]))
    PS = zeros(ComplexF64, np)
    for i = 1:np
        t = split(strip(lines[idx+i]); keepempty = false)
        PS[i] = ComplexF64(baseparse(Float64, t[1]), baseparse(Float64, t[2]))
    end
    idx = findfirst(x -> contains(x, "CONSTANT"), lines)
    d["POLES"] = PS
    # constant
    C = baseparse(Float64, strip(lines[idx][9:end]))
    d["CONSTANT"] = C
    return d
end

function parse_combine(lines::Vector{T}) where {T<:AbstractString}
    el = findall(x -> contains(x, "CONSTANT"), lines)
    if length(el) > 1
        bl = [1; el[1:end-1] .+ 1]
        return map((x, y) -> parse_one_station(lines[x:y]), bl, el)
    else
        return parse_one_station(lines)
    end
end

parse(lines::Vector{T}) where {T<:AbstractString} = parse_combine(lines)
parsefile(path::String) = readlines(path) |> parse_combine

function pzfmt2dict(d::Dict)
    nd = Dict{String,Any}()
    for h in PZFile.HEAD
        nd[h] = (typeof(d[h]) <: Vector) ? join(d[h], ' ') : d[h]
    end
    nd["CONSTANT"] = d["CONSTANT"]
    nd["POLES_RE"] = real.(d["POLES"])
    nd["POLES_IM"] = imag.(d["POLES"])
    nd["ZEROS_RE"] = real.(d["ZEROS"])
    nd["ZEROS_IM"] = imag.(d["ZEROS"])
    return nd
end

function print(io::IO, x::Dict{String,Any})
    Printf.@printf(io, "* **********************************\n")
    for ih = 1:length(HEAD)
        h = HEAD[ih]
        fmt = FMT[ih]
        if h ∈ keys(x)
            if h in ("CREATED", "START", "END")
                s = Printf.format(io, Printf.Format(fmt), (x[h] .|> [year, month, day, hour, minute, second])...)
            elseif h in ("INSTGAIN", "SENSITIVITY")
                s = Printf.format(io, Printf.Format(fmt), x[h]...)
            else
                s = Printf.format(io, Printf.Format(fmt), x[h])
            end
        end
    end
    Printf.@printf(io, "* **********************************\n")
    Printf.@printf(io, "ZEROS\t%d\n", length(x["ZEROS"]))
    for i = 1:length(x["ZEROS"])
        Printf.@printf(io, "\t%+.6e\t%+.6e\n", real(x["ZEROS"][i]), imag(x["ZEROS"][i]))
    end
    Printf.@printf(io, "POLES\t%d\n", length(x["POLES"]))
    for i = 1:length(x["POLES"])
        Printf.@printf(io, "\t%+.6e\t%+.6e\n", real(x["POLES"][i]), imag(x["POLES"][i]))
    end
    Printf.@printf(io, "CONSTANT\t%+.6e\n", x["CONSTANT"])
end
end
