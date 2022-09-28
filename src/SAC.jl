module SAC

using Dates, Printf, Statistics

include("macros.jl")

SAC_HEADER_LIST = ("delta", "depmin", "depmax", "scale", "odelta", "b", "e", "o", "a", "internal1", "t0", "t1", "t2",
                   "t3", "t4", "t5", "t6", "t7", "t8", "t9", "f", "resp0", "resp1", "resp2", "resp3", "resp4", "resp5",
                   "resp6", "resp7", "resp8", "resp9", "stla", "stlo", "stel", "stdp", "evla", "evlo", "evel", "evdp",
                   "mag", "user0", "user1", "user2", "user3", "user4", "user5", "user6", "user7", "user8", "user9",
                   "dist", "az", "baz", "gcarc", "internal2", "internal3", "depmen", "cmpaz", "cmpinc", "xminimun",
                   "xmaximum", "yminimum", "ymaximum", "unused1", "unused2", "unused3", "unused4", "unused5", "unused6",
                   "unused7", "nzyear", "nzjday", "nzhour", "nzmin", "nzsec", "nzmsec", "nvhdr", "norid", "nevid",
                   "npts", "internal4", "nwfid", "nxsize", "nysize", "unused8", "iftype", "idep", "iztype", "unused9",
                   "iinst", "istreg", "ievreg", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc", "unused10",
                   "unused11", "unused12", "unused13", "unused14", "unused15", "unused16", "unused17", "leven",
                   "lpspol", "lovrok", "lcalda", "unused18", "kstnm", "kevnm", "khole", "ko", "ka", "kt0", "kt1", "kt2",
                   "kt3", "kt4", "kt5", "kt6", "kt7", "kt8", "kt9", "kf", "kuser0", "kuser1", "kuser2", "kcmpnm",
                   "knetwk", "kdatrd", "kinst")

SPECIAL_DEFAULT_VALUE = (nvhdr = 6, iftype = "ITIME", idep = "IVEL", iztype = "IB", ievtyp = "IUNKN", iqual = "IOTHER",
                         isynth = "IRLDA", imagsrc = "IUNKNOWN", leven = true, lpspol = false, lovrok = true,
                         lcalda = false)

ENUMERATE_VAR_LIST = ("iftype", "idep", "iztype", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc")

LOGICAL_VAR_LIST = ("leven", "lpspol", "lovrok", "lcalda", "unused18")

VERSION_7_FOOT_VAR = ("delta", "b", "e", "o", "a", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "f",
                      "evlo", "evla", "stlo", "stla", "sb", "sdelta")

AUTO_CALC_VAR_LIST = ("depmin", "depmax", "depmen", "dist", "az", "baz", "gcarc")

indextable = Dict([(1, "ITIME"), (2, "IRLIM"), (3, "IAMPH"), (4, "IXY"), (5, "IUNKN"), (6, "IDISP"), (7, "IVEL"),
                   (8, "IACC"), (9, "IB"), (10, "IDAY"), (11, "IO"), (12, "IA"), (13, "IT0"), (14, "IT1"), (15, "IT2"),
                   (16, "IT3"), (17, "IT4"), (18, "IT5"), (19, "IT6"), (20, "IT7"), (21, "IT8"), (22, "IT9"),
                   (37, "INUCL"), (38, "IPREN"), (39, "IPOSTN"), (40, "IQUAKE"), (41, "IPREQ"), (42, "IPOSTQ"),
                   (43, "ICHEM"), (44, "IOTHER"), (45, "IGOOD"), (46, "IGLCH"), (47, "IDROP"), (48, "ILOSN"),
                   (49, "IRLDA"), (50, "IVOLTS"), (51, "IXYZ"), (52, "IMB"), (53, "IMS"), (54, "IML"), (55, "IMW"),
                   (56, "IMD"), (57, "IMX"), (58, "INEIC"), (59, "IPDE"), (60, "IISC"), (61, "IREB"), (62, "IUSGS"),
                   (63, "IBRK"), (64, "ICALTECH"), (65, "ILLNL"), (66, "IEVLOC"), (67, "IJSOP"), (68, "IUSER"),
                   (69, "IUNKNOWN"), (70, "IQB"), (71, "IQB1"), (72, "IQB2"), (73, "IQBX"), (74, "IQMT"), (75, "IEQ"),
                   (76, "IEQ1"), (77, "IEQ2"), (78, "IME"), (79, "IEX"), (80, "INU"), (81, "INC"), (82, "IO_"),
                   (83, "IL"), (84, "IR"), (85, "IT"), (86, "IU")])

hashtable = Dict([(i.second, i.first) for i in collect(indextable)])

function simplify(x::Vector{UInt8})
    idx = findfirst(==(0x00), x)
    if isnothing(idx)
        return String(strip(String(Char.(x)), [' ', '\0']))
    else
        return String(strip(String(Char.(x[1:idx-1])), [' ', '\0']))
    end
end

function DateTime(meta::Dict)
    return Dates.DateTime(meta["nzyear"], 1, 1, meta["nzhour"], meta["nzmin"], meta["nzsec"], meta["nzmsec"]) +
           Day(meta["nzjday"] - 1)
end

"""
readhead(io::IO) -> Dict

read header from opened io
"""
function readhead(io::IO)
    hf = zeros(Float32, 70)
    hi = zeros(Int32, 40)
    hc = zeros(UInt8, 192)
    read!(io, hf)
    read!(io, hi)
    read!(io, hc)
    hf = Float64.(hf)
    hf[hf.==-12345.0] .= NaN
    hi = Float64.(hi)
    hi[hi.==-12345.0] .= NaN
    htu = Vector{Tuple{String,Any}}(undef, 0)
    for i = 1:70
        push!(htu, (SAC_HEADER_LIST[i], hf[i]))
    end
    for i = 1:40
        push!(htu, (SAC_HEADER_LIST[70+i], isnan(hi[i]) ? NaN : Int(hi[i])))
    end
    push!(htu, (SAC_HEADER_LIST[111], simplify(hc[1:8])))
    push!(htu, (SAC_HEADER_LIST[112], simplify(hc[9:24])))
    for i = 3:23
        push!(htu, (SAC_HEADER_LIST[110+i], simplify(hc[(i*8+1):(i*8+8)])))
    end
    head = Dict(htu)
    for i in ENUMERATE_VAR_LIST
        if !isnan(head[i])
            head[i] = indextable[head[i]]
        else
            head[i] = "IUNKN"
        end
    end
    for i in LOGICAL_VAR_LIST
        if !isnan(head[i])
            head[i] = Bool(head[i])
        end
    end
    if head["nvhdr"] == 7
        @debug "SAC version: 7"
        foot = zeros(Float64, 22)
        read!(io, foot)
        for i in eachindex(VERSION_7_FOOT_VAR)
            if foot[i] != -12345.0
                head[VERSION_7_FOOT_VAR[i]] = foot[i]
            end
        end
    end
    return head
end

"""
read(io::IO) -> (hdr=head, data=data)

read sac data from opened io
"""
function read(io::IO)
    head = readhead(io)
    td = zeros(Float32, Int(head["npts"]))
    if head["npts"] < 1
        @error "header: npts less than 1"
        return (hdr = head, data = Float64[])
    end
    if isnan(head["leven"])
        @error "header: leven is NaN"
        return (hdr = head, data = Float64[])
    end
    if head["iftype"] ∉ ("IXY", "IRLIM", "IAMPH", "IXYZ", "ITIME")
        @error "header: iftype is illegal"
        return (hdr = head, data = Float64[])
    end
    @debug "Data type: $(head["iftype"])"
    if head["iftype"] == "ITIME"
        read!(io, td)
        return (hdr = head, data = td)
    elseif head["iftype"] in ("IXY", "IRLIM", "IAMPH")
        read!(io, td)
        td1 = deepcopy(td)
        read!(io, td)
        return (hdr = head, data = (td1, td))
    else
        read!(io, td)
        td1 = deepcopy(td)
        read!(io, td)
        td2 = deepcopy(td)
        read!(io, td)
        return (hdr = head, data = (td1, td2, td))
    end
end

"""
read(path::String) -> (hdr=head, data=data)

read sac data from opened io
"""
function read(path::AbstractString)
    return open(read, path, "r")
end

function autocal!(hdr, data)
    if hdr["iftype"] == "ITIME"
        hdr["depmax"] = maximum(data)
        hdr["depmen"] = mean(data)
        hdr["depmin"] = minimum(data)
    end
    return nothing
end

"""
write(io::IO, hdr::Dict, data; autocal::Bool=true)

write WaveFrame to file with SAC format. If `autocal` is true, the header variable:
`depmax`, `depmen`, `depmin` will be update before writting.
"""
function write!(io::IO, hdr::Dict, data; autocalc::Bool = true)
    if autocalc
        autocal!(hdr, data)
    end
    for i = 1:110
        hname = SAC_HEADER_LIST[i]
        if isequal(hdr[hname], NaN)
            if i <= 70
                Base.write(io, Float32(-12345.0))
            else
                Base.write(io, Int32(-12345))
            end
        else
            if i <= 70
                Base.write(io, Float32(hdr[hname]))
            else
                if hname in ENUMERATE_VAR_LIST
                    if typeof(hdr[hname]) <: AbstractString
                        Base.write(io, Int32(hashtable[hdr[hname]]))
                    else
                        Base.write(io, Int32(-12345))
                    end
                else
                    Base.write(io, Int32(hdr[hname]))
                end
            end
        end
    end
    hname = SAC_HEADER_LIST[111]
    v = hdr[hname]
    for i = 1:8
        if i <= length(v)
            Base.write(io, UInt8(v[i]))
        else
            Base.write(io, UInt8(' '))
        end
    end
    hname = SAC_HEADER_LIST[112]
    v = hdr[hname]
    for i = 1:16
        if i <= length(v)
            Base.write(io, UInt8(v[i]))
        else
            Base.write(io, UInt8(' '))
        end
    end
    for i = 113:length(SAC_HEADER_LIST)
        hname = SAC_HEADER_LIST[i]
        v = hdr[hname]
        for j = 1:8
            if j <= length(v)
                Base.write(io, UInt8(v[j]))
            else
                Base.write(io, UInt8(' '))
            end
        end
    end
    if hdr["iftype"] == "ITIME"
        Base.write(io, Float32.(data))
    else
        for i in data
            Base.write(io, Float32.(i))
        end
    end
    return nothing
end

function write(io::IO, hdr::Dict, data; autocalc::Bool = true)
    thdr = deepcopy(hdr)
    tdata = deepcopy(data)
    write!(io, thdr, tdata; autocalc = autocalc)
    return nothing
end

"""
write(io::IO, frame::WaveFrame; autocal::Bool=true)
"""
function write(path::AbstractString, hdr::Dict, data; autocalc::Bool = true)
    open(path, "w") do io
        write(io, hdr, data; autocalc = autocalc)
    end
    return nothing
end

"""
emptyheader(; hdrvars...)

return a Dict filled with sac header variables. Default value is -12345, and user can also values.
"""
function emptyheader(; hdrvars...)
    th = Dict{String,Any}()
    setvars = string.(keys(SPECIAL_DEFAULT_VALUE))
    for i in eachindex(SAC_HEADER_LIST)
        var = SAC_HEADER_LIST[i]
        if var in setvars
            th[var] = SPECIAL_DEFAULT_VALUE[Symbol(var)]
        else
            if i <= 70
                th[var] = -12345.0
            elseif i <= 110
                th[var] = -12345
            else
                th[var] = "-12345"
            end
        end
    end
    for var in keys(hdrvars)
        if string(var) in SAC_HEADER_LIST
            th[string(var)] = hdrvars[var]
        end
    end
    return th
end

"""
standardname(frame::WaveFrame; standard::AbstractString = "iris",
quality::AbstractString = "D", kwargs...)
"""
function standardname(hdr::Dict; standard::AbstractString = "iris", quality::AbstractString = "D", kwargs...)
    if standard == "iris"
        return @sprintf("%s.%s.%s.%s.%s.%04d.%02d.%02d.%02d.%02d.%02d.%03d.SAC", hdr["knetwk"], hdr["kstnm"],
                        hdr["khole"], hdr["kcmpnm"], quality, hdr["nzyear"], month(DateTime(hdr)), day(DateTime(hdr)),
                        hdr["nzhour"], hdr["nzmin"], hdr["nzsec"], hdr["nzmsec"])
    else
        return @sprintf("%04d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.%s.SAC", hdr["nzyear"], hdr["nzjday"], hdr["nzhour"],
                        hdr["nzmin"], hdr["nzsec"], hdr["nzmsec"], hdr["knetwk"], hdr["kstnm"], hdr["khole"],
                        hdr["kcmpnm"], quality)
    end
end
end
