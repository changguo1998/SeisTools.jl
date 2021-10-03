module SAC

using Dates, Printf, Statistics

sacheadlist = ("delta", "depmin", "depmax", "scale", "odelta", "b", "e", "o", "a", "internal1",
               "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "f", "resp0", "resp1",
               "resp2", "resp3", "resp4", "resp5", "resp6", "resp7", "resp8", "resp9", "stla",
               "stlo", "stel", "stdp", "evla", "evlo", "evel", "evdp", "mag", "user0", "user1",
               "user2", "user3", "user4", "user5", "user6", "user7", "user8", "user9", "dist", "az",
               "baz", "gcarc", "internal2", "internal3", "depmen", "cmpaz", "cmpinc", "xminimun",
               "xmaximum", "yminimum", "ymaximum", "unused1", "unused2", "unused3", "unused4",
               "unused5", "unused6", "unused7", "nzyear", "nzjday", "nzhour", "nzmin", "nzsec",
               "nzmsec", "nvhdr", "norid", "nevid", "npts", "internal4", "nwfid", "nxsize",
               "nysize", "unused8", "iftype", "idep", "iztype", "unused9", "iinst", "istreg",
               "ievreg", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc", "unused10", "unused11",
               "unused12", "unused13", "unused14", "unused15", "unused16", "unused17", "leven",
               "lpspol", "lovrok", "lcalda", "unused18", "kstnm", "kevnm", "khole", "ko", "ka",
               "kt0", "kt1", "kt2", "kt3", "kt4", "kt5", "kt6", "kt7", "kt8", "kt9", "kf", "kuser0",
               "kuser1", "kuser2", "kcmpnm", "knetwk", "kdatrd", "kinst")

specialdefaultvalue = (nvhdr = 6, iftype = "ITIME", idep = "IVEL", iztype = "IB", ievtyp = "IUNKN",
                       iqual = "IOTHER", isynth = "IRLDA", imagsrc = "IUNKNOWN", leven = true,
                       lpspol = false, lovrok = true, lcalda = false)

enumeratevars = ("iftype", "idep", "iztype", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc")

logicalvars = ("leven", "lpspol", "lovrok", "lcalda", "unused18")

ver7foot = ("delta", "b", "e", "o", "a", "t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9",
            "f", "evlo", "evla", "stlo", "stla", "sb", "sdelta")

autofillvars = ("depmin", "depmax", "depmen", "dist", "az", "baz", "gcarc")

indextable = Dict([(1, "ITIME"), (2, "IRLIM"), (3, "IAMPH"), (4, "IXY"), (5, "IUNKN"), (6, "IDISP"),
                   (7, "IVEL"), (8, "IACC"), (9, "IB"), (10, "IDAY"), (11, "IO"), (12, "IA"),
                   (13, "IT0"), (14, "IT1"), (15, "IT2"), (16, "IT3"), (17, "IT4"), (18, "IT5"),
                   (19, "IT6"), (20, "IT7"), (21, "IT8"), (22, "IT9"), (37, "INUCL"), (38, "IPREN"),
                   (39, "IPOSTN"), (40, "IQUAKE"), (41, "IPREQ"), (42, "IPOSTQ"), (43, "ICHEM"),
                   (44, "IOTHER"), (45, "IGOOD"), (46, "IGLCH"), (47, "IDROP"), (48, "ILOSN"),
                   (49, "IRLDA"), (50, "IVOLTS"), (51, "IXYZ"), (52, "IMB"), (53, "IMS"),
                   (54, "IML"), (55, "IMW"), (56, "IMD"), (57, "IMX"), (58, "INEIC"), (59, "IPDE"),
                   (60, "IISC"), (61, "IREB"), (62, "IUSGS"), (63, "IBRK"), (64, "ICALTECH"),
                   (65, "ILLNL"), (66, "IEVLOC"), (67, "IJSOP"), (68, "IUSER"), (69, "IUNKNOWN"),
                   (70, "IQB"), (71, "IQB1"), (72, "IQB2"), (73, "IQBX"), (74, "IQMT"), (75, "IEQ"),
                   (76, "IEQ1"), (77, "IEQ2"), (78, "IME"), (79, "IEX"), (80, "INU"), (81, "INC"),
                   (82, "IO_"), (83, "IL"), (84, "IR"), (85, "IT"), (86, "IU")])

hashtable = Dict([(i.second, i.first) for i in colle(indextable)])

simplify(x::Vector{UInt8}) = strip(String(Char.(x)), [' ', '\0'])

"""
readhead(io::IO)

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
    hi = Int64.(hi)
    hi[hi.==-12345] .= NaN
    htu = Vector{Tuple{String,Any}}(undef, 0)
    for i = 1:70
        push!(htu, (varname[i], hf[i]))
    end
    for i = 1:40
        push!(htu, (varname[70+i], hi[i]))
    end
    push!(htu, (varname[111], simplify(hc[1:8])))
    push!(htu, (varname[112], simplify(hc[9:24])))
    for i = 3:23
        push!(htu; (varname[110+i], simplify(hc[(i*8+1):(i*8+8)])))
    end
    head = Dict(htu)
    for i in enumeratevars
        if !isnan(head[i])
            head[i] = indextable[head[i]]
        end
    end
    for i in logicalvars
        if !isnan(head[i])
            head[i] = Bool(head[i])
        end
    end
    if head["nvhdr"] == 7
        foot = zeros(Float64, 22)
        read!(io, foot)
        for i = 1:length(ver7foot)
            if foot[i] != -12345.0
                head[ver7foot[i]] = foot[i]
            end
        end
    end
    return head
end

"""
read(io::IO) -> frame::WaveFrame

read sac data from opened io
"""
function read(io::IO)
    head = readhead(io)
    td = zeros(Float32, Int(head["npts"]))
    if head["npts"] >= 1 &&
       (head["iftype"] in ["IXY", "IRLIM", "IAMPH", "IXYZ", "ITIME"]) &&
       !isnan(head["leven"])
        read!(io, td)
        if !head["leven"] || (head["iftype"] in ["IAMPH", "IRLIM"])
            td1 = Float64.(td)
            read!(io, td)
            data = [td1, Float64.(td)]
        else
            data = Float64.(td)
        end
    else
        data = Float64[]
    end
    return WaveFrame(head["knetwk"], head["kstnm"], head["khole"], head["kcmpnm"][end],
                     DateTime(head["nzyear"], 1, 1, head["nzhour"], head["nzmin"], head["nzsec"],
                              head["nzmsec"]) + Day(head["nzjday"] - 1), head["delta"], data, head)
end

"""
read(path::String) -> frame::WaveFrame

read sac data from opened io
"""
function read(path::AbstractString)
    return open(readsac, path)
end

function autocal!(frame::WaveFrame)
    frame.meta["depmax"] = maximum(meta.data)
    frame.meta["depmen"] = mean(meta.data)
    frame.meta["depmim"] = minimum(meta.data)
    return nothing
end

"""
write(io::IO, frame::WaveFrame; autocal::Bool=true)

write WaveFrame to file with SAC format. If `autocal` is true, the header variable:
`depmax`, `depmen`, `depmin` will be update before writting.
"""
function write(io::IO, frame::WaveFrame; autocalc::Bool = true)
    tf = deepcopy(frame)
    if autocalc
        autocal!(tf)
    end
    for i = 1:110
        hname = sacheadlist[i]
        if isequal(tf.meta[hname], NaN)
            if i <= 70
                write(io, Float32(-12345.0))
            else
                write(io, Int32(-12345))
            end
        else
            if i <= 70
                write(io, Float32(tf.meta[hname]))
            else
                if hname in enumeratevars
                    write(io, Int32(hashtable[tf.meta[hname]]))
                else
                    write(io, Int32(tf.meta[hname]))
                end
            end
        end
    end
    hname = sacheadlist[111]
    v = tf.meta[hname]
    for i = 1:8
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    hname = sacheadlist[112]
    v = tf.meta[hname]
    for i = 1:16
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    for i = 113:length(sacheadlist)
        hname = sacheadlist[i]
        v = tf.meta[hname]
        for j = 1:8
            if j <= length(v)
                write(io, UInt8(v[j]))
            else
                write(io, UInt8(' '))
            end
        end
    end
    for i in tf.data
        write(io, Float32.(i))
    end
    return nothing
end

"""
write(io::IO, frame::WaveFrame; autocal::Bool=true)
"""
function write(path::AbstractString, frame::WaveFrame; autocalc::Bool = true)
    open(path, "w") do io
        write(io, frame; autocalc = autocalc)
    end
    return nothing
end

"""
emptyheader(; hdrvars...)

return a Dict filled with sac header variables. Default value is -12345, and user can also values.
"""
function emptyheader(; hdrvars...)
    th = Dict{String,Any}()
    setvars = string.(keys(specialdefaultvalue))
    for i = 1:length(sacheadlist)
        var = sacheadlist[i]
        if var in setvars
            th[var] = specialdefaultvalue[Symbol(var)]
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
        if string(var) in sacheadlist
            th[string(var)] = hdrvars[var]
        end
    end
    return th
end

"""
standardname(frame::WaveFrame; standard::AbstractString = "iris",
                      quality::AbstractString = "D", kwargs...)
"""
function standardname(frame::WaveFrame; standard::AbstractString = "iris",
                      quality::AbstractString = "D", kwargs...)
    if standard == "iris"
        return @sprintf("%s.%s.%s.%s.%s.%04d.%02d.%02d.%02d.%02d.%02d.%03d.SAC",
                        frame.meta["knetwk"], frame.meta["kstnm"], frame.meta["khole"],
                        frame.meta["kcmpnm"], quality, year(frame.begintime),
                        month(frame.begintime), day(frame.begintime), hour(frame.begintime),
                        minute(frame.begintime), second(frame.begintime),
                        millisecond(frame.begintime))
    else
        return @sprintf("%04d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.%s.SAC", year(frame.begintime),
                        dayofyear(frame.begintime), hour(frame.begintime), minute(frame.begintime),
                        second(frame.begintime), millisecond(frame.begintime), frame.meta["knetwk"],
                        frame.meta["kstnm"], frame.meta["khole"], frame.meta["kcmpnm"], quality)
    end
end
end
