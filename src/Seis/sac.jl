# sac related functions

"""
readsachead( path::IO )
=====================================
read head of sac binary file
"""
function readsachead(io::IO)
    hf = zeros(Float32, 70)
    hi = zeros(Int32, 40)
    hc = zeros(UInt8, 192)
    read!(io, hf)
    read!(io, hi)
    read!(io, hc)
    # change undefined
    htu = [];
    hf = Float64.(hf)
    hf[hf .== -12345.0] .= NaN
    hi = Float64.(hi)
    hi[hi .== -12345.0] .= NaN
    varname = HEADER.sac.headlist
    # construct head
    for i = 1:70
        htu = [htu; (varname[i], hf[i])]
    end
    for i = 1:40
        htu = [htu; (varname[70 + i], hi[i])]
    end
    htu = [htu; (varname[111], Char.(hc[1:8]) |> String |> x -> strip(x, [' ', '\0']))]
    htu = [htu; (varname[112], Char.(hc[9:24]) |> String |> x -> strip(x, [' ', '\0']))]
    for i = 3:23
        htu = [htu; (varname[110 + i], Char.(hc[(i * 8 + 1):(i * 8 + 8)]) |> String |> x -> strip(x, [' ', '\0']))]
    end
    head = Dict(htu)
    # numerator and logical
    for i = HEADER.sac.enumeratevars
        if !isnan(head[i])
            head[i] = HEADER.sac.headtrans.int2other[i][head[i]]
        end
    end
    for i = HEADER.sac.logicalvars
        if !isnan(head[i])
            head[i] = Bool(head[i])
        end
    end
    return head
end

"""
readsac

    read sac binary file

    =========================

    readsac(path::AbstractString)

    =========================

    readsac(io::IO)
"""
function readsac(io::IO)
    # read head
    head = readsachead(io)
    # read data
    td = zeros(Float32, Int(head["npts"]))
    if head["npts"] >= 1 && (head["iftype"] in ["IXY", "IRLIM", "IAMPH", "IXYZ", "ITIME"]) && !isnan(head["leven"])
        read!(io, td)
        data = [td]
        if !head["leven"] || (head["iftype"] in ["IAMPH", "IRLIM"])
            read!(io, td)
            data = [data, td]
        end
    else
        data = []
    end
    return SACFrame(head, data)
end

function readsac(path::AbstractString)
    return open(x -> readsac(x), path, "r")
end

"""
writesac

    this function won't check if the head is right!!!

    =======================

    writesac(io::IO, data::SACFrame)

    =======================

    writesac(path:AbstractString="./seismogram.sac", data::SACFrame)
"""
# todo   add function to check header and data
# !  this function don't check right now!!
function writesac(io::IO, data::SACFrame=SACFrame(Dict(), []))
    for i = 1:110
        hname = HEADER.sac.headlist[i]
        if isequal(data.head[hname], NaN)
            if i <= 70
                write(io, Float32(-12345.0))
            else
                write(io, Int32(-12345))
            end
        else
            if i <= 70
                write(io, Float32(data.head[hname]))
            else
                if hname in HEADER.sac.enumeratevars
                    write(io, Int32(HEADER.sac.headtrans.other2int[hname][data.head[hname]]))
                else
                    write(io, Int32(data.head[hname]))
                end
            end
        end
    end
    hname = HEADER.sac.headlist[111]
    v = data.head[hname]
    for i = 1:8
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    hname = HEADER.sac.headlist[112]
    v = data.head[hname]
    for i = 1:16
        if i <= length(v)
            write(io, UInt8(v[i]))
        else
            write(io, UInt8(' '))
        end
    end
    for i = 113:length(HEADER.sac.headlist)
        hname = HEADER.sac.headlist[i]
        v = data.head[hname]
        for j = 1:8
            if j <= length(v)
                write(io, UInt8(v[j]))
            else
                write(io, UInt8(' '))
            end
        end
    end
    for i = data.data
        write(io, i)
    end
end

function writesac(path::AbstractString, data::SACFrame=SACFrame(Dict(), []))
    open(x -> writesac(x, data), path, "w")
end

function writesac(data::SACFrame=SACFrame(Dict(), []))
    open(x -> writesac(x, data), "./seismogram.SAC", "w")
end

"""
newsachead()
=========================
return new sac file head
"""
function newsachead()
    htu = []
    varname = HEADER.sac.headlist
    for i = 1:105
        htu = [htu; (varname[i], NaN32)]
    end
    for i = 106:110
        htu = [htu; (varname[i], false)]
    end
    for i = 1:23
        htu = [htu; (varname[110 + i], "-12345")]
    end
    head = Dict(htu)
    return head
end

function WaveFrame(s::SACFrame)
    return WaveFrame("sac", s.head, s.data)
end

function SACFrame(s::WaveFrame)
    headkeys = keys(s.head)
    h = newsachead()
    for k in HEADER.sac.headlist
        if k in headkeys
            h[k] = s.head[k]
        end
    end
    return SACFrame(h, s.data)
end

function stdname(s::SACFrame, type::Int=1, qtag::AbstractChar='D')
    if type == 1
        return @sprintf("%04d.%03d.%02d.%02d.%02d.%04d.%s.%s.%s.%s.%s.SAC",
        s.head["nzyear"], s.head["nzjday"], s.head["nzhour"], s.head["nzmin"], s.head["nzsec"], s.head["nzmsec"],
        s.head["knetwk"], s.head["kstnm"], s.head["khole"], s.head["kcmpnm"], qtag)
    elseif type == 2
        return @sprintf("%s.%s.%s.%s.%s.%04d.%03d.%02d%02d%02d.SAC",
        s.head["knetwk"], s.head["kstnm"], s.head["khole"], s.head["kcmpnm"], qtag,
        s.head["nzyear"], s.head["nzjday"], s.head["nzhour"], s.head["nzmin"], s.head["nzsec"])
    end
end