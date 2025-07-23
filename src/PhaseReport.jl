module PhaseReport

using Dates, Printf
import Base: isequal, isless, print

abstract type Phase <: Any end

# = =================================================================
# China Earthquake Network Center
# = =================================================================

# GX 2013/01/16 23:19:55.0  25.555  105.708   7  4.0     1  19 eq 52 贵州贞丰
# GZ ZFT   BHZ I U Pg      1.0 V  23:19:57.72  -0.75   20.0 203.1
#          BHN     Sg      1.0 V  23:20:00.61  -0.29
#          BHN     SMN     1.0 D  23:20:01.12                       49362.2   0.27
#          BHE     SME     1.0 D  23:20:01.16                      132841.0   0.34 ML   4.3
# GZ AST   BHZ     Pg      1.0 V  23:20:07.75   0.27   74.0  27.5
#          BHE     Sg      1.0 V  23:20:17.75   1.59

const REF_EVENT_CENC = [(1, 2, 's'), (4, 24, 'd'), (26, 33, 'f'), (34, 42, 'f'), (43, 46, 'f'),
                   (47, 51, 'f'), (52, 56, 'i'), (57, 60, 'i'), (62, 63, 's'), (64, 66, 'i'),
                   (67, 0, 's')]
const REF_PHASE_CENC = [(1, 2, 's'), (4, 8, 's'), (10, 12, 's'), (14, 14, 'c'), (16, 16, 'c'),
                   (18, 24, 's'), (26, 28, 'f'), (30, 30, 'c'), (33, 43, 't'), (45, 51, 'f'),
                   (53, 58, 'f'), (59, 63, 'f'), (64, 73, 'f'), (74, 80, 'f'), (82, 83, 's'),
                   (84, 0, 'f')]
const RERANGE_EVENT_CENC = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
const RERANGE_PHASE_CENC = [1, 2, 3, 6, 9, 5, 15, 16, 4, 8, 7, 10, 11, 12, 13, 14]

struct CENCPhase
    network::String
    station::String
    channel::String
    type::String
    time::Time
    polarity::Char
    magtype::String
    mag::Float64
    code1::Char
    code2::Char
    value1::Float64
    value2::Float64
    value3::Float64
    value4::Float64
    value5::Float64
    value6::Float64
end

CENCEvent = NamedTuple{(:region,:par,:phase),Tuple{Tuple{String,String},Tuple{DateTime,Float64,Float64,Float64,Float64},Vector{CENCPhase}}}

function CENCPhase(network::String, station::String, channel::String, type::String, time::Time;
               polarity::Char='\0', mtype::String="", mag::Float64=NaN, code1::Char='\0', code2::Char='\0',
               v1::Float64=NaN, v2::Float64=NaN, v3::Float64=NaN, v4::Float64=NaN, v5::Float64=NaN, v6::Float64=NaN)
    return CENCPhase(network, station, channel, type, time, polarity, mtype, mag, code1, code2, v1, v2, v3, v4, v5, v6)
end

function _equal_must(x, y)
    return x === y
end

function _equal_default(x, y, d)
    return (x === y) || (x === d) || (y === d)
end

export isless

function isless(a::CENCPhase, b::CENCPhase)
    if a.network < b.network
        return true
    elseif a.network > b.network
        return false
    else
        if a.station < b.station
            return true
        elseif a.station > b.station
            return false
        else
            if a.channel < b.channel
                return true
            elseif a.channel > b.channel
                return false
            else
                if a.type < b.type
                    return true
                elseif a.type > b.type
                    return false
                else
                    if a.time < b.time
                        return true
                    else
                        return false
                    end
                end
            end
        end
    end
end

export isequal

function isequal(a::CENCPhase, b::CENCPhase)
    f = true
    for fld in (:network, :station, :channel, :type, :time)
        f &= _equal_must(getfield(a, fld), getfield(b, fld))
    end
    for fld in (:code1, :polarity, :code2)
        f &= _equal_default(getfield(a, fld), getfield(b, fld), '\0')
    end
    for fld in (:value1, :value2, :value3, :value4, :value5, :value6, :mag)
        f &= _equal_default(getfield(a, fld), getfield(b, fld), NaN)
    end
    for fld in (:magtype,)
        f &= _equal_default(getfield(a, fld), getfield(b, fld), "")
    end
    return f
end

function _select_default(x, y, d)
    if !(x === d)
        return x
    elseif !(y === d)
        return y
    else
        return d
    end
end

export mergephase

function mergephase(a::CENCPhase, b::CENCPhase)
    if !isequal(a, b)
        error("phase a and phase b not equal")
    end
    network = a.network
    station = a.station
    channel = a.channel
    type = a.type
    time = a.time
    polarity = _select_default(a.polarity, b.polarity, '\0')
    mtype = _select_default(a.magtype, b.magtype, "")
    mag = _select_default(a.mag, b.mag, NaN)
    code1 = _select_default(a.code1, b.code1, '\0')
    code2 = _select_default(a.code2, b.code2, '\0')
    v1 = _select_default(a.value1, b.value1, NaN)
    v2 = _select_default(a.value2, b.value2, NaN)
    v3 = _select_default(a.value3, b.value3, NaN)
    v4 = _select_default(a.value4, b.value4, NaN)
    v5 = _select_default(a.value5, b.value5, NaN)
    v6 = _select_default(a.value6, b.value6, NaN)
    return CENCPhase(network, station, channel, type, time, polarity, mtype, mag, code1, code2, v1, v2, v3, v4, v5, v6)
end

function _split_parse_CENC(l::AbstractString, ref::Vector{Tuple{Int,Int,Char}}, rr::Vector{Int})
    v = Vector{Any}(undef, length(ref))
    for ir in eachindex(ref)
        r = ref[ir]
        if iszero(r[2])
            b = strip(l[r[1]:end])
        else
            b = strip(l[r[1]:r[2]])
        end
        if (r[3] != 'd') && contains(b, ' ')
            error("error while parsing:\n$(l)\nsegment: $(b)")
        end
        if r[3] == 'i'
            v[ir] = isempty(b) ? NaN : parse(Int, b)
        elseif r[3] == 'f'
            v[ir] = isempty(b) ? NaN : parse(Float64, b)
        elseif r[3] == 'd'
            v[ir] = DateTime(b, dateformat"Y/m/d H:M:S.s")
        elseif r[3] == 't'
            v[ir] = Time(b, dateformat"H:M:S.s")
        elseif r[3] == 'c'
            v[ir] = isempty(b) ? '\0' : b[1]
        elseif r[3] == 's'
            v[ir] = String(b)
        end
    end
    return v[rr]
end

function _phasereport_event_CENC(l::Vector{<:AbstractString})
    (evtregionCode, evtdate, evtlat, evtlon, evtdep,
    evtmag, _, _, _, _, evtregionName) = _split_parse_CENC(l[1], REF_EVENT_CENC, RERANGE_EVENT_CENC)

    curNetwork = ""
    curStation = ""
    phases = CENCPhase[]
    for il in eachindex(l)
        if il == 1
            continue
        end
        v = _split_parse_CENC(l[il], REF_PHASE_CENC, RERANGE_PHASE_CENC)
        if !isempty(v[1])
            curNetwork = v[1]
        end
        if !isempty(v[2])
            curStation = v[2]
        end
        v[1] = curNetwork
        v[2] = curStation
        push!(phases, CENCPhase(v...))
    end
    return (region=(evtregionCode, evtregionName),
            par=(evtdate, evtlat, evtlon, evtdep, evtmag),
            phase=phases)
end

export scan_CENC

function scan_CENC(io::IO)
    lines = readlines(io; keep=true)
    filter!(!isempty, lines)
    starts = findall(contains("eq"), lines)
    ends = [starts[2:end] .- 1; length(lines)]
    evts = Vector{CENCEvent}(undef, length(starts))
    for il in eachindex(starts)
        evts[il] = _phasereport_event_CENC(lines[starts[il]:ends[il]])
    end
    return evts
end

function scan_CENC(fname::AbstractString)
    return open(scan_CENC, fname, "r")
end

function _fill_buffer_CENC!(b::Vector{Char}, i::Int, j::Int, info::String)
    k = min(j, i+length(info)-1)
    q = collect(info)
    for p = i:k
        b[p] = Char(q[p-i+1])
    end
    return nothing
end

function print(io::IO, pr::Vector{CENCEvent})
    buffer = Char.(collect(" "^90))
    for e in pr
        buffer .= ' '
        _fill_buffer_CENC!(buffer, 1, 2, e.region[1])
        et = e.par[1]
        _fill_buffer_CENC!(buffer, 4, 24, @sprintf("%04d/%02d/%02d %02d:%02d:%02d.%01d", year(et), month(et), day(et),
            hour(et), minute(et), second(et), round(Int, millisecond(et)/100)))
        if !isnan(e.par[2])
            _fill_buffer_CENC!(buffer, 26, 32, @sprintf("%7.3f", e.par[2]))
        end
        if !isnan(e.par[3])
            _fill_buffer_CENC!(buffer, 34, 41, @sprintf("%8.3f", e.par[3]))
        end
        if !isnan(e.par[4])
            _fill_buffer_CENC!(buffer, 43, 45, @sprintf("%3d", round(Int, e.par[4])))
        end
        if !isnan(e.par[5])
            _fill_buffer_CENC!(buffer, 47, 50, @sprintf("%4.1f", e.par[5]))
        end
        _fill_buffer_CENC!(buffer, 62, 63, "eq")
        _fill_buffer_CENC!(buffer, 68, 89, e.region[2])
        buffer[90] = '\n'
        Base.print(io, String(buffer))
        cp = [""]
        for p in sort(e.phase, lt=isless)
            buffer .= ' '
            tag = p.network*"."*p.station
            if tag != cp[1]
                cp[1] = tag
                _fill_buffer_CENC!(buffer, 1, 2, p.network)
                _fill_buffer_CENC!(buffer, 4, 8, p.station)
            end
            _fill_buffer_CENC!(buffer, 10, 12, p.channel)
            if p.code1 != '\0'
                buffer[14] = p.code1
            end
            if p.polarity != '\0'
                buffer[16] = p.polarity
            end
            _fill_buffer_CENC!(buffer, 18, 24, p.type)
            _fill_buffer_CENC!(buffer, 26, 28, @sprintf("%3.1f", p.value1))
            if p.code2 != '\0'
                buffer[30] = p.code2
            end
            _fill_buffer_CENC!(buffer, 33, 43, @sprintf("%02d:%02d:%02d.%02d", hour(p.time), minute(p.time), second(p.time),
                round(Int, millisecond(p.time)/10)))
            if !isnan(p.value2)
                _fill_buffer_CENC!(buffer, 45, 50, @sprintf("%6.2f", p.value2))
            end
            if !isnan(p.value3)
                _fill_buffer_CENC!(buffer, 52, 57, @sprintf("%6.1f", p.value3))
            end
            if !isnan(p.value4)
                _fill_buffer_CENC!(buffer, 59, 63, @sprintf("%5.1f", p.value4))
            end
            if !isnan(p.value5)
                _fill_buffer_CENC!(buffer, 65, 73, @sprintf("%9.1f", p.value5))
            end
            if !isnan(p.value6)
                _fill_buffer_CENC!(buffer, 75, 80, @sprintf("%6.2f", p.value6))
            end
            if !isempty(p.magtype)
                _fill_buffer_CENC!(buffer, 82, 83, p.magtype)
            end
            if !isnan(p.mag)
                _fill_buffer_CENC!(buffer, 85, 89, @sprintf("%5.1f", p.mag))
            end
            buffer[90] = '\n'
            Base.print(io, String(buffer))
        end
    end
    return nothing
end

end # module PhaseReport
