module SEGY
import Base.write

include("macros.jl")

const EBCDIC_TABLE = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                      27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                      51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 74, 75, 76, 77, 78, 79, 80, 90, 91, 92,
                      93, 94, 95, 96, 97, 106, 107, 108, 109, 110, 111, 121, 122, 123, 124, 125, 126, 127, 129, 130,
                      131, 132, 133, 134, 135, 136, 137, 145, 146, 147, 148, 149, 150, 151, 152, 153, 161, 162, 163,
                      164, 165, 166, 167, 168, 169, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 208, 209, 210,
                      211, 212, 213, 214, 215, 216, 217, 224, 226, 227, 228, 229, 230, 231, 232, 233, 240, 241, 242,
                      243, 244, 245, 246, 247, 248, 249, 255, 0, 1, 2, 3, 55, 45, 46, 47, 47, 22, 5, 37, 11, 12, 13, 16,
                      17, 18, 19, 60, 61, 50, 38, 24, 63, 39, 28, 29, 29, 30, 31, 7, 64, 90, 127, 123, 91, 108, 80, 125,
                      77, 93, 92, 78, 107, 96, 75, 97, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 122, 94, 76,
                      126, 110, 111, 124, 0, 224, 0, 0, 109, 121, 192, 79, 208, 161)

const ASCII_TABLE = (0, 1, 2, 3, 156, 9, 134, 127, 151, 141, 142, 11, 12, 13, 14, 15, 16, 17, 18, 19, 157, 133, 8, 135,
                     24, 25, 146, 143, 28, 29, 30, 31, 128, 129, 130, 131, 132, 10, 23, 27, 136, 137, 138, 139, 140, 5,
                     6, 7, 144, 145, 22, 147, 148, 149, 150, 4, 152, 153, 154, 155, 20, 21, 158, 26, 32, 162, 46, 60,
                     40, 43, 124, 38, 33, 36, 42, 41, 59, 172, 45, 47, 166, 44, 37, 95, 62, 63, 96, 58, 35, 64, 39, 61,
                     34, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 126,
                     115, 116, 117, 118, 119, 120, 121, 122, 123, 65, 66, 67, 68, 69, 70, 71, 72, 73, 125, 74, 75, 76,
                     77, 78, 79, 80, 81, 82, 92, 83, 84, 85, 86, 87, 88, 89, 90, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
                     159, 0, 1, 2, 3, 4, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27,
                     28, 29, 29, 30, 31, 127, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                     50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 91, 92, 93, 94, 95, 96, 123, 124, 125,
                     126)
const FILE_HEADER_VAR_LIST = ("job", "line", "reel", "dataTracePerEnsemble", "auxiliaryTracePerEnsemble", "dt",
                              "dtOrig", "ns", "nsOrig", "dataSampleFormat", "ensembleFold", "traceSorting",
                              "verticalSumCode", "sweepFrequencyStart", "sweepFrequencyEnd", "sweepLength", "sweepType",
                              "sweepChannel", "sweepTaperLengthStart", "sweepTaperLengthEnd", "taperType",
                              "correlatedDataTraces", "binaryGain", "amplitudeRecoveryMethod", "measurementSystem",
                              "impulseSignalPolarity", "vibratoryPolarityCode", "segyFormatRevisionNumber",
                              "fixedLengthTraceFlag", "numberOfExtTextualHeaders")
const FILE_HEADER_VAR_TYPE_LIST = (Int32, Int32, Int32, Int16, UInt16, UInt16, UInt16, UInt16, UInt16, UInt16, Int16,
                                   UInt16, UInt16, UInt16, UInt16, Int16, UInt16, UInt16, UInt16, UInt16, Int16, UInt16,
                                   UInt16, UInt16, UInt16, UInt16, UInt16, UInt16, Int16, UInt16)
const TRACE_HEADER_VAR_LIST = ("traceSequenceLine", "traceSequenceFile", "fileRecord", "traceNumber",
                               "energySourcePoing", "cdp", "cdpTrace", "traceIdenitifactionCode", "nSummedTraces",
                               "nStackedTraces", "dataUse", "offset", "receiverGroupElevation",
                               "sourceSurfaceElevation", "sourceDepth", "receiverDatumElevation",
                               "sourceDatumElevation", "sourceWaterDepth", "groupWaterDepth", "elevationScalar",
                               "sourceGroupScalar", "sourceX", "sourceY", "groupX", "groupY", "coordinateUnits",
                               "weatheringVelocity", "subWeatheringVelocity", "sourceUpholeTime", "groupUpholeTime",
                               "sourceStaticCorrection", "groupStaticCorrection", "totalStaticApplied", "lagTimeA",
                               "lagTimeB", "delayRecordingTime", "muteTimeStart", "muteTimeEnd", "ns", "dt", "gainType",
                               "instrumentGainConstant", "instrumentInitialGain", "correlated", "sweepFrequencyStart",
                               "sweepFrequencyEnd", "sweepLength", "sweepType", "sweepTraceTaperLengthStart",
                               "sweepTraceTaperLengthEnd", "taperType", "aliasFilterFrequency", "aliasFilterSlope",
                               "notchFilterFrequency", "notchFilterSlope", "lowCutFrequency", "highCutFrequency",
                               "lowCutSlope", "highCutSlope", "yearDataRecorded", "dayOfYear", "hourOfDay",
                               "minuteOfHour", "secondOfMinute", "timeBaseCode", "traceWeightningFactor",
                               "geophoneGroupNumberOfRoll1", "geophoneGroupNumberFirstTraceOrigField",
                               "geophoneGroupNumberLastTraceOrigField", "gapSize", "overTravel", "cdpX", "cdpY",
                               "inline3D", "crossline3D", "shotPoint", "shotPointScalar", "traceValueMeasurementUnit",
                               "transductionConstantMantissa", "transductionConstantPower", "transductionUnit",
                               "traceIdentifier", "scalarTraceHeader", "sourceType", "sourceEnergyDirectionMantissa",
                               "sourceEnergyDirectionExponent", "sourceMeasurementMantissa",
                               "sourceMeasurementExponent", "sourceMeasurementUnit", "unassignedInt1", "unassignedInt2")
const TRACE_HEADER_VAR_TYPE_LIST = (Int32, Int32, Int32, Int32, Int32, Int32, Int32, Int16, Int16, Int16, Int16, Int32,
                                    Int32, Int32, Int32, Int32, Int32, Int32, Int32, Int16, Int16, Int32, Int32, Int32,
                                    Int32, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16,
                                    Int16, Int16, UInt16, UInt16, Int16, Int16, Int16, Int16, Int16, Int16, Int16,
                                    Int16,
                                    Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16,
                                    Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int32,
                                    Int32, Int32, Int32, Int32, UInt16, UInt16, Int32, Int16, Int16, Int16, Int16,
                                    Int16,
                                    Int32, Int16, Int32, Int16, Int16, Int32, Int32)

function ebcdic2ascii(c::Char)
    I = findfirst(EBCDIC_TABLE .== Int(c))
    if isnothing(I)
        I = findfirst(EBCDIC_TABLE .== 0)
    end
    return Char(ASCII_TABLE[I])
end

function filldict!(d::Dict, k, t)
    for i in eachindex(k)
        d[k[i]] = t[i]
    end
end

function readfilehead(io::IO)
    # taperLabel = Base.read(io, 128) |> String
    taperLabel                                  = ""
    textualFileHead                             = Base.read(io, 3200) |> String
    binaryFileHead                              = Dict{String,Any}()
    binaryFileHead["job"]                       = Base.read(io, Int32) |> ntoh |> Int
    binaryFileHead["line"]                      = Base.read(io, Int32) |> ntoh |> Int
    binaryFileHead["reel"]                      = Base.read(io, Int32) |> ntoh |> Int
    binaryFileHead["dataTracePerEnsemble"]      = Base.read(io, Int16) |> ntoh |> Int
    binaryFileHead["auxiliaryTracePerEnsemble"] = Base.read(io, UInt16) |> ntoh |> Int
    binaryFileHead["dt"]                        = Base.read(io, UInt16) |> ntoh |> Int
    binaryFileHead["dtOrig"]                    = Base.read(io, UInt16) |> ntoh |> Int
    binaryFileHead["ns"]                        = Base.read(io, UInt16) |> ntoh |> Int
    binaryFileHead["nsOrig"]                    = Base.read(io, UInt16) |> ntoh |> Int
    t                                           = zeros(Int16, 138)
    read!(io, t)
    t = Int.(ntoh.(t[1:18]))
    k = ("dataSampleFormat", "ensembleFold", "traceSorting", "verticalSumCode", "sweepFrequencyStart",
         "sweepFrequencyEnd", "sweepLength", "sweepType", "sweepChannel", "sweepTaperLengthStart",
         "sweepTaperLengthEnd", "taperType", "correlatedDataTraces", "binaryGain", "amplitudeRecoveryMethod",
         "measurementSystem", "impulseSignalPolarity", "vibratoryPolarityCode")
    filldict!(binaryFileHead, k, t)
    binaryFileHead["unassigned1"]               = t[19:end]
    binaryFileHead["segyFormatRevisionNumber"]  = Base.read(io, UInt16) |> ntoh |> Int
    binaryFileHead["fixedLengthTraceFlag"]      = Base.read(io, Int16) |> ntoh |> Int
    binaryFileHead["numberOfExtTextualHeaders"] = Base.read(io, UInt16) |> ntoh |> Int
    t                                           = zeros(Int16, 47)
    read!(io, t)
    t = Int.(ntoh.(t))
    binaryFileHead["unassigned2"] = t
    extendedTextualFileHead = String[]
    if binaryFileHead["numberOfExtTextualHeaders"] > 0
        for i = 1:binaryFileHead["numberOfExtTextualHeaders"]
            s = Base.read(io, 3200) |> String
            push!(extendedTextualFileHead, s)
        end
    end
    return (taperlabel = taperLabel, txtfh = textualFileHead, bfh = binaryFileHead, etxtfh = extendedTextualFileHead)
end

function readtracehead(io::IO)
    th = Dict{String,Real}()

    t = zeros(Int32, 7)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("traceSequenceLine", "traceSequenceFile", "fileRecord", "traceNumber", "energySourcePoing", "cdp", "cdpTrace")
    filldict!(th, k, t)

    t = zeros(Int16, 4)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("traceIdenitifactionCode", "nSummedTraces", "nStackedTraces", "dataUse")
    filldict!(th, k, t)

    t = zeros(Int32, 8)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("offset", "receiverGroupElevation", "sourceSurfaceElevation", "sourceDepth", "receiverDatumElevation",
         "sourceDatumElevation", "sourceWaterDepth", "groupWaterDepth")
    filldict!(th, k, t)

    th["elevationScalar"] = Base.read(io, Int16) |> ntoh |> Int
    th["sourceGroupScalar"] = Base.read(io, Int16) |> ntoh |> Int

    t = zeros(Int32, 4)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("sourceX", "sourceY", "groupX", "groupY")
    filldict!(th, k, t)

    t = zeros(Int16, 13)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("coordinateUnits", "weatheringVelocity", "subWeatheringVelocity", "sourceUpholeTime", "groupUpholeTime",
         "sourceStaticCorrection", "groupStaticCorrection", "totalStaticApplied", "lagTimeA", "lagTimeB",
         "delayRecordingTime", "muteTimeStart", "muteTimeEnd")
    filldict!(th, k, t)

    th["ns"] = Base.read(io, UInt16) |> ntoh |> Int
    th["dt"] = Base.read(io, UInt16) |> ntoh |> Int

    t = zeros(Int16, 31)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("gainType", "instrumentGainConstant", "instrumentInitialGain", "correlated", "sweepFrequencyStart",
         "sweepFrequencyEnd", "sweepLength", "sweepType", "sweepTraceTaperLengthStart", "sweepTraceTaperLengthEnd",
         "taperType", "aliasFilterFrequency", "aliasFilterSlope", "notchFilterFrequency", "notchFilterSlope",
         "lowCutFrequency", "highCutFrequency", "lowCutSlope", "highCutSlope", "yearDataRecorded", "dayOfYear",
         "hourOfDay", "minuteOfHour", "secondOfMinute", "timeBaseCode", "traceWeightningFactor",
         "geophoneGroupNumberOfRoll1", "geophoneGroupNumberFirstTraceOrigField",
         "geophoneGroupNumberLastTraceOrigField", "gapSize", "overTravel")
    filldict!(th, k, t)

    t = zeros(Int32, 5)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("cdpX", "cdpY", "inline3D", "crossline3D", "shotPoint")
    filldict!(th, k, t)

    th["shotPointScalar"] = Base.read(io, Int16) |> ntoh |> Int
    th["traceValueMeasurementUnit"] = Base.read(io, Int16) |> ntoh |> Int
    th["transductionConstantMantissa"] = Base.read(io, Int32) |> ntoh |> Int
    t = zeros(Int16, 5)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("transductionConstantPower", "transductionUnit", "traceIdentifier", "scalarTraceHeader", "sourceType")
    filldict!(th, k, t)

    th["sourceEnergyDirectionMantissa"] = Base.read(io, Int32) |> ntoh |> Int
    th["sourceEnergyDirectionExponent"] = Base.read(io, Int16) |> ntoh |> Int
    th["sourceMeasurementMantissa"]     = Base.read(io, Int32) |> ntoh |> Int
    th["sourceMeasurementExponent"]     = Base.read(io, Int16) |> ntoh |> Int
    th["sourceMeasurementUnit"]         = Base.read(io, Int16) |> ntoh |> Int
    th["unassignedInt1"]                = Base.read(io, Int32) |> ntoh |> Int
    th["unassignedInt2"]                = Base.read(io, Int32) |> ntoh |> Int

    return th
end

function readtrace(io::IO, fhdr::Dict)
    hdr = readtracehead(io)
    if fhdr["dataSampleFormat"] == 5
        t = zeros(Float32, fhdr["ns"])
    elseif fhdr["dataSampleFormat"] == 6
        t = zeros(Float64, fhdr["ns"])
    else
        @error "data sample format not supported now."
    end
    read!(io, t)
    t = Float64.(ntoh.(t))
    return (hdr = hdr, data = t)
end

function read(io::IO)
    fh = readfilehead(io)
    traces = []
    while !eof(io)
        t = readtrace(io, fh.bfh)
        push!(traces, t)
    end
    return (filehead = fh, traces = traces)
end

function read(p::AbstractString)
    return open(read, p, "r")
end

function writefilehead(io::IO; taperlabel::String = "", textualFileHead::String = "",
                       binaryFileHead::Dict = Dict(), extendedTextualFileHead::Vector{String} = String[])
    sbuf = zeros(UInt8, 3200)
    sbuf[1:length(textualFileHead)] .= UInt8.(collect(textualFileHead))
    Base.write(io, sbuf)
    for i = 1:27
        T = FILE_HEADER_VAR_TYPE_LIST[i]
        Base.write(io, T(binaryFileHead[FILE_HEADER_VAR_LIST[i]]))
    end
    unas = zeros(Int8, 240)
    Base.write(io, unas)
    for i = 28:30
        T = FILE_HEADER_VAR_TYPE_LIST[i]
        Base.write(io, T(binaryFileHead[FILE_HEADER_VAR_LIST[i]]))
    end
    unas = zeros(Int8, 47)
    Base.write(io, unas)
    for i in eachindex(extendedTextualFileHead)
        sbuf .= '\0'
        sbuf[1:length(extendedTextualFileHead[i])] .= UInt8.(collect(extendedTextualFileHead[i]))
        Base.write(io, sbuf)
    end
    return nothing
end

function writetrace(io::IO, theader::Dict, data::Vector{<:Real})
    for i = 1:91
        T = TRACE_HEADER_VAR_TYPE_LIST[i]
        Base.write(io, T(theader[TRACE_HEADER_VAR_LIST[i]]))
    end
    Base.write(io, data)
    return nothing
end

function write(io::IO, binaryFileHead::Dict, theader::Vector{Dict}, data::Vector{Vector{<:Real}},
               textualFileHead::String = "", extendedTextualFileHead::Vector{String} = String[])
    writefilehead(io; textualFileHead = textualFileHead, binaryFileHead = binaryFileHead,
                  extendedTextualFileHead = extendedTextualFileHead)
    for i in eachindex(theader)
        writetrace(io, theader[i], data[i])
    end
end

function write(io::IO, binaryFileHead::Dict, theader::Vector{Dict}, data::Matrix{<:Real},
               textualFileHead::String = "", extendedTextualFileHead::Vector{String} = String[])
    writefilehead(io; textualFileHead = textualFileHead, binaryFileHead = binaryFileHead,
                  extendedTextualFileHead = extendedTextualFileHead)
    for i in eachindex(theader)
        writetrace(io, theader[i], data[:, i])
    end
end

end
