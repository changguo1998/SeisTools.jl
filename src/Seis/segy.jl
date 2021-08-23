const EBCDIC_TABLE = (0, 1, 2, 3, 4, 5, 6, 7,
8, 9, 10, 11, 12, 13, 14, 15,
16, 17, 18, 19, 20, 21, 22, 23,
24, 25, 26, 27, 28, 29, 30, 31,
32, 33, 34, 35, 36, 37, 38, 39,
40, 41, 42, 43, 44, 45, 46, 47,
48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63,
64, 74, 75, 76, 77, 78, 79, 80,
90, 91, 92, 93, 94, 95, 96, 97,
106, 107, 108, 109, 110, 111, 121, 122,
123, 124, 125, 126, 127, 129, 130, 131,
132, 133, 134, 135, 136, 137, 145, 146,
147, 148, 149, 150, 151, 152, 153, 161,
162, 163, 164, 165, 166, 167, 168, 169,
192, 193, 194, 195, 196, 197, 198, 199,
200, 201, 208, 209, 210, 211, 212, 213,
214, 215, 216, 217, 224, 226, 227, 228,
229, 230, 231, 232, 233, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 255,
0, 1, 2, 3, 55, 45, 46, 47,
47, 22, 5, 37, 11, 12, 13, 16,
17, 18, 19, 60, 61, 50, 38, 24,
63, 39, 28, 29, 29, 30, 31, 7,
64, 90, 127, 123, 91, 108, 80, 125,
77, 93, 92, 78, 107, 96, 75, 97,
240, 241, 242, 243, 244, 245, 246, 247,
248, 249, 122, 94, 76, 126, 110, 111,
124, 0, 224, 0, 0, 109, 121, 192, 79, 208, 161)

const ASCII_TABLE = (0, 1, 2, 3, 156, 9, 134, 127,
151, 141, 142, 11, 12, 13, 14, 15,
16, 17, 18, 19, 157, 133, 8, 135,
24, 25, 146, 143, 28, 29, 30, 31,
128, 129, 130, 131, 132, 10, 23, 27,
136, 137, 138, 139, 140, 5, 6, 7,
144, 145, 22, 147, 148, 149, 150, 4,
152, 153, 154, 155, 20, 21, 158, 26,
32, 162, 46, 60, 40, 43, 124, 38,
33, 36, 42, 41, 59, 172, 45, 47,
166, 44, 37, 95, 62, 63, 96, 58,
35, 64, 39, 61, 34, 97, 98, 99,
100, 101, 102, 103, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 126,
115, 116, 117, 118, 119, 120, 121, 122,
123, 65, 66, 67, 68, 69, 70, 71,
72, 73, 125, 74, 75, 76, 77, 78,
79, 80, 81, 82, 92, 83, 84, 85,
86, 87, 88, 89, 90, 48, 49, 50,
51, 52, 53, 54, 55, 56, 57, 159,
0, 1, 2, 3, 4, 5, 6, 7,
7, 8, 9, 10, 11, 12, 13, 16,
17, 18, 19, 20, 21, 22, 23, 24,
26, 27, 28, 29, 29, 30, 31, 127,
32, 33, 34, 35, 36, 37, 38, 39,
40, 41, 42, 43, 44, 45, 46, 47,
48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63,
64, 91, 92, 93, 94, 95, 96, 123,
124, 125, 126)

function ebcdic2ascii(c::Char)
    I = findfirst(EBCDIC_TABLE .== Int(c))
    if isnothing(I)
        I = findfirst(EBCDIC_TABLE .== 0)
    end
    return Char(ASCII_TABLE[I])
end

function filldict!(d::Dict, k, t)
    for i = 1:length(k)
        d[k[i]] = t[i]
    end
end

function readsegyfilehead(io::IO)
    # taperLabel = read(io, 128) |> String
    taperLabel      = []
    textualFileHead = read(io, 3200) |> String
    binaryFileHead  = Dict()
    binaryFileHead["job"]                       = read(io, Int32) |> ntoh |> Int
    binaryFileHead["line"]                      = read(io, Int32) |> ntoh |> Int
    binaryFileHead["reel"]                      = read(io, Int32) |> ntoh |> Int
    binaryFileHead["dataTracePerEnsemble"]      = read(io, Int16) |> ntoh |> Int
    binaryFileHead["auxiliaryTracePerEnsemble"] = read(io, UInt16) |> ntoh |> Int
    binaryFileHead["dt"]                        = read(io, UInt16) |> ntoh |> Int
    binaryFileHead["dtOrig"]                    = read(io, UInt16) |> ntoh |> Int
    binaryFileHead["ns"]                        = read(io, UInt16) |> ntoh |> Int
    binaryFileHead["nsOrig"]                    = read(io, UInt16) |> ntoh |> Int
    t = zeros(Int16, 138)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("dataSampleFormat", "ensembleFold", "traceSorting", "verticalSumCode", "sweepFrequencyStart", "sweepFrequencyEnd",
    "sweepLength", "sweepType", "sweepChannel", "sweepTaperLengthStart", "sweepTaperLengthEnd", "taperType",
    "correlatedDataTraces", "binaryGain", "amplitudeRecoveryMethod", "measurementSystem", "impulseSignalPolarity",
    "vibratoryPolarityCode")
    filldict!(binaryFileHead, k, t)
    binaryFileHead["unassigned1"]               = t[19:end]
    binaryFileHead["segyFormatRevisionNumber"]  = read(io, UInt16) |> ntoh |> Int
    binaryFileHead["fixedLengthTraceFlag"]      = read(io, Int16) |> ntoh |> Int
    binaryFileHead["numberOfExtTextualHeaders"] = read(io, UInt16) |> ntoh |> Int
    t = zeros(Int16, 47)
    read!(io, t)
    t = Int.(ntoh.(t))
    binaryFileHead["unassigned2"] = t
    extendedTextualFileHead = []
    if binaryFileHead["numberOfExtTextualHeaders"] > 0
        for i = 1:binaryFileHead["numberOfExtTextualHeaders"]
            s = read(io, 3200) |> String
            push!(extendedTextualFileHead, s)
        end
    end
    return (taperlabel = taperLabel, txtfh = textualFileHead, bfh = binaryFileHead, etxtfh = extendedTextualFileHead)
end

function readsegytracehead(io::IO)
    th = Dict()

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

    th["elevationScalar"] = read(io, Int16) |> ntoh |> Int
    th["sourceGroupScalar"] = read(io, Int16) |> ntoh |> Int

    t = zeros(Int32, 4)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("sourceX", "sourceY", "groupX", "groupY")
    filldict!(th, k, t)

    t = zeros(Int16, 13)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("coordinateUnits", "weatheringVelocity", "subWeatheringVelocity", "sourceUpholeTime", "groupUpholeTime",
    "sourceStaticCorrection", "groupStaticCorrection", "totalStaticApplied", "lagTimeA", "lagTimeB", "delayRecordingTime",
    "muteTimeStart", "muteTimeEnd")
    filldict!(th, k, t)

    th["ns"] = read(io, UInt16) |> ntoh |> Int
    th["dt"] = read(io, UInt16) |> ntoh |> Int

    t = zeros(Int16, 31)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("gainType", "instrumentGainConstant", "instrumentInitialGain", "correlated", "sweepFrequencyStart", "sweepFrequencyEnd",
    "sweepLength", "sweepType", "sweepTraceTaperLengthStart", "sweepTraceTaperLengthEnd", "taperType", "aliasFilterFrequency",
    "aliasFilterSlope","notchFilterFrequency", "notchFilterSlope", "lowCutFrequency", "highCutFrequency", "lowCutSlope", "highCutSlope",
    "yearDataRecorded", "dayOfYear", "hourOfDay", "minuteOfHour", "secondOfMinute", "timeBaseCode", "traceWeightningFactor",
    "geophoneGroupNumberOfRoll1", "geophoneGroupNumberFirstTraceOrigField", "geophoneGroupNumberLastTraceOrigField", "gapSize", "overTravel")
    filldict!(th, k, t)

    t = zeros(Int32, 5)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("cdpX", "cdpY", "inline3D", "crossline3D", "shotPoint")
    filldict!(th, k, t)

    th["shotPointScalar"] = read(io, Int16) |> ntoh |> Int
    th["traceValueMeasurementUnit"] = read(io, Int16) |> ntoh |> Int
    th["transductionConstantMantissa"] = read(io, Int32) |> ntoh |> Int
    t = zeros(Int16, 5)
    read!(io, t)
    t = Int.(ntoh.(t))
    k = ("transductionConstantPower", "transductionUnit", "traceIdentifier", "scalarTraceHeader", "sourceType")
    filldict!(th, k, t)

    th["sourceEnergyDirectionMantissa"] = read(io, Int32) |> ntoh |> Int
    th["sourceEnergyDirectionExponent"] = read(io, Int16) |> ntoh |> Int
    th["sourceMeasurementMantissa"]     = read(io, Int32) |> ntoh |> Int
    th["sourceMeasurementExponent"]     = read(io, Int16) |> ntoh |> Int
    th["sourceMeasurementUnit"]         = read(io, Int16) |> ntoh |> Int
    th["unassignedInt1"]                = read(io, Int32) |> ntoh |> Int
    th["unassignedInt2"]                = read(io, Int32) |> ntoh |> Int

    return th
end

function readsegytrace(io::IO, fhdr::Dict)
    hdr = readsegytracehead(io)
    if fhdr["dataSampleFormat"] == 5
        t = zeros(Float32, fhdr["ns"])
    elseif fhdr["dataSampleFormat"] == 6
        t = zeros(Float64, fhdr["ns"])
    else
        @error "data sample format not supported now."
    end
    read!(io, t)
    t = Float64.(ntoh.(t))
    return SEGYFrame(hdr, t)
end

function readsegy(io::IO)
    fh = readsegyfilehead(io)
    traces = []
    while !eof(io)
        t = readsegytrace(io, fh.bfh)
        push!(traces, t)
    end
    return (filehead=fh, traces = traces)
end

function readsegy(p::AbstractString)
    return open(readsegy, p, "r")
end
