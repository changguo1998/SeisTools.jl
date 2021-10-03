abstract type Frame <: Any end

struct WaveFrame <: Frame
    network::AbstractString
    station::AbstractString
    device::AbstractString
    component::AbstractString
    id::AbstractString
    begintime::DateTime
    delta::Peroid
    npts::Unsigned
    meta::Dict{String,Any}
    data::AbstractVector
    function WaveFrame(n::AbstractString = "", s::AbstractString = "", dv::AbstractString = "",
                       cmp::AbstractString = "", bt::DateTime = DateTime(2021), dt::Peroid = Millisecond(10),
                       data::AbstractVector = [], npts::Unsigned = 0, meta = Dict{String,Any}();
                       id::AbstractString = "", fill_meta::Bool = false)
        if isempty(id)
            id = join([n, s, dv, cmp], '.')
        end
        if npts == 0
            npts = length(data)
        end
        if isempty(meta) || fill_meta
            meta["network"] = n
            meta["staiton"] = s
            meta["device"] = dv
            meta["component"] = cmp
            meta["begintime"] = bt
            meta["delta"] = dt
            meta["npts"] = npts
            meta["id"] = id
        end
        return new(n, s, dv, cmp, id, bt, dt, npts, meta, data)
    end
end
module HEADER
segyfilehead = ("job", "line", "reel", "dataTracePerEnsemble", "AuxiliaryTracePerEnsemble", "dt",
                "dtOrig", "ns", "nsOrig", "dataSampleFormat", "ensembleFold", "traceSorting",
                "verticalSumCode", "sweepFrequencyStart", "sweepFrequencyEnd", "sweepLength",
                "sweepType", "sweepChannel", "sweepTaperLengthStart", "sweepTaperLengthEnd",
                "taperType", "correlatedDataTraces", "binaryGain", "amplitudeRecoveryMethod",
                "measurementSystem", "impulseSignalPolarity", "vibratoryPolarityCode",
                "unassigned1", "segyFormatRevisionNumber", "fixedLengthTraceFlag",
                "numberOfExtTextualHeaders", "unassigned2")
segyfileheadtype = (Int32, Int32, Int32, Int16, UInt16, UInt16, UInt16, UInt16, UInt16, Int16,
                    Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16,
                    Int16, Int16, Int16, Int16, Int16, Int16)

end
