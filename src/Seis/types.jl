module HEADER
segyfilehead = ("job", "line", "reel", "dataTracePerEnsemble", "AuxiliaryTracePerEnsemble", "dt", "dtOrig", "ns",
                "nsOrig", "dataSampleFormat", "ensembleFold", "traceSorting", "verticalSumCode", "sweepFrequencyStart",
                "sweepFrequencyEnd", "sweepLength", "sweepType", "sweepChannel", "sweepTaperLengthStart",
                "sweepTaperLengthEnd", "taperType", "correlatedDataTraces", "binaryGain", "amplitudeRecoveryMethod",
                "measurementSystem", "impulseSignalPolarity", "vibratoryPolarityCode", "unassigned1",
                "segyFormatRevisionNumber", "fixedLengthTraceFlag", "numberOfExtTextualHeaders", "unassigned2")
segyfileheadtype = (Int32, Int32, Int32, Int16, UInt16, UInt16, UInt16, UInt16, UInt16, Int16, Int16, Int16, Int16,
                    Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16, Int16)

end
