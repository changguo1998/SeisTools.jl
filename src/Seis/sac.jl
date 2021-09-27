module SAC
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
enumeratevars = ("iftype", "idep", "iztype", "ievtyp", "iqual", "isynth", "imagtyp", "imagsrc")
logicalvars = ("leven", "lpspol", "lovrok", "lcalda", "unused18")

iftype = Dict([(1, "ITIME"), (2, "IRLIM"), (3, "IAMPH"), (4, "IXY"), (51, "IXYZ")])
idep = Dict([(5, "IUMKN"), (6, "IDISP"), (7, "IVEL"), (8, "IACC"), (50, "IVOLTS")])
iztype = Dict([(5, "IUNKN"), (9, "IB"), (10, "IDAY"), (11, "IO"), (12, "IA"), (13, "IT0"),
               (14, "IT1"), (15, "IT2"), (16, "IT3"), (17, "IT4"), (18, "IT5"), (19, "IT6"),
               (20, "IT7"), (21, "IT8"), (22, "IT9")])
ievtyp = Dict([(5, "IUNKN"), (37, "INUCL"), (38, "IPREN"), (39, "IPOSTN"), (40, "IQUAKE"),
               (41, "IPREQ"), (42, "IPOSTQ"), (43, "ICHEM"), (44, "IOTHER"), (70, "IQB"),
               (71, "IQB1"), (72, "IQB2"), (73, "IQBX"), (74, "IQMT"), (75, "IEQ"), (76, "IEQ1"),
               (77, "IEQ2"), (78, "IME"), (79, "IEX"), (80, "INU"), (81, "INC"), (82, "IO_"),
               (83, "IL"), (84, "IR"), (85, "IT"), (86, "IU")])
iqual = Dict([(44, "IOTHER"), (45, "IGOOD"), (46, "IGLCH"), (47, "IDROP"), (48, "ILOSN")])
isynth = Dict([(49, "IRLDA")])
imagtyp = Dict([(52, "IMB"), (53, "IMS"), (54, "IML"), (55, "IMW"), (56, "IMD"), (57, "IMX")])
imagsrc = Dict([(58, "INEIC"), (59, "IPDE"), (60, "IISC"), (61, "IREB"), (62, "IUSGS"),
                (63, "IBRK"), (64, "ICALTECH"), (65, "ILLNL"), (66, "IEVLOC"), (67, "IJSOP"),
                (68, "IUSER"), (69, "IUNKNOWN")])

function readhead(io::IO) end
function read(io::IO) end
function read(path::AbstractString)
    return open(readsac, path)
end
function write(io::IO, frame::WaveFrame) end
function write(path::AbstractString, frame::WaveFrame)
    open(path, "w") do io
        return write(io, frame)
    end
    return nothing
end

function emptyheader() end

function standardname(frame::WaveFrame; standard::AbstractString="iris")
    if standard == ""
        return @sprintf("%04d.%03d")
    else
    end
end
end
