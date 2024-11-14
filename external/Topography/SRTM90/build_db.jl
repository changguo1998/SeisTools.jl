using ArchGDAL, Mmap

exit(0)

function _read_asc_header(io::IO)
    seek(io, 0)

end

mkpath("buffer")

zip_files = readdir("raw_data")

for zf in zip_files
    cmd = Cmd(["unzip", joinpath("raw_data", zf), "-d", "buffer"])
    run(cmd)
end

ascfiles = readdir("buffer")
for af in ascfiles
    local ds = ArchGDAL.readraster(joinpath("buffer", af))
    local dat = ds[:,:,1]
    local el = zeros(size(dat))
    el .= dat
    el[dat .== -32768] .= NaN
end
