wkdir = abspath(@__DIR__)
@info "run in directory: $wkdir"
modeldirs = filter(readdir(wkdir)) do mdir
    if !isdir(joinpath(wkdir, mdir))
        flag = false
    end
    files = readdir(joinpath(wkdir, mdir))
    flag = length(files) > 2
    flag
end

@info "$(length(modeldirs)) model(s) are found:"
for m in modeldirs
    @info "    $(m)"
end

@info "Start converting..."
for m in modeldirs
    @info "    $(m) ..."
    cmd = Cmd(["julia", "build_db.jl"]; dir=joinpath(wkdir, m))
    @info "    $(m) done"
end

@info "Finish converting"
