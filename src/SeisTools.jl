module SeisTools

# * modules
# include("CommonFunction.jl")
include("DataProcess.jl")
include("Geodesy.jl")
include("QualityControl.jl")
include("SAC.jl")
include("SACPZ.jl")
include("SEGY.jl")
include("Source.jl")
include("Topography.jl")
include("VelocityModel.jl")

function download_external_data()
    external_path = abspath(@__DIR__, "..", "external")
    cmd1 = Cmd(Cmd(["julia", "download_model.jl"]); dir = joinpath(external_path, "VelocityModel"))
    run(cmd1)
end

function build_external_database()
    external_path = abspath(@__DIR__, "..", "external")
    cmd1 = Cmd(Cmd(["julia", "convert_format.jl"]); dir = joinpath(external_path, "VelocityModel"))
    run(cmd1)
    cmd2 = Cmd(Cmd(["julia", "convert_format.jl"]); dir = joinpath(external_path, "Topography"))
    run(cmd2)
end

function install_path()
    return abspath(@__DIR__, "..")
end

end # module
