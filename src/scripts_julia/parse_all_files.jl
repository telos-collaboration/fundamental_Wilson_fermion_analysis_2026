using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using DelimitedFiles: readdlm
using ScatteringI1: isospin1_to_hdf5

function parse_all_files(path,h5file,inputfiles)
    println("Parse correlator data from raw log:")
    info  = readdlm(inputfiles,',',skipstart=1)
    for (name,file,run) in eachrow(info)

        file = joinpath(path,file)
        ens  = joinpath(name,run)

        dir = dirname(h5file)
        ispath(dir) || mkpath(dir) 
        isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true,sort=true,deduplicate=true)
    end
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the parsed data"
        required = true
        "--inputfiles"
        help = "CSV file containing the files to be parsed and metadata"
        required = true
        "--path"
        help = "Path to the root directory of all files"
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    parse_all_files(args["path"],args["h5file"],args["inputfiles"])
end
main()
