using Pkg; Pkg.activate("scattering/src/src_jl")
Pkg.instantiate()
using HiRepParsing
using ArgParse
using HDF5

# This script parses the log files in the directory 'dir', and saves them as an hdf5-file 
# in the location provided by 'h5file'.

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`
function main(h5file,files;setup=true,filter_channels=false,channels=nothing)
    for file in files
        # only try parsing if the filesize is non-vanishing
        if filesize(file) > 0
            ensemble = match(r"Lt[0-9]+Ls[0-9]+beta[0-9]+.[0-9]+mas-[0-9]+.[0-9]+FUN",file).match
            smearing_regex = r"source_N[0-9]+_sink_N[0-9]+"
            writehdf5_spectrum_with_regexp(file,h5file,smearing_regex;mixed_rep=false,h5group=ensemble,setup,filter_channels,channels,sort=true,deduplicate=true)
        end
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the parsed data"
        required = true
        "files"
        help = "Log file(s) to be parsed"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    h5file = args["h5file"]
    files = args["files"]
    main(h5file,files)
end
main()


