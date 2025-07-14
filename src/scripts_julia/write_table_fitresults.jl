using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using HDF5
using DelimitedFiles: readdlm, writedlm
using LatticeUtils: errorstring

function fitresult_table(h5file,outfile)
    fid = h5open(h5file)
    ensembles = keys(fid)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))

    io = open(outfile,"w")
    write(io,"name;momentum;groundstate;scatterstate\n")
    for ens in ensembles
        T,L  = read(fid[ens],"lattice")[1:2]
        moms = filter(startswith("p"),keys(fid[ens]))
        for p in moms
            E0 = fid[ens][p]["A1"]["E0"][1]
            E1 = fid[ens][p]["A1"]["E1"][1]
            ΔE0 = fid[ens][p]["A1"]["Delta_E0"][1]
            ΔE1 = fid[ens][p]["A1"]["Delta_E1"][1]
            write(io,"$ens;$p;$(errorstring(E0,ΔE0));$(errorstring(E1,ΔE1))\n")
        end
    end
    close(io)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the corrfitter results"
        required = true
        "--outfile"
        help = "CSV tablew containing the fitted energioes"
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    fitresult_table(args["h5file"],args["outfile"])
end
main()
