using Pkg; Pkg.activate("scattering/src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using HDF5
using DelimitedFiles: readdlm, writedlm
using LatticeUtils: errorstring

function _get_fit_engery_str(hdset,label)
    if haskey(hdset,label)
        E = hdset[label]["E"][1]
        ΔE = hdset[label]["Delta_E"][1]
        return errorstring(E,ΔE)
    else
        return "-"
    end
end

function fitresult_table(h5file,outfile)
    fid = h5open(h5file)
    ensembles = keys(fid)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))

    io = open(outfile,"w")
    write(io,"name;momentum;pipi_groundstate;pipi_scatterstate;pi;rho_T1;rho_E;rho_B1\n")
    for ens in ensembles
        moms = filter(startswith("p"),keys(fid[ens]))
        for p in moms
            if p == "p(0,0,0)"
                str0 = "-"
                str1 = "-"
            else
                E0 = fid[ens][p]["A1"]["E"][1]
                E1 = fid[ens][p]["A1"]["E"][2]
                ΔE0 = fid[ens][p]["A1"]["Delta_E"][1]
                ΔE1 = fid[ens][p]["A1"]["Delta_E"][2]
                str0 = errorstring(E0,ΔE0)
                str1 = errorstring(E1,ΔE1)
            end
            str_π = _get_fit_engery_str(fid[ens][p],"pi")
            str_ρE = _get_fit_engery_str(fid[ens][p],"E")
            str_ρB1 = _get_fit_engery_str(fid[ens][p],"B1")
            str_ρT1 = _get_fit_engery_str(fid[ens][p],"T1")
            write(io,"$ens;$p;$str0;$str1;$str_π;$str_ρT1;$str_ρE;$str_ρB1\n")
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
