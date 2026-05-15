using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using HDF5: h5open
using DelimitedFiles: readdlm, writedlm

function p_label(p0)
    rx = r"p\(([0-9]+),([0-9]+),([0-9]+)\)"
    sx = s"(\1\2\3)"
    return replace(p0,rx=>sx)
end
function all_runs_table(h5file,outfile;ensembles_list=nothing)
    fid = h5open(h5file)
    ensembles = keys(fid)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    only_ens = isnothing(ensembles_list) ? nothing : vec(readdlm(ensembles_list))

    io = open(outfile,"w")
    for ens in ensembles
        for run in keys(fid[ens])
            if !isnothing(only_ens) && ens ∉ only_ens
                continue
            end
            rid  = fid[ens][run]
            T, L = read(rid,"lattice")[1:2]
            beta = read(rid,"beta")
            m0   = read(rid,"quarkmasses")[1]
            Ncfg = read(rid,"Nconf")
            Nsrc = read(rid,"Nsrc")
            p_etxernal = filter(p->p!="p(0,0,0)",read(rid,"p_external"))
            p_string = prod(p_label, p_etxernal) 
            write(io,"$T,$L,$beta,$m0,$Ncfg,$run,$Nsrc,$p_string\n")
        end
    end
    close(io)
    data = readdlm(outfile,',')
    data = sortslices(data,dims=1,by=x->(x[3],-x[4],x[2],-x[7]))
    writedlm(outfile,data,',')
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the parsed data"
        required = true
        "--outfile"
        help = "CSV output containing the list of all parsed runs"
        required = true
        "--ensembles_list"
        help = "CSV file containing the ensembles to analyse "
        default = nothing
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    ensembles_list = args["ensembles_list"]
    all_runs_table(args["h5file"],args["outfile"];ensembles_list)
end
main()
