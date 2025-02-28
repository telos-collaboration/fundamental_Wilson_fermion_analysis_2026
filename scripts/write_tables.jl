using DelimitedFiles
using HDF5

function p_label(p0)
    rx = r"p\(([0-9]+),([0-9]+),([0-9]+)\)"
    sx = s"(\1\2\3)"
    return replace(p0,rx=>sx)
end
function all_runs_table(h5file,outfile)
    fid = h5open(h5file)
    ensembles = keys(fid)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))

    io = open(outfile,"w")
    for ens in ensembles
        for run in keys(fid[ens])
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

h5file  = "output/data/isospin1_sorted.hdf5"
outfile = "output/tables/all_runs.csv"
all_runs_table(h5file,outfile)
