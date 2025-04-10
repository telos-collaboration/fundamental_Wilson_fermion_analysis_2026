using DelimitedFiles
using HDF5
using ScatteringI1

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
function table_yannick(h5fitresults,infvolumefile,outfile)

    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    csv_io = open(outfile,"w")

    fid = h5open(h5fitresults)
    rex = r"Lt([0-9]+)Ls([0-9]+)beta([0-9]+.[0-9]+)m(-[0-9]+.[0-9]+)"
    inf_vol  = readdlm(infvolumefile,',',skipstart=1)

    for ens in keys(fid)

        T, L, β, m0 = parse.(Float64,match(rex,ens).captures)
        ind = findfirst(i -> [β,m0] == inf_vol[i,1:2],1:first(size(inf_vol)))
        mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]
    
        for p in keys(fid[ens])
            px, py, pz = ScatteringI1._parse_momentum(p)
            E1, ΔE1 = read(fid[ens][p],"E0")[1], read(fid[ens][p],"Delta_E0")[1]
            E2, ΔE2 = read(fid[ens][p],"E1")[1], read(fid[ens][p],"Delta_E1")[1]
            println(csv_io,"$(Int(L)) $px $py $pz 1 $E1 $ΔE1 $ΔE1 $mπ $mρ")
            println(csv_io,"$(Int(L)) $px $py $pz 2 $E2 $ΔE2 $ΔE2 $mπ $mρ")
        end
    end
    close(csv_io)
end
