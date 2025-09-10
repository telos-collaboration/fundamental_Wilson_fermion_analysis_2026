using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using HDF5: h5open, h5write
using ScatteringI1
using LaTeXStrings: @L_str
using Plots: gr, plot, plot!, scatter!, savefig, backend_name
using PDFmerger: append_pdf!
using DelimitedFiles: readdlm
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function _get_title(h5dset,ens,p)
    momenta  = read(h5dset,joinpath(ens,p,"A1","momenta"))
    sources  = read(h5dset,joinpath(ens,p,"A1","sources"))
    gevp     = read(h5dset,joinpath(ens,p,"A1","gevp"))
    t0       = read(h5dset,joinpath(ens,p,"A1","t0"))
    T, L  = read(h5dset,joinpath(ens,"lattice"))[1:2]
    m0    = only(read(h5dset,joinpath(ens,"quarkmasses")))
    ncfg  = read(h5dset,joinpath(ens,"Nconf"))
    if gevp
        title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg, t_0 = %$(t0)"
    else
        title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg"
    end
    return title
end

function plot_eigenvalues(file,plotname,metadata)
    h5dset = h5open(file)
    ensembles = keys(h5dset)

    ispath(dirname(plotname)) || mkpath(dirname(plotname))
    isfile(plotname) && rm(plotname)

    @showprogress desc="Plot eigenvalues:" for ens in ensembles
 
        p0           = read(h5dset,"$ens/p_external")
        p_external   = unique_momenta(p0)
        
        for p in p_external
            p == "p(0,0,0)" && continue

            # get metadate for specific momentum
            data = readdlm(metadata,',',skipstart=1)
            metadata_ind = findfirst(i -> isequal(joinpath(ens,p),joinpath(data[i,1:2]...)),1:first(size(data)))
            use3x3 = data[metadata_ind,12]
            three_by_three = haskey(h5dset[ens][p],"A1/Corr3x3") && use3x3

            eigvals  = read(h5dset,joinpath(ens,p,"A1","eigvals"))
            Δeigvals = read(h5dset,joinpath(ens,p,"A1","Delta_eigvals"))
            gevp     = read(h5dset,joinpath(ens,p,"A1","gevp"))
            deriv    = read(h5dset,joinpath(ens,p,"A1","deriv"))
            t0       = read(h5dset,joinpath(ens,p,"A1","t0"))
            T, L     = read(h5dset,joinpath(ens,"lattice"))[1:2]

            if three_by_three
                Δeigvals_3x3 = read(h5dset,joinpath(ens,p,"A1","Delta_eigvals_3x3"))
                eigvals_3x3  = read(h5dset,joinpath(ens,p,"A1","eigvals_3x3"))
            end
            
            t  = deriv ? filter(!isequal(T÷2+1),1:T) : 1:T
            t1 = filter(x->!iszero(eigvals[1,x]),t)
            t2 = filter(x->!iszero(eigvals[2,x]),t)
            f  = deriv ? abs : identity
            
            title = _get_title(h5dset,ens,p)
            plt = plot(yscale=:log10,legend=:top)
            plot!(plt;ylabel=L"$|C(t)|$",xlabel=L"t",title)
            plot_correlator!(plt,t,f.(eigvals[1,t1]),Δeigvals[1,t1],label="eigval #1")
            plot_correlator!(plt,t,f.(eigvals[2,t2]),Δeigvals[2,t2],label="eigval #2")
            if three_by_three
                plot_correlator!(plt,t,f.(eigvals_3x3[1,t1]),Δeigvals_3x3[1,t1],markersize=3,markershape=:rect,label="eigval #1 (3x3)")
                plot_correlator!(plt,t,f.(eigvals_3x3[2,t2]),Δeigvals_3x3[2,t2],markersize=3,markershape=:rect,label="eigval #2 (3x3)")    
                plot_correlator!(plt,t,f.(eigvals_3x3[3,t2]),Δeigvals_3x3[3,t2],markersize=3,markershape=:rect,label="eigval #3 (3x3)")    
            end
            if gevp
                plot!(plt,[t0]    ,seriestype="vline", color=:black, label="")
                plot!(plt,[T-t0+2],seriestype="vline", color=:black, label="")
            end
            
            savefig(plt,"temp.pdf")
            append_pdf!(plotname,"temp.pdf",cleanup=true)
        end
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_in"
        help = "HDF5 file containing the parsed data"
        required = true
        "--metadata"
        help = "CSV file containing the parameters for the variational analysis"
        required = true
        "--plotname"
        help = "PDF filename for the plot"
        default = "eigenvalues.pdf"
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plot_eigenvalues(args["h5file_in"],args["plotname"],args["metadata"])
end
main()