using Pkg; Pkg.activate(".",io=devnull)
using ScatteringI1
using HDF5
using LatticeUtils
using Plots
using LaTeXStrings
using ProgressMeter
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function write_all_eigenvalues(infile,outfile; t0, deriv, maxhits=typemax(Int), plotting=true, average_equivalent_momenta=true, plotpath=joinpath("./plots/eigenvalues/","t0$(t0)"*(deriv ? "_deriv" : "")))
    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)
    plotting && ispath(plotpath) || mkpath(plotpath)

    ensembles = keys(h5dset)
    @showprogress desc="Write eigenvalues:" enabled=true for ens in ensembles
        _copy_lattice_parameters(outfile,infile;group=ens)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        for p in p_external
            p == "p(0,0,0)" && continue

            Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p;maxhits=typemax(Int),average_equivalent_momenta)
            eigvals, Δeigvals, eigvals_cov = ScatteringI1.variational_analysis(Corr;t0,deriv)
            eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)

            # Save plots of eigenvalues so that they can be visually examined for violations of convexity
            if plotting 
                T, L  = read(h5dset,joinpath(ens,"lattice"))[1:2]
                m0    = only(read(h5dset,joinpath(ens,"quarkmasses")))
                ncfg  = read(h5dset,joinpath(ens,"Nconf"))
                title = L"${%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf p = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg, t_0 = %$(t0)$"
                
                t  = deriv ? filter(!isequal(T÷2+1),1:T) : 1:T
                t1 = filter(x->!iszero(eigvals[1,x]),t)
                t2 = filter(x->!iszero(eigvals[2,x]),t)
                f  = deriv ? abs : identity
                
                plt = plot(yscale=:log10)
                plot!(plt;ylabel=L"$|C(t)|$",xlabel=L"t",title)
                plot_correlator!(plt,t,f.(eigvals[1,t1]),Δeigvals[1,t1],label="eigval #1")
                plot_correlator!(plt,t,f.(eigvals[2,t2]),Δeigvals[2,t2],label="eigval #2")
                plot!(plt,[t0]    ,seriestype="vline", color=:black, label="")
                plot!(plt,[T-t0+2],seriestype="vline", color=:black, label="")
                #annotate!(plt,[t0 + 1,T-t0-3] .+ 1,[ylims(plt)[2]/10,ylims(plt)[2]/10],[L"t_0",L"T - t_0"])
                savefig(plt,joinpath(plotpath,"$(ens)_$(p).pdf"))
                isinteractive() && display(plt)
            end

            h5write(outfile,joinpath(ens,p,"eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"cov_eigvals"),eigvals_cov)
            h5write(outfile,joinpath(ens,p,"t0"),t0)
            h5write(outfile,joinpath(ens,p,"deriv"),deriv)
            h5write(outfile,joinpath(ens,p,"average_equivalent_momenta"),average_equivalent_momenta)
        end
    end
end

outfile = "data/isospin1_eigenvalues_t0_8_deriv.hdf5"
infile  = "data/isospin1_corr.hdf5"
t0      = 8
deriv   = true
write_all_eigenvalues(infile,outfile; t0, deriv)