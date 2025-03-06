function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function write_all_eigenvalues(infile,outfile; t0, deriv, maxhits=typemax(Int), plotting=true, average_equivalent_momenta=true, use3x3=true, plotpath)
    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    if plotting
        plotname = "eigenvalues_t0$(t0)_deriv_$deriv.pdf"
        texpath  = joinpath(plotpath,"eigenvalues_tex")
        ispath(texpath)  || mkpath(texpath)
        ispath(plotpath) || mkpath(plotpath)
        isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
    end

    ensembles = keys(h5dset)
    @showprogress desc="Write eigenvalues:" enabled=true for ens in ensembles
        _copy_lattice_parameters(outfile,infile;group=ens)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        for p in p_external
            p == "p(0,0,0)" && continue

            Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits,average_equivalent_momenta)           
            eigvals, Δeigvals, eigvals_cov = ScatteringI1.variational_analysis(Corr;t0,deriv)
            eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)

            three_by_three = use3x3 && haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
            if three_by_three
                Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits,average_equivalent_momenta)
                Corr3x3[1:2,1:2,:,:] .= Corr
                eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3 = ScatteringI1.variational_analysis(Corr3x3;t0,deriv)
                eigvals_3x3, Δeigvals_3x3 = real.(eigvals_3x3), real.(Δeigvals_3x3), real.(eigvals_cov_3x3)
            end

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
                if three_by_three
                    plot_correlator!(plt,t,f.(eigvals_3x3[1,t1]),Δeigvals_3x3[1,t1],markersize=3,markershape=:rect,label="eigval #1 (3x3)")
                    plot_correlator!(plt,t,f.(eigvals_3x3[2,t2]),Δeigvals_3x3[2,t2],markersize=3,markershape=:rect,label="eigval #2 (3x3)")    
                    plot_correlator!(plt,t,f.(eigvals_3x3[3,t2]),Δeigvals_3x3[3,t2],markersize=3,markershape=:rect,label="eigval #3 (3x3)")    
                end
                plot!(plt,[t0]    ,seriestype="vline", color=:black, label="")
                plot!(plt,[T-t0+2],seriestype="vline", color=:black, label="")
                savefig(plt,"temp.pdf")
                savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p.tex") )
                append_pdf!(joinpath(plotpath,plotname),"temp.pdf",cleanup=true)
                isinteractive() && display(plt)
            end

            h5write(outfile,joinpath(ens,p,"eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"cov_eigvals"),eigvals_cov)
            h5write(outfile,joinpath(ens,p,"t0"),t0)
            h5write(outfile,joinpath(ens,p,"deriv"),deriv)
            h5write(outfile,joinpath(ens,p,"average_equivalent_momenta"),average_equivalent_momenta)
            h5write(outfile,joinpath(ens,p,"Corr2x2"),Corr)
            if three_by_three
                h5write(outfile,joinpath(ens,p,"Corr3x3"),Corr3x3)
            end
        end
    end
end