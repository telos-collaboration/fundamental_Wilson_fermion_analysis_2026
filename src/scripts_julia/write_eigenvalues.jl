using ProgressMeter: @showprogress
using HDF5: h5open, h5write
using ScatteringI1
using LaTeXStrings: @L_str
using Plots: gr, plot, plot!, scatter!, savefig, backend_name
using PDFmerger: append_pdf!
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function _copy_lattice_parameters_eigenvalues(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        entry == "p_external" && continue
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end

function write_all_eigenvalues(infile,outfile; t0, deriv, maxhits=typemax(Int), average_equivalent_momenta, gevp, symmetrise)    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    ensembles = keys(h5dset)
    @showprogress desc="Write eigenvalues:" enabled=true for ens in ensembles
        _copy_lattice_parameters_eigenvalues(outfile,infile;group=ens)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        h5write(outfile,"$ens/p_external",p_external)
        for p in p_external
            p == "p(0,0,0)" && continue

            Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits,average_equivalent_momenta)           
            eigvals, Δeigvals, eigvals_cov = ScatteringI1.variational_analysis(Corr;t0,deriv,gevp,symmetrise)
            eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)
            meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv,gevp,symmetrise)

            three_by_three = haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
            if three_by_three
                Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits,average_equivalent_momenta)
                Corr3x3[1:2,1:2,:,:] .= Corr
                eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3 = ScatteringI1.variational_analysis(Corr3x3;t0,deriv,gevp,symmetrise)
                eigvals_3x3, Δeigvals_3x3 = real.(eigvals_3x3), real.(Δeigvals_3x3), real.(eigvals_cov_3x3)
                meff_3x3, Δmeff_3x3 = ScatteringI1.effective_masses(Corr3x3;t0,deriv,gevp,symmetrise)
            end
   
            h5write(outfile,joinpath(ens,p,"A1","eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"A1","Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"A1","cov_eigvals"),eigvals_cov)
            h5write(outfile,joinpath(ens,p,"t0"),t0)
            h5write(outfile,joinpath(ens,p,"gevp"),gevp)
            h5write(outfile,joinpath(ens,p,"deriv"),deriv)
            h5write(outfile,joinpath(ens,p,"symmetrise"),symmetrise)
            h5write(outfile,joinpath(ens,p,"average_equivalent_momenta"),average_equivalent_momenta)
            h5write(outfile,joinpath(ens,p,"momenta"),ascii(momenta))
            h5write(outfile,joinpath(ens,p,"sources"),[sources...])
            h5write(outfile,joinpath(ens,p,"Corr2x2"),Corr)
            h5write(outfile,joinpath(ens,p,"meff"),meff)
            h5write(outfile,joinpath(ens,p,"Delta_meff"),Δmeff)            
            if three_by_three
                h5write(outfile,joinpath(ens,p,"Corr3x3"),Corr3x3)
                h5write(outfile,joinpath(ens,p,"momenta_3x3"),ascii(momenta3x3))
                h5write(outfile,joinpath(ens,p,"sources_3x3"),[sources3x3...])
                h5write(outfile,joinpath(ens,p,"eigvals_3x3"),eigvals_3x3)
                h5write(outfile,joinpath(ens,p,"Delta_eigvals_3x3"),Δeigvals_3x3)
                h5write(outfile,joinpath(ens,p,"cov_eigvals_3x3"),eigvals_cov_3x3)
                h5write(outfile,joinpath(ens,p,"meff_3x3"),meff_3x3)
                h5write(outfile,joinpath(ens,p,"Delta_meff_3x3"),Δmeff_3x3)            
            end
        end
    end
end
function plot_eigenvalues(file,plotpath)
    h5dset = h5open(file)
    ensembles = keys(h5dset)
    plotting = true

    if plotting
        plotname = "eigenvalues.pdf"
        texpath  = joinpath(plotpath,"eigenvalues_tex")
        ispath(texpath)  || mkpath(texpath)
        ispath(plotpath) || mkpath(plotpath)
        isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
    end

    @showprogress desc="Plot eigenvalues:" for ens in ensembles
 
        p0           = read(h5dset,"$ens/p_external")
        p_external   = unique_momenta(p0)
        
        for p in p_external
            p == "p(0,0,0)" && continue
            
            three_by_three = use3x3 && haskey(h5dset[ens][p],"Corr3x3")

            eigvals  = read(h5dset,joinpath(ens,p,"A1","eigvals"))
            Δeigvals = read(h5dset,joinpath(ens,p,"A1","Delta_eigvals"))
            momenta  = read(h5dset,joinpath(ens,p,"momenta"))
            sources  = read(h5dset,joinpath(ens,p,"sources"))
            if three_by_three
                Δeigvals_3x3 = read(h5dset,joinpath(ens,p,"Delta_eigvals_3x3"))
                eigvals_3x3  = read(h5dset,joinpath(ens,p,"eigvals_3x3"))
                momenta_3x3  = read(h5dset,joinpath(ens,p,"momenta_3x3"))
                sources_3x3  = read(h5dset,joinpath(ens,p,"sources_3x3"))
            end

            if plotting 
                T, L  = read(h5dset,joinpath(ens,"lattice"))[1:2]
                m0    = only(read(h5dset,joinpath(ens,"quarkmasses")))
                ncfg  = read(h5dset,joinpath(ens,"Nconf"))
                if gevp
                    title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg, t_0 = %$(t0)"
                else
                    title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg"
                end
                
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
                if backend_name() == :pgfplotsx
                    ispath(texpath)  || mkpath(texpath)
                    savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p.tex") )
                end
                append_pdf!(joinpath(plotpath,plotname),"temp.pdf",cleanup=true)
                isinteractive() && display(plt)
            end

        end
    end
end