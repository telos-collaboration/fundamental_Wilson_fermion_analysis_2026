using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using HDF5: h5open, h5write
using ScatteringI1
using LaTeXStrings: @L_str
using Plots: gr, plot, plot!, scatter!, savefig, backend_name
using PDFmerger: append_pdf!
using Statistics: mean, std
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function plot_correlation_matrix_elements(file,plotpath)
    h5dset = h5open(file)
    ensembles = keys(h5dset)

    plotname = "correlation_matrix_elements.pdf"
    texpath  = joinpath(plotpath,"eigenvalues_tex")
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))

    @showprogress desc="Plot correlation matrix elements:" for ens in ensembles
 
        p0         = read(h5dset,"$ens/p_external")
        p_external = unique_momenta(p0)
        
        for p in p_external
            p == "p(0,0,0)" && continue
            
            three_by_three = haskey(h5dset[ens][p],"A1/Corr3x3")
            T, L  = read(h5dset,joinpath(ens,"lattice"))[1:2]
            m0    = only(read(h5dset,joinpath(ens,"quarkmasses")))
            ncfg  = read(h5dset,joinpath(ens,"Nconf"))

            Corr = read(h5dset,joinpath(ens,p,"A1","Corr2x2"))
            momenta = read(h5dset,joinpath(ens,p,"A1","momenta"))
            sources = read(h5dset,joinpath(ens,p,"A1","sources"))
            title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg"
            if three_by_three
                Corr = read(h5dset,joinpath(ens,p,"A1","Corr3x3"))
                momenta_3x3 = read(h5dset,joinpath(ens,p,"A1","momenta_3x3"))
                sources_3x3 = read(h5dset,joinpath(ens,p,"A1","sources_3x3"))
                title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources)+(%$(sources_3x3)), n_{cfg}=%$ncfg"
            end

            N = size(Corr,3)
            C = dropdims(mean(Corr,dims=3),dims=3)
            ΔC = dropdims(std(Corr,dims=3),dims=3) ./ sqrt(N)

            t = 1:T
            plt = plot(yscale=:log10,legend=:outerright)
            plot!(plt;ylabel=L"$|C_{ij}(t)|$",xlabel=L"t",title)
            for i in axes(C,1), j in axes(C,2)
                plot_correlator!(plt,t,abs.(C[i,j,:]),abs.(ΔC[i,j,:]),markersize=3,label=L"C_{%$i%$j}")    
            end
            
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
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_in"
        help = "HDF5 file containing the parsed data"
        required = true
        "--plotpath"
        help = "HDF5 output file containing the correlation matrices"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plot_correlation_matrix_elements(args["h5file_in"],args["plotpath"])
end
main()