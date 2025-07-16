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

    plotname = "meson_correlators.pdf"
    texpath  = joinpath(plotpath,"eigenvalues_tex")
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))

    @showprogress desc="Plot single meson correlators:" for ens in ensembles
 
        p0 = read(h5dset,"$ens/p_external")
        p_external = unique_momenta(p0)
        T, L = h5dset["$ens/lattice"][1:2]
        m0 = h5dset["$ens/quarkmasses"][1]
        ncfg  = read(h5dset,joinpath(ens,"Nconf"))
        
        for p in p_external
            p == "p(0,0,0)" && continue
            
            momenta = read(h5dset,joinpath(ens,p,"A1","momenta"))
            sources = read(h5dset,joinpath(ens,p,"A1","sources"))
            title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(momenta), n_{src}=%$(sources), n_{cfg}=%$ncfg"
            plt = plot(;yscale=:log10,legend=:outerright,ylabel=L"$C_(t)$",xlabel=L"t",title)
            t = 1:T 

            if haskey(h5dset,joinpath(ens,p,"B1"))
                C = read(h5dset,joinpath(ens,p,"B1","C"))
                ΔC = read(h5dset,joinpath(ens,p,"B1","Delta_C"))
                plot_correlator!(plt,t,C,ΔC,markersize=3,label=L"C_{B1}")
            end
            if haskey(h5dset,joinpath(ens,p,"E"))
                C = read(h5dset,joinpath(ens,p,"E","C"))
                ΔC = read(h5dset,joinpath(ens,p,"E","Delta_C"))
                plot_correlator!(plt,t,C,ΔC,markersize=3,label=L"C_{E}")
            end
            if haskey(h5dset,joinpath(ens,p,"Cpi"))
                C = read(h5dset,joinpath(ens,p,"Cpi"))
                ΔC = read(h5dset,joinpath(ens,p,"Delta_Cpi"))
                plot_correlator!(plt,t,C,ΔC,markersize=3,label=L"C_{\pi}")
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