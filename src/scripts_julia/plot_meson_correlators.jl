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

function _read_plot_correlator(plt,h5dset,ens,p,irrep;kws...)
    T = h5dset["$ens/lattice"][1]
    t = 1:T
    if haskey(h5dset,joinpath(ens,p)) && haskey(h5dset,joinpath(ens,p,irrep))
        C = read(h5dset,joinpath(ens,p,irrep,"C"))
        ΔC = read(h5dset,joinpath(ens,p,irrep,"Delta_C"))
        t = 1:length(C)
        plot_correlator!(plt,t,C ./ C[1] ,ΔC./ C[1],markersize=3;kws...)
    end
end
function _read_plot_fitresults(plt,fitres,ens,p,irrep;kws...)
    if haskey(fitres,joinpath(ens,p)) && haskey(fitres,joinpath(ens,p,irrep))
        t = read(fitres,joinpath(ens,p,irrep,"tfit")) .+ 1
        C = read(fitres,joinpath(ens,p,irrep,"fit"))
        ΔC = read(fitres,joinpath(ens,p,irrep,"Delta_fit"))
        plot_correlator!(plt,t,C,ΔC,markersize=3; type=:ribbon, kws...)
    end
end
function plot_meson_correlators(file,plotpath,fitresults)
    h5dset = h5open(file)
    ensembles = keys(h5dset)
    fitres = h5open(fitresults)

    plotname = "meson_correlators.pdf"
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))

    irreps =  ["B1", "E", "T1", "pi"]
    irrep_labels =  [L"C_\rho^{B1}", L"C_\rho^{E}", L"C_\rho^{T1}", L"C_{\pi}"]

    @showprogress desc="Plot single meson correlators:" for ens in ensembles
 
        p0 = read(h5dset,"$ens/p_external")
        p_external = unique_momenta(p0)
        T, L = h5dset["$ens/lattice"][1:2]
        m0 = h5dset["$ens/quarkmasses"][1]
        ncfg  = read(h5dset,joinpath(ens,"Nconf"))
        
        for p in p_external
            
            title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(p[2:end])"
            plt = plot(;yscale=:log10,legend=:outerright,ylabel=L"$C_(t)$",xlabel=L"t",title)

            for (irrep,label) in zip(irreps, irrep_labels)
                _read_plot_correlator(plt,h5dset,ens,p,irrep;label)
                _read_plot_fitresults(plt,fitres,ens,p,irrep;label="")
            end
            
            savefig(plt,"temp.pdf")
            append_pdf!(joinpath(plotpath,plotname),"temp.pdf",cleanup=true)
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
        "--fitresults"
        help = "CSV file containing the parameters for the variational analysis"
        default = ""
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plot_meson_correlators(args["h5file_in"],args["plotpath"],args["fitresults"])
end
main()