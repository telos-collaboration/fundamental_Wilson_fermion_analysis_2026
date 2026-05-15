using Pkg; Pkg.activate("src/src_jl")
using HDF5
using Plots
using LaTeXStrings
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,size=(300,300),titlefontsize=11)

function finite_volume_plots_appendix(file,dir,table)
    h5id = h5open(file)

    ens = keys(h5id)
    rx  = r"beta(?<b>[0-9]+.[0-9]+)m-(?<m>[0-9]+.[0-9]+)"
    masses = unique(getindex.(match.(rx,ens),"m"))

    io = open(table,"w")
    println(io,"ens,β,m,T,L,E_PS,ΔE_PS,E_V,ΔE_V")

    for m in masses
        ensembles = filter(contains(m),ens)
        β  = only(String.(unique(getindex.(match.(rx,ensembles),"b"))))
        T  = [ read(h5id[e],"lattice")[1] for e in ensembles ] 
        L  = [ read(h5id[e],"lattice")[end] for e in ensembles ] 
        E  = [ only(read(h5id[e],"p(0,0,0)/pi/E")) for e in ensembles ] 
        ΔE = [ only(read(h5id[e],"p(0,0,0)/pi/Delta_E")) for e in ensembles ]
        Eρ = [ only(read(h5id[e],"p(0,0,0)/T1/E")) for e in ensembles ] 
        ΔEρ= [ only(read(h5id[e],"p(0,0,0)/T1/Delta_E")) for e in ensembles ]

        for i in eachindex(L)
            println(io,"$(ens[i]),$β,$m,$(T[i]),$(L[i]),$(E[i]),$(ΔE[i]),$(Eρ[i]),$(ΔEρ[i])")
        end

        # start plotting
        plt = scatter(L.* E,E,yerr=ΔE,label=L"{\rm PS}")
        scatter!(L.* E,Eρ,yerr=ΔEρ,label=L"{\rm V}")
        plot!(plt,xlabel=L"m_{\rm PS} L",ylabel=L"a E_0^{M}")
        plot!(plt,legend=:right)
        ispath(dir) || mkpath(dir)
        savefig(joinpath(dir,"fv_b$(β)_m$m.pdf"))
    end
    close(io)
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
        "--table_out"
        help = "Where to save a CSV file containing the finite volume meson energy levels"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file_in"]
    dir = args["plotpath"]
    table = args["table_out"]
    finite_volume_plots_appendix(file,dir,table)
end
main()