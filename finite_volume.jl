using Pkg; Pkg.activate("src/src_jl")
using HDF5
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,size=(300,300),titlefontsize=11)

file = "data_assets/isospin1_fitresults.hdf5"
h5id = h5open(file)

ens = keys(h5id)
masses = unique(getindex.(match.(r"m-(?<m>[0-9]+.[0-9]+)",ens),"m"))

for m in masses
    ensembles = filter(contains(m),ens)
    L  = [ read(h5id[e],"lattice")[end] for e in ensembles ] 
    E  = [ only(read(h5id[e],"p(0,0,0)/pi/E")) for e in ensembles ] 
    ΔE = [ only(read(h5id[e],"p(0,0,0)/pi/Delta_E")) for e in ensembles ]
    plt = scatter(inv.(L),E,yerr=ΔE,label="")
    plot!(plt,xlabel=L"1/L",ylabel=L"E_{\rm PS}")
    plot!(plt,xlims=[0,maximum(inv,0.8L)])
    display(plt)
end