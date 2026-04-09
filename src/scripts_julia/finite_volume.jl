using Pkg; Pkg.activate("src/src_jl")
using HDF5
using Plots
using LaTeXStrings
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,size=(300,300),titlefontsize=11)

file = "data_assets/isospin1_fitresults.hdf5"
h5id = h5open(file)
dir = "tmp/finite_volume/"

ens = keys(h5id)
rx  = r"beta(?<b>[0-9]+.[0-9]+)m-(?<m>[0-9]+.[0-9]+)"
masses = unique(getindex.(match.(rx,ens),"m"))

io = open("tmp/fv.csv","w")
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
    plt = scatter(inv.(L),E,yerr=ΔE,label=L"\pi")
    scatter!(inv.(L),Eρ,yerr=ΔEρ,label=L"\rho")
    plot!(plt,xlabel=L"1/L",ylabel=L"E_{\rm PS}")
    plot!(plt,xlims=[0,maximum(inv,0.8L)])
    ispath(dir) || mkpath(dir)
    savefig(joinpath(dir,"fv_b$(β)_m$m.pdf"))
end
close(io)