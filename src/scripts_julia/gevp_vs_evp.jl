using Pkg; Pkg.activate("src/src_jl")
using HDF5
using LatticeUtils
using Plots

file = "data_assets/isospin1_fitresults_evp.hdf5"
runs = keys(h5open(file))

println()
println("Comparison: GEVP vs EVP")
for r in runs
    h5dset = h5open(file)[r]
    ensemble_ids = filter(contains("evp"),keys(h5dset))
    println("-"^88)
    for i in eachindex(ensemble_ids)
        p = "p(0,0,1)"
        irrep = "A1"
        E  = h5dset[ensemble_ids[i]]["$p/$irrep/E"][]
        ΔE = h5dset[ensemble_ids[i]]["$p/$irrep/Delta_E"][]
        println("ens = $(rpad(joinpath(r,p,irrep),36)) $(rpad(ensemble_ids[i],10)): E0 = $(rpad(errorstring(E[1], ΔE[1]),13))E1 = $(errorstring(E[2], ΔE[2]))")
    end
end