using Pkg; Pkg.activate("src/src_jl")
using HDF5
using LatticeUtils
using Plots

file = "data_assets/isospin1_fitresults_evp.hdf5"
runs = keys(h5open(file))

h5dset = h5open(file)[first(runs)]
ensemble_ids = filter(contains("evp"),keys(h5dset))

for i in eachindex(ensemble_ids)
    E  = h5dset[ensemble_ids[i]]["p(0,0,1)/A1/E"][]
    ΔE = h5dset[ensemble_ids[i]]["p(0,0,1)/A1/Delta_E"][]
    @show ensemble_ids[i], errorstring(E[1], ΔE[1]), errorstring(E[2], ΔE[2])
end