using Pkg; Pkg.activate(".")
using ScatteringI1

dir    = expanduser("~/Downloads/Measurements/")
h5file = expanduser("~/Downloads/free_theory_results_v2.hdf5")
pmax   = 1

files = filter(startswith("out_")∘basename,readdir(dir,join=true))
for file in files
    name = split(basename(file),['_','.'])[end-1]
    isospin1_to_hdf5(file,h5file,pmax;ensemble=name)
end

using HDF5
f = h5open(h5file)
