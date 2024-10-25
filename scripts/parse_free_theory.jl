using Pkg; Pkg.activate(".")
using ScatteringI1

dir = "/home/fabian/Downloads/ScatteringTestsUnitMatrix/Measurements/"
pmax = 1
files = filter(startswith("out_")∘basename,readdir(dir,join=true))
h5file = "/home/fabian/Downloads/free_theory_results.hdf5"
for file in files
    name = split(basename(file),['_','.'])[end-1]
    isospin1_to_hdf5(file,h5file,pmax;ensemble=name)
end

using HDF5
f = h5open(h5file)