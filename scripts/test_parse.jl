using Pkg; Pkg.activate(".")
using ScatteringI1

file ="/home/fabian/Downloads/out_scattering_I1"
h5file = "isospin1_raw.hdf5"
isospin1_to_hdf5(file,h5file,ensemble="E1")
