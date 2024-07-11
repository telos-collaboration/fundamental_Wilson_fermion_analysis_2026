using Pkg; Pkg.activate(".")
using ScatteringI1

#pmax = 1
#file = "/home/fabian/Downloads/out_scattering_I1"
#h5file = "isospin1_h4_raw.hdf5"
#isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

#pmax = 1
#file ="/home/fabian/Downloads/out_scattering_I1_h16"
#h5file = "isospin1_h16_raw.hdf5"
#isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

pmax = 1
file ="/home/fabian/Downloads/out_scattering_I1_h4"
h5file = "isospin1_h4_L16_raw.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")
