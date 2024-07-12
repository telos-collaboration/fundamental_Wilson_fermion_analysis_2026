using Pkg; Pkg.activate(".")
using ScatteringI1

pmax = 1
file = "/home/fabian/Downloads/I1_HDF5/out_scattering_I1_L24_h4"
h5file = "isospin1_L24_h4.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

pmax = 1
file ="/home/fabian/Downloads/I1_HDF5/out_scattering_I1_L24_h4"
h5file = "isospin1_L24_h16.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

pmax = 3
file ="/home/fabian/Downloads/I1_HDF5/out_scattering_I1_L24_h1_p23"
h5file = "isospin1_L24_h16_p23.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

pmax = 1
file = "/home/fabian/Downloads/I1_HDF5/out_scattering_I1_L16_h4"
h5file = "isospin1_L16_h4.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")

pmax = 2
file = "/home/fabian/Downloads/I1_HDF5/out_scattering_I1_L16_h2_p2"
h5file = "isospin1_L16_h2_p2.hdf5"
isospin1_to_hdf5(file,h5file,pmax;ensemble="E1")