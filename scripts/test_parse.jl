using Pkg; Pkg.activate(".")
using ScatteringI1
h5file = "isospin1.hdf5"

path = "/home/fabian/Documents/Physics/Data/I1_HDF5/m0.92/"
file = joinpath(path,"out_scattering_I1_T32_L24_h4_p1")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="T32_L24_h4_p1")
file = joinpath(path,"out_scattering_I1_T32_L16_h4_p1")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="T32_L24_h16_p1")
file = joinpath(path,"out_scattering_I1_T32_L16_h4_p1")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="T32_L16_h4_p1")
file = joinpath(path,"out_scattering_I1_T32_L16_h4_p1_PA")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="T32_L16_h4_p1_PA")
file = joinpath(path,"out_scattering_I1_T32_L16_h2_p2")
isospin1_to_hdf5(file,h5file;pmax=2,ensemble="T32_L16_h2_p2")
file = joinpath(path,"out_scattering_I1_T32_L24_h1_p23")
isospin1_to_hdf5(file,h5file;pmax=3,ensemble="T32_L24_h16_p23")

path = "/home/fabian/Documents/Physics/Data/I1_HDF5/m0.90/"
file = joinpath(path,"out_scattering_I1_T24_L12_h4_p1")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="m0.90/T24_L12_h4_p1")
file = joinpath(path,"out_scattering_I1_T32_L12_h4_p1")
isospin1_to_hdf5(file,h5file;pmax=1,ensemble="m0.90/T32_L12_h4_p1")