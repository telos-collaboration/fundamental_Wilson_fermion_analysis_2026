using Pkg
Pkg.activate(".")
using ScatteringI1

file ="/home/fabian/Downloads/out_scattering_I1"

@show ScatteringI1._count_labels(file)
@show ScatteringI1._sources(file)
@show ScatteringI1._label_list(file)
Re, Im = _parse_isospin_one(file)



# HDF5 Struktur wie mit Yannick besprochen
"""" 
    File_001
    "/pi/(000)"  => array(Nconf,Nsrc,T)
    "rho/(001)" => array(Nconf,Nsrc,T)

    File_011
    "/pi/(000)"  => array(Nconf,Nsrc,T)
    "rho/(001)" => array(Nconf,Nsrc,T)
"""
