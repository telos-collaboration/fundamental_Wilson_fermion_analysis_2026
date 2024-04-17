using Pkg; Pkg.activate(".")
using ScatteringI1

file ="/home/fabian/Downloads/out_scattering_I1"

labels = label_list(file)
Re, Im = parse_isospin_one(file)

function _write_lattice_setup(file,h5file;h5group="")
    h5write(h5file,joinpath(h5group,"plaquette"),plaquettes(file))
    h5write(h5file,joinpath(h5group,"configurations"),confignames(file))
    h5write(h5file,joinpath(h5group,"beta"),inverse_coupling(file))
    h5write(h5file,joinpath(h5group,"lattice"),latticesize(file))
    h5write(h5file,joinpath(h5group,"quarkmasses"),fermionmasses(file))
end

# HDF5 Struktur wie mit Yannick besprochen
"""" 
    File_001
    "/pi/(000)"  => array(Nconf,Nsrc,T)
    "rho/(001)" => array(Nconf,Nsrc,T)

    File_011
    "/pi/(000)"  => array(Nconf,Nsrc,T)
    "rho/(001)" => array(Nconf,Nsrc,T)
"""
