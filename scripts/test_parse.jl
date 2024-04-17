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
function _splitlabel(label)
    @show label, contains(label,"p(") 
    label == "pi" && (return label, "p(0,0,0)")
    if contains(label,"p0") 
        l = replace(label,"_p0"=>"")
        return l, "p(0,0,0)"
    elseif contains(label,"p(")
        f = first(findfirst("_p(",label))
        n = findlast(')',label)
        l = label[1:f-1]*label[n+1:end]
        p = label[f+1:n]
        return l, p
    end
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
