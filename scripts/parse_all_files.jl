using Pkg; Pkg.activate(".")
using ScatteringI1
using DelimitedFiles

h5file = "data/isospin1.hdf5"
path   = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"
path   = "/home/fabian/Dokumente/Physics/Data/DataVSC/measurements/"
info   = readdlm("input/input_files.csv",',',skipstart=1)
single_file = true

for (name,dir,file,run) in eachrow(info)

    file = joinpath(path,dir,name,file)
    ens  = joinpath(name,run)

    @show file
    
    if single_file 
        h5file == "test.hdf5" && isfile(h5file) && rm(h5file)
        isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true)
    else
        dir = joinpath(dirname(h5file),"ensembles")
        ispath(dir) || mkpath(dir) 
        f = joinpath(dir,name*".hdf5")
        isospin1_to_hdf5(file,f;ensemble=run,setup=true)
    end
end
