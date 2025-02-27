using Pkg; Pkg.activate(".",io=devnull)
using ScatteringI1
using DelimitedFiles

h5file = "data/isospin1_sorted.hdf5"
path   = "/home/fabian/Dokumente/Physics/Data/DataVSC/measurements/"
path   = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"
info   = readdlm("input/input_files.csv",',',skipstart=1)
single_file = true

println("Parse correlator data from raw log:")
for (name,dir,file,run) in eachrow(info)

    file = joinpath(path,dir,name,file)
    ens  = joinpath(name,run)

    if single_file 
        isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true,sort=true,deduplicate=true)
    else
        dir = joinpath(dirname(h5file),"ensembles")
        ispath(dir) || mkpath(dir) 
        f = joinpath(dir,name*".hdf5")
        isospin1_to_hdf5(file,f;ensemble=run,setup=true,sort=true,deduplicate=true)
    end
end
