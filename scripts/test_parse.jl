using Pkg; Pkg.activate(".")
using ScatteringI1

h5file = "data/isospin1.hdf5"
path   = "/home/fabian/Dokumente/Physics/Data/DataVSC/measurements/"
path   = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"
log    = "out/out_scattering_I1"
names  = [
    #"runsSp4/Lt24Ls14beta6.9m1-0.92m2-0.92",
    "runsSp4/Lt32Ls16beta6.9m1-0.92m2-0.92",
    "Pusan/Lt32Ls24beta6.9m1-0.92m2-0.92",
    #"runsSp4/Lt24Ls14beta7.05m1-0.85m2-0.85",
    #"runsSp4/Lt32Ls16beta7.05m1-0.85m2-0.85",
    #"runsSp4/Lt32Ls16beta7.2m1-0.78m2-0.78",
    #"runsSp4/Lt48Ls16beta7.4m1-0.74m2-0.74",
    #"Pusan/Lt36Ls16beta7.2m1-0.76m2-0.76",
    #"Pusan/Lt32Ls24beta6.9m1-0.924m2-0.924",
    #"runsSp4/Lt48Ls16beta7.4m1-0.75m2-0.75",
    #"runsSp4/Lt32Ls16beta7.2m1-0.794m2-0.794",
    #"runsSp4/Lt36Ls24beta7.05m1-0.867m2-0.867",
    #"Pusan/Lt36Ls20beta7.05m1-0.835m2-0.835",
]

single_file = true

for name in names
    file = joinpath(path,"$name/$log")
    ens  = basename(name)
    @show name

    if single_file 
        h5file == "test.hdf5" && isfile(h5file) && rm(h5file)
        isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true)
    else
        dir = joinpath(dirname(h5file),"ensembles")
        ispath(dir) || mkpath(dir) 
        f = joinpath(dir,ens*".hdf5")
        isospin1_to_hdf5(file,f;ensemble="",setup=true)
    end
end

using HDF5
f = h5open(h5file)
ensembles = keys(f)
entries   = keys(f[ensembles[3]])
println.(entries)
