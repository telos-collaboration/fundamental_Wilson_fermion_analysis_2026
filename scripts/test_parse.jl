using Pkg; Pkg.activate(".")
using ScatteringI1

h5file = "data/isospin1.hdf5"
path   = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"
names  = [
    "runsSp4/Lt32Ls16beta6.9m1-0.92m2-0.92",
    "runsSp4/Lt32Ls16beta7.05m1-0.85m2-0.85",
    "runsSp4/Lt32Ls16beta7.2m1-0.78m2-0.78",
    "runsSp4/Lt32Ls16beta7.2m1-0.794m2-0.794",
    "runsSp4/Lt48Ls16beta7.4m1-0.74m2-0.74",
    "runsSp4/Lt48Ls16beta7.4m1-0.75m2-0.75",
    "Pusan/Lt32Ls24beta6.9m1-0.92m2-0.92",
    "Pusan/Lt36Ls16beta7.2m1-0.76m2-0.76",
]

for name in names
    file = joinpath(path,"$name/out/out_scattering_I1")
    ens  = basename(name)
    h5file == "test.hdf5" && isfile(h5file) && rm(h5file)
    isospin1_to_hdf5(file,h5file;ensemble=ens,setup=true)
end
