using Pkg; Pkg.activate(".")
using ScatteringI1
h5file = "isospin1_finer_lattices.hdf5"

path = "/home/fabian/Documents/Physics/Data/DataVSC/measurements/"

names = [
    "runsSp4/Lt32Ls16beta7.05m1-0.85m2-0.85",
    "runsSp4/Lt32Ls16beta7.2m1-0.78m2-0.78",
    "runsSp4/Lt32Ls16beta7.2m1-0.794m2-0.794",
    #"runsSp4/Lt48Ls16beta7.4m1-0.74m2-0.74",
    "Pusan/Lt36Ls16beta7.2m1-0.76m2-0.76",
]

for name in names
    file = joinpath(path,"$name/out/out_scattering_I1")
    isospin1_to_hdf5(file,h5file;pmax=2,ensemble=splitpath(file)[end-2])
end
