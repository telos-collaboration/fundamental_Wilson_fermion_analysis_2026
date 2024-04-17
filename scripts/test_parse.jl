using Pkg
Pkg.activate(".")
using ScatteringI1

file ="/home/fabian/Downloads/out_scattering_I1"
#parse_isospin_one(file)

@show ScatteringI1._count_labels(file)
@show ScatteringI1._sources(file)
@show ScatteringI1._label_list(file)