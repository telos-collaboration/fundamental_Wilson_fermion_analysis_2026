using Pkg
Pkg.activate(".")
using ScatteringI1

file ="/home/fabian/Downloads/out_scattering_I1"

@show ScatteringI1._count_labels(file)
@show ScatteringI1._sources(file)
@show ScatteringI1._label_list(file)
Re, Im = parse_isospin_one(file)