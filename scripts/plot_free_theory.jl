using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
include("utils.jl")
gr(frame=:box, legend=:topright,legendfontsize=12)

hdf5file = "/home/fabian/Downloads/free_theory_results.hdf5"
h5dset = h5open(hdf5file)

p1  = "(1,1,1)"
p   = 1
ens = "T32L4"
T,L = h5dset["$ens/lattice"][1:2]

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,"T32L4";p)
CD2 = dropdims(mean(CorrD2,dims=2),dims=(1,2))
scatter(CD2,yscale=:log10)
