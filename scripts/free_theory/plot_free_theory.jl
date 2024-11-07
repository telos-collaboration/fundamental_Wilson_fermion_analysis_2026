using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
include("../utils.jl")
include("analytic.jl")
gr(frame=:box,legendfontsize=12)

hdf5file = expanduser("~/Downloads/free_theory_results_v2.hdf5")
h5dset = h5open(hdf5file)

ens = "T16L8"
T,L = h5dset["$ens/lattice"][1:2]
m = -0.6
p = 1
P_phys = [0,0,0,1*2π/L]

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens;p=1)
average_sources(C) = dropdims(mean(C,dims=2),dims=(1,2))

Cπ = average_sources(Corrπ)
C  = analytic_free_pion(T,L,m,P_phys)

plt = plot(legend=:top)
scatter!(plt,abs.(Cπ),yscale=:log10,label="π (measured)")
scatter!(plt,real.(C),yscale=:log10,label="π (analytic)")
