using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
#gr(frame=:box, legend=:topright,legendfontsize=12)
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=11, legendfontsize=8,markersize=5,tickfontsize=12)
plotlyjs(frame=:box)
include("utils.jl")

hdf5file = "isospin1_L24_h16.hdf5"
hdf5file = "isospin1_L24_h16_p23.hdf5"
hdf5file = "isospin1_L16_PA_h4.hdf5"
hdf5file = "isospin1_L16_h4.hdf5"
hdf5file = "isospin1_T24_L12.hdf5"
hdf5file = "isospin1_T32_L12.hdf5"

h5dset = h5open(hdf5file)

p1  = "(0,0,1)"
p   = 1
ens = "E1"
T, L = h5dset["$ens/lattice"][1:2]

Corrπ0, Corrρ0, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp000(h5dset,ens)
Corrπ , Corrρ , CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens;p)

N, nhits, T = size(CorrD1)
L3, L6 = L^3, L^6

CD1, ΔCD1 = _average_correlator(CorrD1/L6)
CD2, ΔCD2 = _average_correlator(CorrD2/L6)
Cπ,  ΔCπ  = _average_correlator(Corrπ/L3)
Cρ,  ΔCρ  = _average_correlator(Corrρ/L3)
CT1, ΔCT1 = _average_correlator(CorrT1/L3)
CT2, ΔCT2 = _average_correlator(CorrT2/L3)
CR1, ΔCR1 = _average_correlator(CorrR1/L3)
CR2, ΔCR2 = _average_correlator(CorrR2/L3)
CR3, ΔCR3 = _average_correlator(CorrR3/L3)
CR4, ΔCR4 = _average_correlator(CorrR4/L3)

Corr2π = dropdims(mean(CorrD1/L6+2CorrR1/L3+CorrD2/L6+2CorrR3/L3,dims=2),dims=2)
Corr2π = correlator_derivative(Corr2π,t_dim=2)

Corr1π0 = dropdims(mean(Corrπ0/L3,dims=2),dims=2)
Corr1π  = dropdims(mean(Corrπ/L3 ,dims=2),dims=2)
Corr1ρ  = dropdims(mean(Corrρ/L3 ,dims=2),dims=2)
Corr1ρ0 = dropdims(mean(Corrρ0/L3 ,dims=2),dims=2)
m2π, Δm2π   = implicit_meff_jackknife(Corr2π',sign=-1)
m1π, Δm1π   = implicit_meff_jackknife(Corr1π')
m1π0, Δm1π0 = implicit_meff_jackknife(Corr1π0')
m1ρ, Δm1ρ   = implicit_meff_jackknife(Corr1ρ')
m1ρ0, Δm1ρ0 = implicit_meff_jackknife(Corr1ρ0')

scatter(m2π,yerr=Δm2π, label="pipi")
scatter!(m1ρ0,yerr=Δm1ρ0,label="rho")
scatter!(m1π+m1π0,yerr=Δm1π+Δm1π0,label="pi(0) + pi(p)")
plot!(ylims=(0,3))