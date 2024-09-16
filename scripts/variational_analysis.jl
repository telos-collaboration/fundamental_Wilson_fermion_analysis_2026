using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
include("utils.jl")

hdf5file = "isospin1_L24_h16_p23.hdf5"
hdf5file = "isospin1_L24_h16.hdf5"
hdf5file = "isospin1_L24_h16.hdf5"
hdf5file = "isospin1_L16_h4.hdf5"
hdf5file = "isospin1_L16_PA_h4.hdf5"
h5dset = h5open(hdf5file)

p1  = "(0,0,1)"
p   = 1
ens = "E1"
T, L = h5dset["$ens/lattice"][1:2]

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens;p)
N, nhits, T = size(CorrD1)

function pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
    L3, L6 = L^3, L^6
    Corr2π = CorrD1/L6+2CorrR1/L3+CorrD2/L6+2CorrR3/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(2,2,N,nhits,T)
    corr[1,1,:,:,:] = Corrρ/L3
    corr[1,2,:,:,:] = CorrT1/L3
    corr[2,1,:,:,:] = CorrT2/L3
    corr[2,2,:,:,:] = Corr2π
    corr = dropdims(mean(corr,dims=4),dims=4)
    return corr
end
Corr2π = pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
Corr   = correlator_derivative(Corr,t_dim=4)
eigvals, Δeigvals = eigenvalues(Corr)
eigvals_resamples = eigenvalues_jackknife_samples(Corr)
meff, Δmeff =  meff_from_jackknife(eigvals_resamples;sign=+1,swap=nothing)

plotlyjs()
plt = plot(yscale=:log10)
scatter!(plt,eigvals[2,:],yerr=Δeigvals[2,:])
scatter!(plt,abs.(eigvals[1,:]),yerr=Δeigvals[1,:])

plt_meff = plot()
scatter!(plt_meff,meff[2,:],yerr=Δmeff[2,:])
scatter!(plt_meff,meff[1,:],yerr=Δmeff[1,:])
