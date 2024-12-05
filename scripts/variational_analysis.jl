using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
using LatticeUtils
include("utils.jl")
gr(frame=:box, legend=:topright,legendfontsize=12)

hdf5file = "data/isospin1.hdf5"
h5dset = h5open(hdf5file)
p   = 1
p1  = "(0,0,$p)"
ens = "Lt48Ls16beta7.4m1-0.75m2-0.75"
ens = "Lt48Ls16beta7.4m1-0.74m2-0.74"
ens = "Lt36Ls16beta7.2m1-0.76m2-0.76"
ens = "Lt32Ls16beta7.2m1-0.794m2-0.794"
ens = "Lt32Ls16beta7.2m1-0.78m2-0.78"
ens = "Lt32Ls16beta7.05m1-0.85m2-0.85"
ens = "Lt32Ls16beta6.9m1-0.92m2-0.92"
ens = "Lt32Ls24beta6.9m1-0.92m2-0.92"
t0  = 2

T, L = h5dset["$ens/lattice"][1:2]
title_corr = L"$%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1$: Eigenvalues from GEVP "
title_meff = L"$%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1$: Effective mass from GEVP"

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp001(h5dset,ens;p)
N, nhits, T = size(CorrD1)

function pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
    L3, L6 = L^3, L^6
    Corr2π = + CorrD1/L6 - CorrD2/L6 + 2CorrR1/L3 - 2CorrR3/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(ComplexF64,(2,2,N,nhits,T))
    corr[1,1,:,:,:] =  @. Corrρ/L3 + 0*im
    corr[1,2,:,:,:] =  @. 0        + im*(CorrT1-CorrT2)/L3
    corr[2,1,:,:,:] =  @. 0        + im*(CorrT2-CorrT1)/L3
    corr[2,2,:,:,:] =  @. Corr2π   + 0*im
    corr = dropdims(mean(corr,dims=4),dims=4)
    return corr
end
Corr2π = pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
sign   = +1
Corr   = correlator_derivative(Corr;t_dim=4)
sign   = -1

eigvals, Δeigvals = eigenvalues(Corr;t0)
eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
meff, Δmeff =  LatticeUtils.implicit_meff_jackknife(real.(eigvals_resamples);sign)

r = 1:15
#plt = plot(;title=title_meff)
scatter!(plt,meff[2,:],yerr=Δmeff[2,:],label="Eigenvalue #1")
scatter!(plt,meff[1,r],yerr=Δmeff[1,r],label="Eigenvalue #2")
plot!(plt,ylims=(0.5,1.5),xlims=(1.5,13),xticks=2:2:T)
display(plt)
savefig(plt,joinpath("gevp.pdf"))