using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
gr(frame=:box)
include("utils.jl")

hdf5file = "isospin1_finer_lattices.hdf5"
h5dset = h5open(hdf5file)

p   = 1
p1  = "(0,0,$p)"
ens = "Lt32Ls16beta7.2m1-0.794m2-0.794"
ens = "Lt32Ls16beta7.05m1-0.85m2-0.85"
ens = "Lt36Ls16beta7.2m1-0.76m2-0.76"

T, L = h5dset["$ens/lattice"][1:2]
title = L"$%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1$: Effective masses"

Corrπ0, Corrρ0 = correlatorsp000(h5dset,ens;p)
Corrπ , Corrρ , CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp001(h5dset,ens;p)

N, nhits, T = size(CorrD1)
L3, L6 = L^3, L^6

Corr2π = dropdims(mean(CorrD1/L6 - CorrD2/L6 + 2CorrR1/L3 - 2CorrR3/L3,dims=2),dims=2) # Isospin - 1

Corr1π0 = dropdims(mean(Corrπ0/L3,dims=2),dims=2)
Corr1π  = dropdims(mean(Corrπ/L3 ,dims=2),dims=2)
Corr1ρ  = dropdims(mean(Corrρ/L3 ,dims=2),dims=2)
Corr1ρ0 = dropdims(mean(Corrρ0/L3 ,dims=2),dims=2)
m2π, Δm2π   = implicit_meff_jackknife(Corr2π')
m1π, Δm1π   = implicit_meff_jackknife(Corr1π')
m1ρ, Δm1ρ   = implicit_meff_jackknife(Corr1ρ')
m1π0, Δm1π0 = implicit_meff_jackknife(Corr1π0')
m1ρ0, Δm1ρ0 = implicit_meff_jackknife(Corr1ρ0')

plt = plot(;title, ylabel=L"m_{\rm eff}", xlabel=L"t")
scatter!(plt,m2π,yerr=Δm2π, label="pipi")
scatter!(plt,m1ρ0,yerr=Δm1ρ0,label="rho(0)")
scatter!(plt,m1ρ ,yerr=Δm1ρ ,label="rho(p)")
scatter!(plt,m1π+m1π0,yerr=Δm1π+Δm1π0,label="pi(0) + pi(p)")
plot!(plt,xlims=(0,T÷2+0.5),xticks=1:2:T)

save=true
if save
    path = "plots/$hdf5file/"
    ispath(path) || mkpath(path)
    savefig(plt,joinpath(path,"effctive_masses.pdf"))
end
display(plt)