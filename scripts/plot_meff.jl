using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
gr(frame=:box,markersize=5,markeralpha=0.75)
include("utils.jl")

hdf5file = "data/isospin1.hdf5"
h5dset = h5open(hdf5file)

p   = 1
p1  = "(0,0,$p)"
ens = "Lt32Ls16beta7.2m1-0.794m2-0.794"
ens = "Lt32Ls16beta7.05m1-0.85m2-0.85"
ens = "Lt36Ls16beta7.2m1-0.76m2-0.76"
ens = "Lt48Ls16beta7.4m1-0.75m2-0.75"
ens = "Lt48Ls16beta7.4m1-0.74m2-0.74"
ens = "Lt32Ls16beta6.9m1-0.92m2-0.92"
ens = "Lt32Ls24beta6.9m1-0.92m2-0.92"

T, L = h5dset["$ens/lattice"][1:2]
title = L"$32 \times 24^3, \beta=6.9, m_0^f=-0.92, \mathbf p = (0,0,1):$ Effective masses"

Corrπ0, Corrρ0 = correlatorsp000(h5dset,ens;p)
Corrπ , Corrρ , CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp001(h5dset,ens;p)

N, nhits, T = size(CorrD1)
L3, L6 = L^3, L^6

Corr2π = dropdims(mean(CorrD1/L6 - CorrD2/L6 + 2CorrR1/L3 - 2CorrR3/L3,dims=2),dims=2) # Isospin - 1
sign   = +1
Corr2π = correlator_derivative(Corr2π;t_dim=2)
sign   = -1

Corr1π0 = dropdims(mean(Corrπ0/L3,dims=2),dims=2)
Corr1π  = dropdims(mean(Corrπ/L3 ,dims=2),dims=2)
Corr1ρ  = dropdims(mean(Corrρ/L3 ,dims=2),dims=2)
Corr1ρ0 = dropdims(mean(Corrρ0/L3 ,dims=2),dims=2)
m2π, Δm2π   = implicit_meff_jackknife(Corr2π';sign)
m1π, Δm1π   = implicit_meff_jackknife(Corr1π')
m1ρ, Δm1ρ   = implicit_meff_jackknife(Corr1ρ')
m1π0, Δm1π0 = implicit_meff_jackknife(Corr1π0')
m1ρ0, Δm1ρ0 = implicit_meff_jackknife(Corr1ρ0')

plt = plot(;title, ylabel=L"m_{\rm eff}", xlabel=L"t")
#scatter!(plt,m2π     ,yerr=Δm2π      ,marker=:circ ,label="pipi")
scatter!(plt,m1ρ     ,yerr=Δm1ρ      ,marker=:rect ,label="rho(p)")
scatter!(plt,m1π+m1π0,yerr=Δm1π+Δm1π0,marker=:pent ,label="pi(0) + pi(p)")
plot!(plt,ylims=(0.5,1.5),xlims=(1.5,13),xticks=2:2:T)

save=false
if save
    path = "plots/$hdf5file/"
    ispath(path) || mkpath(path)
    savefig(plt,joinpath(path,"effctive_masses.pdf"))
end
display(plt)