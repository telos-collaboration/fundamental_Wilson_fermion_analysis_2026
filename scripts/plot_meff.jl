using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
using LatticeUtils
gr(frame=:box,markersize=5,markeralpha=0.75)
include("read_rhopipi_diagrams.jl")

hdf5file = "data/isospin1.hdf5"
h5dset = h5open(hdf5file)

p   = 1
p1  = "(0,0,$p)"
ens = "Lt32Ls24beta6.9m1-0.92m2-0.92"
T, L = h5dset["$ens/lattice"][1:2]

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
m2π, Δm2π   = log_meff(Corr2π')
m1π, Δm1π   = log_meff(Corr1π')
m1ρ, Δm1ρ   = log_meff(Corr1ρ')
m1π0, Δm1π0 = log_meff(Corr1π0')
m1ρ0, Δm1ρ0 = log_meff(Corr1ρ0')

plt = plot(;ylabel=L"m_{\rm eff}", xlabel=L"t")
scatter!(plt,m1π     ,yerr=Δm1π      ,marker=:circ ,label="pi(0)")
scatter!(plt,m2π     ,yerr=Δm2π      ,marker=:circ ,label="pipi")
scatter!(plt,m1ρ     ,yerr=Δm1ρ      ,marker=:rect ,label="rho(p)")
plot!(plt,ylims=(0.0,2),xlims=(1.5,T÷2),xticks=2:2:T)

save=false
if save
    path = "plots/$hdf5file/"
    ispath(path) || mkpath(path)
    savefig(plt,joinpath(path,"effctive_masses.pdf"))
end
display(plt)