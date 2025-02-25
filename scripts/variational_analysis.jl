using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
using LatticeUtils
using DelimitedFiles    
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

function non_interacting_energy_2π(mπ,Δmπ,p2,L)
    E1, ΔE1 = non_interacting_energy_1π(mπ,Δmπ,p2,L)
    E   = E1 + mπ
    ΔE  = sqrt(Δmπ^2 + ΔE1^2)
    return E, ΔE
end
function non_interacting_energy_1π(mπ,Δmπ,p2,L)
    E1  = sqrt(mπ^2 + p2*(2*pi/L)^2)
    ΔE1 = Δmπ*mπ/E1
    return E1, ΔE1
end
function effective_masses(Corr;t0,maxhits=typemax(Int))

    nhits = size(Corr)[4]
    h     = min(nhits,maxhits)
    Corr  = dropdims(mean(Corr[:,:,:,1:h,:],dims=4),dims=4)

    sign   = +1
    Corr   = correlator_folding(Corr;t_dim=4,sign)
    Corr   = correlator_derivative(Corr;t_dim=4)
    sign   = -1

    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    meff1, Δmeff1 = LatticeUtils.implicit_meff_jackknife(real.(eigvals_resamples);sign)
    meff2, Δmeff2 = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))
    
    return meff1, Δmeff1, meff2, Δmeff2, h
end
function plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p; t1_max=T÷2, t2_max=T÷2,all_non_interacting=false)
    plot!(plt,ylabel=L"effective mass $[a^{-1}]$",xlabel=L"t",title=L"${%$T} \times {%$L}^3: am^f_0={%$m0}, J^P = 1^-$, ops$ = \pi(\mathbf p)\pi(\mathbf 0), \rho(\mathbf p), %$(p), n_{src}=%$h, t_0 = %$(t0)$")
    scatter!(plt,meff[2,1:t1_max],yerr=Δmeff[2,1:t1_max],label="eigenvalue #1")
    scatter!(plt,meff[1,1:t2_max],yerr=Δmeff[1,1:t2_max],label="eigenvalue #2")
    plot!(plt,ylims=(0.2,1.5),xlims=(1.5,T÷2),xticks=2:2:T)
    p2 = sum(x->x^2,[parse(Int,c) for c in filter(isdigit,p)])
    label2π  = L"n.i. $E[\pi(\mathbf p)\pi(\mathbf 0)]$ (err x10)" 
    add_mass_band!(plt,non_interacting_energy_2π(mπ,10Δmπ,p2,L)...;label=label2π)
    if all_non_interacting
        label1π  = L"n.i. $E[\pi(\mathbf p)]$ (err x10)" 
        label1π0 = L"n.i. $E[\pi(\mathbf 0)]$ (err x10)" 
        label1ρ  = L"n.i. $E[\rho(\mathbf p)]$ (err x2)" 
        label1ρ0 = L"n.i. $E[\rho(\mathbf 0)]$ (err x2)" 
        add_mass_band!(plt,non_interacting_energy_1π(mπ,10Δmπ,p2,L)...;color=:black,label=label1π)
        add_mass_band!(plt,non_interacting_energy_1π(mπ,10Δmπ,0 ,L)...;color=:black,label=label1π0)
        add_mass_band!(plt,non_interacting_energy_1π(mρ,2Δmρ,p2,L)...;color=:orange,label=label1ρ)
        add_mass_band!(plt,non_interacting_energy_1π(mρ,2Δmρ,0 ,L)...;color=:orange,label=label1ρ0)
    end
end

hdf5file = "data/isospin1_corr.hdf5"
h5dset = h5open(hdf5file)
ens = "Lt32Ls16beta6.9m1-0.92m2-0.92/"
t0 = 8

p    = "p(0,0,1)"
T, L = h5dset[joinpath("$ens","lattice")][1:2]
m0   = h5dset[joinpath("$ens","quarkmasses")][1]
β    = h5dset[joinpath("$ens","beta")][]
Corr = h5dset[joinpath(ens,p,"correlation_matrix")][]

inf_vol = readdlm("input/infinite_volume.csv",',',skipstart=1)
ind = findfirst(i -> [β,m0] == inf_vol[i,1:2],eachindex(inf_vol))
mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]

meff, Δmeff, meff2, Δmeff2, h = effective_masses(Corr;maxhits=typemax(Int),t0)

plt = plot(legend=:outerright)
plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p; t1_max=T÷2,t2_max=T÷2)
display(plt)