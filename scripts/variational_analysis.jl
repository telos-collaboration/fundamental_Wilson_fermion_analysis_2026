using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
using LatticeUtils
using DelimitedFiles    
pgfplotsx(frame=:box,markersize=5,labelfontsize=16,tickfontsize=14,legendfontsize=14,legend=:bottomleft,markeralpha=0.7)

function non_interacting_energy(mπ,Δmπ,p2,L)
    E1  = sqrt(mπ^2 + p2*(2*pi/L)^2)
    ΔE1 = Δmπ*mπ/E1
    E   = E1 + mπ
    ΔE  = sqrt(Δmπ^2 + ΔE1^2)
    return E, ΔE
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
    meff, Δmeff = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))

    return meff, Δmeff, h
end
function plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, p; t1_max=T÷2, t2_max=T÷2)
    plot!(plt,ylabel="effective mass",xlabel=L"t",title=L"${%$T} \times {%$L}^3: am^f_0={%$m0}, J^P = 1^-$, ops$ = \pi(\mathbf p)\pi(\mathbf 0), \rho(\mathbf p), %$(p), n_{src}=%$h, t_0 = %$(t0)$")
    scatter!(plt,meff[2,1:t1_max],yerr=Δmeff[2,1:t1_max],label="eigenvalue #1")
    scatter!(plt,meff[1,1:t2_max],yerr=Δmeff[1,1:t2_max],label="eigenvalue #2")
    plot!(plt,ylims=(0.2,1.5),xlims=(1.5,T÷2),xticks=2:2:T)
    for p2 in 1:2
        label = p2 == 1 ? L"non-interacting $E_{\pi(\mathbf p)\pi(\mathbf 0)}$ with $|p|=1,\sqrt{2}$" : "" 
        add_mass_band!(plt,non_interacting_energy(mπ,10Δmπ,p2,L)...;label)
    end
end

hdf5file = "data/isospin1_corr.hdf5"
hdf5file = "data/isospin1_corr_allruns.hdf5"
h5dset = h5open(hdf5file)
ens = "Lt24Ls14beta6.9m1-0.92m2-0.92/"
ens = "Lt24Ls14beta6.9m1-0.92m2-0.92/run1"
t0 = 1

p    = "p(0,0,1)"
T, L = h5dset[joinpath("$ens","lattice")][1:2]
m0   = h5dset[joinpath("$ens","quarkmasses")][1]
β    = h5dset[joinpath("$ens","beta")][]
Corr = h5dset[joinpath(ens,p,"correlation_matrix")][]

inf_vol = readdlm("input/infinite_volume.csv",',')
mπ, Δmπ = NaN, NaN

CorrAVG = dropdims(mean(Corr,dims=4),dims=4)
meff, Δmeff, h = effective_masses(Corr;maxhits=typemax(Int),t0)

plt = plot(legend=:outerright)
plot_diag = true
if plot_diag
    for i in 1:2
        m, Δm =  implicit_meff(CorrAVG[i,i,:,:]')
        scatter!(plt,m[1:T÷2], yerr = Δm[1:T÷2], label="diagonal correlator #$i")
    end
end
plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, p; t1_max=T÷2,t2_max=T÷2)
display(plt)