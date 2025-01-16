using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
using LatticeUtils
include("read_rhopipi_diagrams.jl")
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

    eigvals, Δeigvals = eigenvalues(Corr;t0)
    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    meff, Δmeff = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))

    return meff, Δmeff, h
end
function plot_effective_masses(meff, Δmeff, h, T, L, m0, t0, t1_max,t2_max, mπ, Δmπ)
    plt = plot(ylabel="effective mass",xlabel=L"t",title=L"${%$T} \times {%$L}^3: am^f_0={%$m0}, J^P = 1^-$, ops$ = \pi(\mathbf p)\pi(\mathbf 0), \rho(\mathbf p), \mathbf p=(0,0,1), n_{src}=%$h, t_0 = %$(t0)$")
    scatter!(plt,meff[2,1:t1_max],yerr=Δmeff[2,1:t1_max],label="eigenvalue #1")
    scatter!(plt,meff[1,1:t2_max],yerr=Δmeff[1,1:t2_max],label="eigenvalue #2")
    plot!(plt,ylims=(0.2,1.5),xlims=(1.5,T÷2),xticks=2:2:T)
    for p2 in 1:2
        label = p2 == 1 ? L"non-interacting $E_{\pi(\mathbf p)\pi(\mathbf 0)}$ with $|p|=1,\sqrt{2}$" : "" 
        add_mass_band!(plt,non_interacting_energy(mπ,10Δmπ,p2,L)...;label)
    end
    return plt
end

ens = "Lt48Ls16beta7.4m1-0.74m2-0.74"   # (more hits)
ens = "Lt48Ls16beta7.4m1-0.75m2-0.75"   # (more hits)
ens = "Lt32Ls16beta7.2m1-0.794m2-0.794" # (more configs?)
mπ,Δmπ = 0.2532, 0.0007
ens = "Lt32Ls24beta6.9m1-0.924m2-0.924"
mπ,Δmπ = 0.33880, 0.00120
ens = "Lt36Ls20beta7.05m1-0.835m2-0.835" # t0=8; deriv
mπ,Δmπ = 0.43800, 0.00100
ens = "Lt24Ls14beta7.05m1-0.85m2-0.85"   # t0=8; deriv; (bad signal,larger L?)
ens = "Lt32Ls16beta7.05m1-0.85m2-0.85"   # t0=8; deriv; (bad signal,larger L?)
mπ,Δmπ = 0.33076, 0.00097
ens = "Lt32Ls16beta7.2m1-0.78m2-0.78"    # t0=8; deriv (below?) 
mπ,Δmπ = 0.36963, 0.00039
ens = "Lt36Ls24beta7.05m1-0.867m2-0.867" # t0=8; deriv
mπ,Δmπ =  0.14810, 0.00090
ens = "Lt36Ls16beta7.2m1-0.76m2-0.76"    # t0=8; deriv (bad signal,larger L?)
mπ,Δmπ = 0.45700, 0.00130
ens = "Lt32Ls24beta6.9m1-0.92m2-0.92"    # t0=8; deriv
ens = "Lt32Ls16beta6.9m1-0.92m2-0.92"    # t0=8; deriv (bad signal,larger L?)
ens = "Lt24Ls14beta6.9m1-0.92m2-0.92"    # t0=7; deriv (ok signal, better meff)
mπ,Δmπ = 0.38649, 0.00051

hdf5file = "data/isospin1_corr_p001.hdf5"
h5dset = h5open(hdf5file)
t1_max  = 18 
t2_max  = 15
maxhits = 64
t0      = 8

p    = "p(0,0,1)"
T, L = h5dset["$ens/lattice"][1:2]
m0   = -parse(Float64,last(split(ens,'-')))
Corr = h5dset[joinpath(ens,p,"correlation_matrix")][]

meff, Δmeff, h = effective_masses(Corr;maxhits,t0)
plt = plot_effective_masses(meff, Δmeff, h, T, L, m0, t0, t1_max,t2_max, mπ, Δmπ)
plot!(plt,legend=:outerright)
display(plt)
savefig(plt,"gevp_$(ens).pdf")

#=
anim_t0 = Animation()
for t0 in 1:12
    meff, Δmeff, h = effective_masses(CorrD1,CorrR1,CorrD2,CorrR3,Corrρ,CorrT1,CorrT2,L;maxhits,t0)
    plt = plot_effective_masses(meff, Δmeff, h, T, L, m0, t0, t1_max,t2_max, mπ, Δmπ)
    frame(anim_t0,plt)
end
webm(anim_t0,"gevp_$(ens)_t0.webm",fps=1)

anim_h  = Animation()
for maxhits in [2^i for i in 1:7]
    meff, Δmeff, h = effective_masses(CorrD1,CorrR1,CorrD2,CorrR3,Corrρ,CorrT1,CorrT2,L;maxhits,t0)
    plt = plot_effective_masses(meff, Δmeff, h, T, L, m0, t0, t1_max,t2_max, mπ, Δmπ)
    frame(anim_h,plt)
end
webm(anim_h,"gevp_$(ens)_hits.webm",fps=1)
=#