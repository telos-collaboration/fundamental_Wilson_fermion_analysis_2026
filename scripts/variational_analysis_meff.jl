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
function swap_eigval_numbering(old,swap)
    new = copy(old)
    @. new[1,:,1:swap] = old[2,:,1:swap]
    @. new[2,:,1:swap] = old[1,:,1:swap]
    return new
end
function effective_masses(Corr;t0,deriv,maxhits=typemax(Int))

    nhits = size(Corr)[4]
    h     = min(nhits,maxhits)
    Corr  = dropdims(mean(Corr[:,:,:,1:h,:],dims=4),dims=4)

    sign   = +1
    Corr   = correlator_folding(Corr;t_dim=4,sign)
    Corr   = correlator_derivative(Corr;t_dim=4)
    sign   = -1

    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0 -1 )
    meff1, Δmeff1 = LatticeUtils.implicit_meff_jackknife(real.(eigvals_resamples);sign)
    meff2, Δmeff2 = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))
    
    return meff1, Δmeff1, meff2, Δmeff2, h
end
function plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p, ncfg; t1_max=T÷2, t2_max=T÷2,all_non_interacting=false)
    plot!(plt,ylabel=L"effective mass $[a^{-1}]$",xlabel=L"t",title=L"${%$T} \times {%$L}^3: am^f_0={%$m0}, J^P = 1^-$, ops$ = \pi(\mathbf p)\pi(\mathbf 0), \rho(\mathbf p), %$(p), n_{src}=%$h, n_{cfg}=%$ncfg, t_0 = %$(t0)$")
    scatter!(plt,meff[2,1:t1_max],yerr=Δmeff[2,1:t1_max],label="eigenvalue #1")
    scatter!(plt,meff[1,1:t2_max],yerr=Δmeff[1,1:t2_max],label="eigenvalue #2")
    plot!(plt,ylims=(0.0,π/2),xlims=(1.5,T÷2+0.5),xticks=2:2:T)
    p2 = sum(x->x^2,[parse(Int,c) for c in filter(isdigit,p)])
    label2π  = L"n.i. $E[\pi(\mathbf p)\pi(\mathbf 0)]$" 
    label1ρ  = L"n.i. $E[\rho(\mathbf p)]$" 
    add_mass_band!(plt,non_interacting_energy_2π(mπ,Δmπ,p2,L)...;color=:black,label=label2π)
    add_mass_band!(plt,non_interacting_energy_1π(mρ,Δmρ,p2,L)... ;color=:black,label=label1ρ)
    if all_non_interacting
        label1π  = L"n.i. $E[\pi(\mathbf p)]$" 
        label1π0 = L"n.i. $E[\pi(\mathbf 0)]$" 
        label1ρ0 = L"n.i. $E[\rho(\mathbf 0)]$" 
        add_mass_band!(plt,non_interacting_energy_1π(mπ,Δmπ,p2,L)...;color=:black,label=label1π)
        add_mass_band!(plt,non_interacting_energy_1π(mπ,Δmπ,0 ,L)...;color=:black,label=label1π0)
        add_mass_band!(plt,non_interacting_energy_1π(mρ,Δmρ,0 ,L)...;color=:black,label=label1ρ0)
    end
end
function meff_from_gevp(h5dset,ens,p,t0,deriv,inf_vol)
    T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
    m0   = only(read(h5dset,joinpath(ens,"quarkmasses")))
    β    = read(h5dset,joinpath(ens,"beta"))
    ncfg = read(h5dset,joinpath(ens,"Nconf"))
    Corr = read(h5dset,joinpath(ens,p,"correlation_matrix"))

    ind = findfirst(i -> [β,m0] == inf_vol[i,1:2],1:first(size(inf_vol)))
    mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]
    meff, Δmeff, meff2, Δmeff2, h = effective_masses(Corr;maxhits=typemax(Int),t0,deriv)

    plt = plot(legend=:outerright)
    plot_effective_masses!(plt, meff2, Δmeff2, h, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p, ncfg; t1_max=T÷2,t2_max=T÷2)
    return plt
end
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, fitparam; t0, deriv)
    h5dset  = h5open(corr_file)
    res     = h5open(fitresults)

    inf_vol  = readdlm(infvolfile,',',skipstart=1)
    fittable = readdlm(fitparam,',',skipstart=1)

    ispath(plotpath) || mkpath(plotpath)

    for ens in keys(h5dset)
        for p in read(h5dset[ens],"p_external")
            
            joinpath(ens,p) ∉ fittable[:,5] && continue
            p == "p(0,0,0)" && continue

            plt = meff_from_gevp(h5dset,ens,p,t0,deriv,inf_vol)
            
            if haskey(res,joinpath(ens,p))
                r = res[joinpath(ens,p)]
                E0, ΔE0 = read(r,"E0")[1], read(r,"Delta_E0")[1] 
                E1, ΔE1 = read(r,"E1")[1], read(r,"Delta_E1")[1]
                add_mass_band!(plt,E0, ΔE0;label="fit #1")
                add_mass_band!(plt,E1, ΔE1;label="fit #2")
            end
            
            display(plt)
            @show "$(ens)_$p.pdf"
            savefig(joinpath(plotpath,"$(ens)_$p.pdf"))
        end
    end
end

corr_file  = "data/isospin1_corr.hdf5"
fitresults = "data/isospin1_fitresults_t0_3_deriv.hdf5"
plotpath   = "plots/effective_masses"
infvolfile = "input/infinite_volume.csv"
fitparam   = "input/pipi_fitintervals.csv"

deriv = true
t0    = 3

plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, fitparam; t0, deriv)