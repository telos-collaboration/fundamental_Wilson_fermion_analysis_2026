function plot_effective_masses!(plt, meff, Δmeff, h, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p, ncfg, p_label; t1_max=T÷2, t2_max=T÷2,all_non_interacting=false,kws...)
    plot!(plt,ylabel=L"effective mass $[a^{-1}]$",xlabel=L"t",title=L"${%$T} \times {%$L}^3: am^f_0={%$m0}, J^P = 1^-$, ops$ = \pi(\mathbf p)\pi(\mathbf 0), \rho(\mathbf p), \mathbf p = %$(p_label), n_{src}=%$(h), n_{cfg}=%$ncfg, t_0 = %$(t0)$")
    Nev = size(meff)[1]
    for i in 1:Nev
        scatter!(plt,meff[Nev+1-i,1:t1_max],yerr=Δmeff[Nev+1-i,1:t1_max],label="eigenvalue #$i";kws...)
    end
    plot!(plt,ylims=(0.0,π/2),xlims=(1.5,T÷2+0.5),xticks=2:2:T)
    px,py,pz = [parse(Int,c) for c in filter(isdigit,p)]
    label2π  = L"n.i. $E[\pi(\mathbf p)\pi(\mathbf 0)]$" 
    label1ρ  = L"n.i. $E[\rho(\mathbf p)]$" 
    add_mass_band!(plt,non_interacting_energy_2P_lattice(mπ,Δmπ,px,py,pz,L)...;color=:black,label=label2π)
    add_mass_band!(plt,non_interacting_energy_1P_lattice(mρ,Δmρ,px,py,pz,L)... ;color=:black,label=label1ρ)
    if all_non_interacting
        label1π  = L"n.i. $E[\pi(\mathbf p)]$" 
        label1π0 = L"n.i. $E[\pi(\mathbf 0)]$" 
        label1ρ0 = L"n.i. $E[\rho(\mathbf 0)]$" 
        add_mass_band!(plt,non_interacting_energy_1P_lattice(mπ,Δmπ,px,py,pz,L)...;color=:black,label=label1π)
        add_mass_band!(plt,non_interacting_energy_1P(mπ,Δmπ,0 ,L)...;color=:black,label=label1π0)
        add_mass_band!(plt,non_interacting_energy_1P(mρ,Δmρ,0 ,L)...;color=:black,label=label1ρ0)
    end
end
function plot_meff_from_gevp(h5dset,ens,p,t0,deriv,inf_vol;average_equivalent_momenta=true)
    T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
    m0   = only(read(h5dset,joinpath(ens,"quarkmasses")))
    β    = read(h5dset,joinpath(ens,"beta"))
    ncfg = read(h5dset,joinpath(ens,"Nconf"))

    Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits=typemax(Int),average_equivalent_momenta)
    meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv)
    
    ind = findfirst(i -> [β,m0] == inf_vol[i,1:2],1:first(size(inf_vol)))
    mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]
    
    plt = plot(legend=:outerright)
    plot_effective_masses!(plt, meff, Δmeff, sources, T, L, m0, t0, mπ, Δmπ, mρ, Δmρ, p, ncfg, momenta; t1_max=T÷2,t2_max=T÷2)
    three_by_three = haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
    if three_by_three
        Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits=typemax(Int),average_equivalent_momenta)
        Corr3x3[1:2,1:2,:,:] = Corr
        meff_3x3, Δmeff_3x3 = ScatteringI1.effective_masses(Corr3x3;t0,deriv)
        plot_effective_masses!(plt, meff_3x3, Δmeff_3x3, sources, T, L, m0, t0, NaN, Δmπ, NaN, Δmρ, p, ncfg, momenta; t1_max=T÷2,t2_max=T÷2,markersize=3,markershape=:rect)
    end
    return plt
end
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, fitparam; t0, deriv, average_equivalent_momenta=true)
    h5dset  = h5open(corr_file)
    res     = h5open(fitresults)

    plotname = "effective_masses_t0$(t0)_deriv_$deriv.pdf"
    inf_vol  = readdlm(infvolfile,',',skipstart=1)
    fittable = readdlm(fitparam,',',skipstart=1)

    texpath  = joinpath(plotpath,"effective_masses_tex")
    ispath(texpath)  || mkpath(texpath)
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))

    @showprogress desc="Plot effective masses" for ens in keys(h5dset)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        for p in p_external
            
            #joinpath(ens,p) ∉ fittable && continue
            p == "p(0,0,0)" && continue

            plt = plot_meff_from_gevp(h5dset,ens,p,t0,deriv,inf_vol;average_equivalent_momenta)
            
            if haskey(res,joinpath(ens,p))
                r = res[joinpath(ens,p)]
                E0, ΔE0 = read(r,"E0")[1], read(r,"Delta_E0")[1] 
                E1, ΔE1 = read(r,"E1")[1], read(r,"Delta_E1")[1]
                add_mass_band!(plt,E0, ΔE0;label="fit #1")
                add_mass_band!(plt,E1, ΔE1;label="fit #2")
            end
            
            isinteractive() && display(plt)
            savefig(plt,"temp.pdf")
            savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p.tex") )
            append_pdf!(joinpath(plotpath,plotname), "temp.pdf", cleanup=true)
        end
    end
end