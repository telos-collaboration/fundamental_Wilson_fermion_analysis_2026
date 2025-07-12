function plot_effective_masses!(plt, meff, Δmeff, sources ;kws...)
    Nev,T = size(meff)
    tmax  = T÷2
    for i in 1:Nev
        scatter!(plt,meff[Nev+1-i,1:tmax],yerr=Δmeff[Nev+1-i,1:tmax],label=L"\textrm{eigenvalue }~%$i~~(n_{src}=%$(sources))";kws...)
    end
    plot!(plt,ylims=(0.0,π/2),xlims=(1.5,T÷2+0.5),xticks=2:2:T)
end
function plot_non_interacting_levels!(plt,h5dset,ens,p,inf_vol)
    T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
    m0   = only(read(h5dset,joinpath(ens,"quarkmasses")))
    β    = read(h5dset,joinpath(ens,"beta"))
    
    ind = findfirst(i -> [β,m0] == inf_vol[i,1:2],1:first(size(inf_vol)))
    if !isnothing(ind)
        mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]
        px,py,pz = [parse(Int,c) for c in filter(isdigit,p)]
        label2π  = L"\textrm{n.i.} E[\pi(\mathbf{p})\pi(\mathbf{0})]" 
        label1ρ  = L"\textrm{n.i.} E[\rho(\mathbf{p})]" 
        add_mass_band!(plt,non_interacting_energy_2P_lattice(mπ,Δmπ,px,py,pz,L)...;color=:black,label=label2π)
        add_mass_band!(plt,non_interacting_energy_1P_lattice(mρ,Δmρ,px,py,pz,L)... ;color=:black,label=label1ρ)
    end
end
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath; use3x3=true)
    h5dset  = h5open(corr_file)
    if isfile(fitresults)
        res = h5open(fitresults)
    end

    plotname = "effective_masses_(g)evp.pdf"
    inf_vol  = readdlm(infvolfile,',',skipstart=1)
    texpath  = joinpath(plotpath,"effective_masses_tex")
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))

    @showprogress desc="Plot effective masses" for ens in keys(h5dset)
        p_external = read(h5dset,"$ens/p_external")
        for p in p_external
            
            p == "p(0,0,0)" && continue

            # write title and axis labels
            T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
            ncfg = read(h5dset,joinpath(ens,"Nconf"))
            m0 = only(read(h5dset,joinpath(ens,"quarkmasses")))
            β  = read(h5dset,joinpath(ens,"beta"))
            t0 = read(h5dset,joinpath(ens,p,"t0"))
            gevp = read(h5dset,joinpath(ens,p,"gevp"))
            deriv = read(h5dset,joinpath(ens,p,"deriv"))
            symmetrise = read(h5dset,joinpath(ens,p,"symmetrise"))

            if gevp
                title  = L"{%$T} \times {%$L}^3: \beta=%$β, am^f_0={%$m0}, \mathbf{p} = %$(p), n_{cfg}=%$ncfg, gevp, t_0 = %$(t0)"
            else
                title  = L"{%$T} \times {%$L}^3: \beta=%$β, am^f_0={%$m0}, \mathbf{p} = %$(p), n_{cfg}=%$ncfg, evp"
            end
            plt1 = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"\textrm{effective mass } [a^{-1}]")

            meff = read(h5dset[ens][p],"meff")
            Δmeff = read(h5dset[ens][p],"Delta_meff")
            sources = read(h5dset[ens][p],"sources")
            plot_effective_masses!(plt1, meff, Δmeff, sources; markershape=:rect)
            plot_non_interacting_levels!(plt1,h5dset,ens,p,inf_vol)
            
            has3x3 = haskey(h5dset[ens][p],"meff_3x3")
            if use3x3 && has3x3
                meff_3x3 = read(h5dset[ens][p],"meff_3x3")
                Δmeff_3x3 = read(h5dset[ens][p],"Delta_meff_3x3")
                sources_3x3 = read(h5dset[ens][p],"sources_3x3")
                plot_effective_masses!(plt1, meff_3x3, Δmeff_3x3, sources_3x3)
            end
              
            for plt in [plt1]
                if isfile(fitresults) && haskey(res,joinpath(ens,p))
                    r = res[joinpath(ens,p,"A1")]
                    E0, ΔE0 = read(r,"E0")[1], read(r,"Delta_E0")[1] 
                    E1, ΔE1 = read(r,"E1")[1], read(r,"Delta_E1")[1]
                    add_mass_band!(plt,E0, ΔE0;label="fit #1")
                    add_mass_band!(plt,E1, ΔE1;label="fit #2")
                end
                
                isinteractive() && display(plt)
                savefig(plt,"temp.pdf")
                if backend_name() == :pgfplotsx
                    ispath(texpath)  || mkpath(texpath)
                    savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p.tex") )
                end
                append_pdf!(joinpath(plotpath,plotname), "temp.pdf", cleanup=true)
            end
        end
    end
end