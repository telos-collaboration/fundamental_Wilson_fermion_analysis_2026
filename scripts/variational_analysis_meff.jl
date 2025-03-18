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
    mπ, Δmπ, mρ, Δmρ = inf_vol[ind,3:6]
    px,py,pz = [parse(Int,c) for c in filter(isdigit,p)]
    label2π  = L"\textrm{n.i.} E[\pi(\mathbf{p})\pi(\mathbf{p})]" 
    label1ρ  = L"\textrm{n.i.} E[\rho(\mathbf{p})]" 
    add_mass_band!(plt,non_interacting_energy_2P_lattice(mπ,Δmπ,px,py,pz,L)...;color=:black,label=label2π)
    add_mass_band!(plt,non_interacting_energy_1P_lattice(mρ,Δmρ,px,py,pz,L)... ;color=:black,label=label1ρ)
end
function construct_3x3_correlation_matrix(h5dset,ens,p;maxhits=typemax(Int))
    Corr2x2, sources2x2, momenta2x2 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits)
    Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits)
    Corr3x3[1:2,1:2,:,:] = Corr2x2
    return Corr3x3, sources3x3, momenta3x3
end
function plot_meff_from_gevp!(plot2x2,h5dset,ens,p,t0,deriv;use3x3=true,half_sources=false)
    three_by_three = haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
    if three_by_three && use3x3
        Corr, sources, momenta = construct_3x3_correlation_matrix(h5dset,ens,p;maxhits=typemax(Int))
        meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv)
        plot_effective_masses!(plot2x2, meff, Δmeff, sources; markershape=:rect)
        if half_sources
            Corr, sources, momenta = construct_3x3_correlation_matrix(h5dset,ens,p;maxhits=maximum(sources)÷2)
            meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv)
            plot_effective_masses!(plot2x2, meff, Δmeff, sources; markershape=:rect)    
        end
    elseif !use3x3
        Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits=typemax(Int))
        meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv)
        plot_effective_masses!(plot2x2, meff, Δmeff, sources)
        if half_sources
            Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits=maximum(sources)÷2)
            meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv)
            plot_effective_masses!(plot2x2, meff, Δmeff, sources)    
        end
    end
end
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, fitparam; t0, deriv, average_equivalent_momenta=true, use3x3=true)
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

            # write title and axis labels
            T, L   = read(h5dset,joinpath(ens,"lattice"))[1:2]
            m0     = only(read(h5dset,joinpath(ens,"quarkmasses")))
            β      = read(h5dset,joinpath(ens,"beta"))
            ncfg   = read(h5dset,joinpath(ens,"Nconf"))
            title  = L"{%$T} \times {%$L}^3: \beta=%$β, am^f_0={%$m0}, \mathbf{p} = %$(p), n_{cfg}=%$ncfg, t_0 = %$(t0)"
            plt1  = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"\textrm{effective mass } [a^{-1}]")
            plt2 = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"\textrm{effective mass } [a^{-1}]")
        
            plot_meff_from_gevp!(plt1,h5dset,ens,p,t0,deriv;use3x3=false)
            plot_meff_from_gevp!(plt1,h5dset,ens,p,t0,deriv;use3x3)
            plot_non_interacting_levels!(plt1,h5dset,ens,p,inf_vol)
            plot_meff_from_gevp!(plt2,h5dset,ens,p,t0,deriv;use3x3=false,half_sources=true)
            
            for plt in [plt1,plt2]
                if haskey(res,joinpath(ens,p))
                    r = res[joinpath(ens,p)]
                    E0, ΔE0 = read(r,"E0")[1], read(r,"Delta_E0")[1] 
                    E1, ΔE1 = read(r,"E1")[1], read(r,"Delta_E1")[1]
                    add_mass_band!(plt,E0, ΔE0;label="fit #1")
                    add_mass_band!(plt,E1, ΔE1;label="fit #2")
                end
                
                isinteractive() && display(plt)
                savefig(plt,"temp.pdf")
                if backend_name() == :pgfplotsx
                    savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p.tex") )
                end
                append_pdf!(joinpath(plotpath,plotname), "temp.pdf", cleanup=true)
            end
        end
    end
end