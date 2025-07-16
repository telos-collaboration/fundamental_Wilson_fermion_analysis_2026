using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using HDF5: h5open, h5write
using ScatteringI1
using LaTeXStrings: @L_str
using Plots: gr, plot, plot!, scatter!, savefig, backend_name
using PDFmerger: append_pdf!
using DelimitedFiles: readdlm
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function plot_effective_mass!(plt, meff, Δmeff ;kws...)
    T = length(meff)
    tmax = findfirst(t->abs(Δmeff[t]/meff[t]) > 0.5, 1:T÷2)
    tmax = isnothing(tmax) ? T÷2 : tmax - 1
    t = 1:tmax
    scatter!(plt,t,meff[t],yerr=Δmeff[t];kws...)
end
function plot_effective_masses!(plt, meff, Δmeff, sources ;kws...)
    Nev,T = size(meff)
    tmax  = 
    for i in 1:Nev
        n = Nev+1-i
        tmax = findfirst(t->abs(Δmeff[n,t]/meff[n,t]) > 0.5, 1:T÷2)
        tmax = isnothing(tmax) ? T÷2 : tmax - 1
        t = 1:tmax
        scatter!(plt,t,meff[n,t],yerr=Δmeff[n,t],label=L"\textrm{eigenvalue }~%$i~~(n_{src}=%$(sources))";kws...)
    end
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
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath; plot2x2 = false)
    h5dset  = h5open(corr_file)
    if isfile(fitresults)
        res = h5open(fitresults)
    end

    plotname = "effective_masses_(g)evp.pdf"
    plotname_mesons = "effective_masses_mesons.pdf"
    inf_vol  = readdlm(infvolfile,',',skipstart=1)
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
    isfile(joinpath(plotpath,plotname_mesons)) && rm(joinpath(plotpath,plotname_mesons))

    @showprogress desc="Plot effective masses" for ens in keys(h5dset)
        p_external = read(h5dset,"$ens/p_external")
        for p in p_external
            
            p == "p(0,0,0)" && continue

            # write title and axis labels
            T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
            ncfg = read(h5dset,joinpath(ens,"Nconf"))
            m0 = only(read(h5dset,joinpath(ens,"quarkmasses")))
            β  = read(h5dset,joinpath(ens,"beta"))
            t0 = read(h5dset,joinpath(ens,p,"A1","t0"))
            gevp = read(h5dset,joinpath(ens,p,"A1","gevp"))
            deriv = read(h5dset,joinpath(ens,p,"A1","deriv"))
            symmetrise = read(h5dset,joinpath(ens,p,"A1","symmetrise"))

            if gevp
                title  = L"{%$T} \times {%$L}^3: \beta=%$β, am^f_0={%$m0}, \mathbf{p} = %$(p), n_{cfg}=%$ncfg, gevp, t_0 = %$(t0)"
            else
                title  = L"{%$T} \times {%$L}^3: \beta=%$β, am^f_0={%$m0}, \mathbf{p} = %$(p), n_{cfg}=%$ncfg, evp"
            end
            plt = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"\textrm{effective mass } [a^{-1}]")
            plt_mesons = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"\textrm{effective mass } [a^{-1}]")

            
            has3x3 = haskey(h5dset[ens][p]["A1"],"meff_3x3")
            if has3x3
                meff_3x3 = read(h5dset[ens][p]["A1"],"meff_3x3")
                Δmeff_3x3 = read(h5dset[ens][p]["A1"],"Delta_meff_3x3")
                sources_3x3 = read(h5dset[ens][p]["A1"],"sources_3x3")
                plot_effective_masses!(plt, meff_3x3, Δmeff_3x3, sources_3x3)
            end
            if !has3x3 || plot2x2
                meff = read(h5dset[ens][p]["A1"],"meff")
                Δmeff = read(h5dset[ens][p]["A1"],"Delta_meff")
                sources = read(h5dset[ens][p]["A1"],"sources")
                plot_effective_masses!(plt, meff, Δmeff, sources; markershape=:rect)
            end
            plot_non_interacting_levels!(plt,h5dset,ens,p,inf_vol)
            
            if isfile(fitresults) && haskey(res,joinpath(ens,p))
                r = res[joinpath(ens,p,"A1")]
                E0, ΔE0 = read(r,"E0")[1], read(r,"Delta_E0")[1] 
                E1, ΔE1 = read(r,"E1")[1], read(r,"Delta_E1")[1]
                add_mass_band!(plt,E0, ΔE0;label="fit #1")
                add_mass_band!(plt,E1, ΔE1;label="fit #2")

                if haskey(res,joinpath(ens,p,"B1"))
                    r = res[joinpath(ens,p,"B1")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    add_mass_band!(plt_mesons,E, ΔE;label="")
                end
                if haskey(res,joinpath(ens,p,"E"))
                    r = res[joinpath(ens,p,"E")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    add_mass_band!(plt_mesons,E, ΔE;label="")
                end
                if haskey(res,joinpath(ens,p,"pi"))
                    r = res[joinpath(ens,p,"pi")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    add_mass_band!(plt_mesons,E, ΔE;label="")
                end

            end
            
            if haskey(h5dset[ens][p],"meff_pi")
                meff = read(h5dset[ens][p],"meff_pi")
                Δmeff = read(h5dset[ens][p],"Delta_meff_pi")
                plot_effective_mass!(plt_mesons, meff, Δmeff, label=L"\pi")
            end

            if haskey(h5dset[ens][p],"E")
                meff = read(h5dset[ens][p]["E"],"meff")
                Δmeff = read(h5dset[ens][p]["E"],"Delta_meff")
                plot_effective_mass!(plt_mesons, meff, Δmeff, label=L"\rho (E)")
            end

            if haskey(h5dset[ens][p],"B1")
                meff = read(h5dset[ens][p]["B1"],"meff")
                Δmeff = read(h5dset[ens][p]["B1"],"Delta_meff")
                plot_effective_mass!(plt_mesons, meff, Δmeff, label=L"\rho (B1)")
            end

            plot!(plt,ylims=(0.0,π/2),xlims=(1.5,T÷2+0.5),xticks=2:2:T)
            savefig(plt,"temp.pdf")
            savefig(plt_mesons ,"temp_mesons.pdf")
            append_pdf!(joinpath(plotpath,plotname), "temp.pdf", cleanup=true)
            append_pdf!(joinpath(plotpath,plotname_mesons), "temp_mesons.pdf", cleanup=true)
        end
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_eig"
        help = "HDF5 file containing the eigenvalues"
        required = true
        "--h5file_fit"
        help = "HDF5 file containing the fit results"
        required = true
        "--infinite_volume"
        help = "CSV file containing the infinite volume results"
        required = true
        "--plotpath"
        help = "HDF5 output file containing the correlation matrices"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plot_effective_masses(args["h5file_eig"], args["h5file_fit"], args["infinite_volume"], args["plotpath"])
end
main()