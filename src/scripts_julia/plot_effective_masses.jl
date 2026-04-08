using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using HDF5: h5open, h5write
using ScatteringI1
using LaTeXStrings: @L_str
using Plots: gr, plot, plot!, scatter!, savefig, backend_name, ylims
using PDFmerger: append_pdf!
using DelimitedFiles: readdlm
gr(fontfamily="Computer Modern",frame=:box, legend=:topright, size = (400,300), markeralpha=0.7,titlefontsize=9)

function plot_effective_mass!(plt, meff, Δmeff ;kws...)
    T = length(meff)
    tmax1 = findfirst(t->abs(Δmeff[t]/meff[t]) > 0.4, 1:T÷2)
    tmax1 = isnothing(tmax1) ? T÷2 : tmax1 - 1
    tmax2 = findfirst(t->abs(Δmeff[t]/meff[t]) > 0.4, T:-1:T÷2)
    tmax2 = isnothing(tmax2) ? T÷2 : tmax2 - 1
    t = vcat(1:tmax1,T:-1:T-tmax2)
    scatter!(plt,t,meff[t],yerr=Δmeff[t];kws...)
end
function plot_effective_masses!(plt, meff, Δmeff, sources ;kws...)
    Nev,T = size(meff)
    Nev_max = 2
    for n in Nev_max:-1:1
        tmax1 = findfirst(t->abs(Δmeff[n,t]/meff[n,t]) > 0.5, 1:T÷2)
        tmax1 = isnothing(tmax1) ? T÷2 : tmax1 - 1
        tmax2 = findfirst(t->abs(Δmeff[n,t]/meff[n,t]) > 0.5, T:-1:T÷2)
        tmax2 = isnothing(tmax2) ? T÷2 : tmax2 - 1
        t = vcat(1:tmax1,T:-1:T-tmax2)
        scatter!(plt,t,meff[n,t],yerr=Δmeff[n,t],label=L"\textrm{eigenvalue }~%$n";kws...)
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
function plot_effective_masses(corr_file, fitresults, infvolfile, plotpath, metadata, basename, plot_mesons; plot2x2 = false)
    h5dset  = h5open(corr_file)
    if isfile(fitresults)
        res = h5open(fitresults)
    end

    plotname = "$(basename)_(g)evp.pdf"
    plotname_mesons = "$(basename)_mesons.pdf"
    plotname_mesons_p0 = "$(basename)_mesons_p0.pdf"
    
    inf_vol  = readdlm(infvolfile,',',skipstart=1)
    data = readdlm(metadata,',',skipstart=1)
    ispath(plotpath) || mkpath(plotpath)
    isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
    isfile(joinpath(plotpath,plotname_mesons)) && rm(joinpath(plotpath,plotname_mesons))
    isfile(joinpath(plotpath,plotname_mesons_p0)) && rm(joinpath(plotpath,plotname_mesons_p0))
    
    if plot_mesons
        for ens in unique(data[:,1])
            p = "p(0,0,0)"
            T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
            ncfg = read(h5dset,joinpath(ens,"Nconf"))
            m0 = only(read(h5dset,joinpath(ens,"quarkmasses")))
            β  = read(h5dset,joinpath(ens,"beta"))
            p_fmt = replace(p,"p"=>"")
            title  = L"{%$T} \times {%$L}^3, \beta=%$β, am_0={%$m0}, \vec{d} = %$p_fmt)"
            plt_mesons = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"am_\textrm{eff}")
            if haskey(h5dset[ens][p],"T1")
                meff = read(h5dset[ens][p]["T1"],"meff")
                Δmeff = read(h5dset[ens][p]["T1"],"Delta_meff")
                plot_effective_mass!(plt_mesons, meff, Δmeff, label=L"\rho (T_1)")
            end
            if isfile(fitresults) haskey(res,joinpath(ens,p,"T1"))
                r = res[joinpath(ens,p,"T1")]
                E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1]
                tmin = read(r,"tmin") + 1
                tmax = read(r,"tmax") + 1
                add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
            end
            if isfile(fitresults) && haskey(res,joinpath(ens,p))
                if haskey(res,joinpath(ens,p,"pi"))
                    r = res[joinpath(ens,p,"pi")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    tmin = read(r,"tmin") + 1
                    tmax = read(r,"tmax") + 1
                    add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
                end
            end
            if haskey(h5dset[ens][p],"meff_pi")
                meff = read(h5dset[ens][p],"meff_pi")
                Δmeff = read(h5dset[ens][p],"Delta_meff_pi")
                plot_effective_mass!(plt_mesons, meff, Δmeff, label=L"\pi")
            end
            savefig(plt_mesons,"temp_p0.pdf")
            append_pdf!(joinpath(plotpath,plotname_mesons_p0), "temp_p0.pdf", cleanup=true)
        end
    end

    @showprogress desc="Plot effective masses:" for row in eachrow(data)
        
        ens, p, irrep, id = row[1], row[2], row[3], row[14]
        t0, deriv, gevp, use3x3 = Int(row[8]), Bool(row[9]), Bool(row[10]), Bool(row[12])

        T, L = read(h5dset,joinpath(ens,"lattice"))[1:2]
        ncfg = read(h5dset,joinpath(ens,"Nconf"))
        m0 = only(read(h5dset,joinpath(ens,"quarkmasses")))
        β  = read(h5dset,joinpath(ens,"beta"))
        p_fmt = replace(p,"p"=>"")
        title  = L"{%$T} \times {%$L}^3, \beta=%$β, am_0={%$m0}, \vec{d} = %$(p_fmt)"

        if plot_mesons
            plt_mesons = plot(;title,legend=:bottomleft,xlabel=L"t",ylabel=L"am_\textrm{eff }")
            if isfile(fitresults) && haskey(res,joinpath(ens,p))
                if haskey(res,joinpath(ens,p,"pi"))
                    r = res[joinpath(ens,p,"pi")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    tmin = read(r,"tmin") + 1
                    tmax = read(r,"tmax") + 1
                    add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
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
            if isfile(fitresults) && haskey(res,joinpath(ens,p))
                r = res[joinpath(ens,p,"A1")]
                if haskey(res,joinpath(ens,p,"B1"))
                    r = res[joinpath(ens,p,"B1")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1]
                    tmin = read(r,"tmin") + 1
                    tmax = read(r,"tmax") + 1
                    add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
                end
                if haskey(res,joinpath(ens,p,"E"))
                    r = res[joinpath(ens,p,"E")]
                    E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                    tmin = read(r,"tmin") + 1
                    tmax = read(r,"tmax") + 1
                    add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
                end
            end
            savefig(plt_mesons ,"temp_mesons.pdf")
            append_pdf!(joinpath(plotpath,plotname_mesons), "temp_mesons.pdf", cleanup=true)
        end

        # write title and axis labels
        t0 = read(h5dset,joinpath(ens,p,irrep,id,"t0"))
        gevp = read(h5dset,joinpath(ens,p,irrep,id,"gevp"))
        deriv = read(h5dset,joinpath(ens,p,irrep,id,"deriv"))
        symmetrise = read(h5dset,joinpath(ens,p,irrep,id,"symmetrise"))
        if gevp
            title  = L"{%$T} \times {%$L}^3, \beta=%$β, am_0={%$m0}, \vec{d} = %$(p_fmt), \mathrm{GEVP}, t_0 = %$(t0)"
        else
            title  = L"{%$T} \times {%$L}^3, \beta=%$β, am_0={%$m0}, \vec{d} = %$(p_fmt), \mathrm{EVP}"
        end
        plt = plot(;title,xlabel=L"t",ylabel=L"am_\textrm{eff}")
        # get metadate for specific momentum
        data = readdlm(metadata,',',skipstart=1)
        metadata_ind = findfirst(i -> isequal(joinpath(ens,p),joinpath(data[i,1:2]...)),1:first(size(data)))
        use3x3 = data[metadata_ind,12]
        has3x3 = haskey(h5dset[ens][p]["$irrep/$id"],"meff_3x3")
        if has3x3 && use3x3
            meff_3x3 = read(h5dset[ens][p]["$irrep/$id"],"meff_3x3")
            Δmeff_3x3 = read(h5dset[ens][p]["$irrep/$id"],"Delta_meff_3x3")
            sources_3x3 = read(h5dset[ens][p]["$irrep/$id"],"sources_3x3")
            plot_effective_masses!(plt, meff_3x3, Δmeff_3x3, sources_3x3)
        end
        if !(has3x3 && use3x3) || plot2x2
            meff = read(h5dset[ens][p]["$irrep/$id"],"meff")
            Δmeff = read(h5dset[ens][p]["$irrep/$id"],"Delta_meff")
            sources = read(h5dset[ens][p]["$irrep/$id"],"sources")
            plot_effective_masses!(plt, meff, Δmeff, sources; markershape=:rect)
        end
        #plot_non_interacting_levels!(plt,h5dset,ens,p,inf_vol)
        # use default x- and y-limits (override if we have fits)
        y_min = +Inf
        y_max = -Inf
        if isfile(fitresults) && haskey(res,joinpath(ens,id,p))
            r = res[joinpath(ens,id,p,"A1")]
            E0, ΔE0 = read(r,"E")[1], read(r,"Delta_E")[1] 
            E1, ΔE1 = read(r,"E")[2], read(r,"Delta_E")[2]
            tmin1 = read(r,"tmin1") + 1
            tmax1 = read(r,"tmax1") + 1
            tmin2 = read(r,"tmin2") + 1
            tmax2 = read(r,"tmax2") + 1
            add_fit_range!(plt, tmin1, tmax1, E0, ΔE0;label="fit #1")
            add_fit_range!(plt, tmin2, tmax2, E1, ΔE1;label="fit #2")

            # update plot limits
            y_min = min(0.75*(E0-ΔE0),y_min)
            y_max = max(1.25*(E1+ΔE1),y_max)
            if has3x3 && use3x3
                meff_3x3 = read(h5dset[ens][p]["$irrep/$id"],"meff_3x3")
                Δmeff_3x3 = read(h5dset[ens][p]["$irrep/$id"],"Delta_meff_3x3")
                y_min = min(y_min,minimum(meff_3x3[end,  tmin2:tmax2]))
                y_max = min(y_max,maximum(meff_3x3[end-1,tmin2:tmax2]))
            end
            if !(has3x3 && use3x3) || plot2x2
                meff = read(h5dset[ens][p]["$irrep/$id"],"meff")
                Δmeff = read(h5dset[ens][p]["$irrep/$id"],"Delta_meff")
                y_min = min(y_min,minimum(meff[end,  tmin2:tmax2]))
                y_max = min(y_max,maximum(meff[end-1,tmin2:tmax2]))
            end

            if haskey(res,joinpath(ens,p,"B1"))
                r = res[joinpath(ens,p,"B1")]
                E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1]
                tmin = read(r,"tmin") + 1
                tmax = read(r,"tmax") + 1
                add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
            end
            if haskey(res,joinpath(ens,p,"E"))
                r = res[joinpath(ens,p,"E")]
                E, ΔE = read(r,"E")[1], read(r,"Delta_E")[1] 
                tmin = read(r,"tmin") + 1
                tmax = read(r,"tmax") + 1
                add_fit_range!(plt_mesons, tmin, tmax, E, ΔE;label="")
            end
        end
        plot!(plt,xlims=(1.5,T÷2 - 0.5),xticks=2:2:T÷2)
        plot!(plt,ylims=(0,π/2))

        tmpfile = tempname()*".pdf"
        savefig(plt,tmpfile)
        append_pdf!(joinpath(plotpath,plotname), tmpfile, cleanup=true)
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
        "--plotbasename"
        help = "Naming scheme for the PDF files geneate by this script"
        default = "effective_masses"
        "--plot_mesons"
        help = "Also plot effective masses for single meson operators"
        arg_type = Bool
        default = true
        "--metadata"
        help = "CSV file containing the parameters for the variational analysis"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plot_effective_masses(args["h5file_eig"], args["h5file_fit"], args["infinite_volume"], args["plotpath"], args["metadata"], args["plotbasename"],args["plot_mesons"])
end
main()