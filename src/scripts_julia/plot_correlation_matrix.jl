using ProgressMeter: @showprogress
using HDF5: h5open
using ScatteringI1
using LaTeXStrings: @L_str
using Statistics: mean, std
using Plots: gr, plot, scatter!, savefig, backend_name
using PDFmerger: append_pdf!
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

pseudolog10(x,C=1E+4) = sign(x)*log10(abs(C*x)+1)
function plot_correlation_matrices(file_in,plotpath;only_ens=nothing)
    fid = h5open(file_in)
    # TODO: Rewrite using the filter! function
    ensembles = keys(fid)
    ensembles = isnothing(only_ens) ? ensembles : intersect(ensembles,only_ens)

    if plotting
        plotname = "diagrams_v2.pdf"
        plotname3x3 = "diagrams_3x3_v2.pdf"
        texpath  = joinpath(plotpath,"diagrams_tex")
        ispath(plotpath) || mkpath(plotpath)
        isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
        isfile(joinpath(plotpath,plotname3x3)) && rm(joinpath(plotpath,plotname3x3))
    end

    @showprogress desc="Plot correlation matrix elements" for ens in ensembles

        T, L = read(fid,joinpath(ens,"lattice"))[1:2]
        p_external = read(fid,"$ens/p_external")

        for p0 in p_external 
            if p0 == "p(0,0,0)"
                continue
            else 
                pv = ScatteringI1._parse_momentum(p0)
                three_by_three = haskey(fid,"$ens/p($(pv[1]),$(pv[2]),$(pv[3]))/rho_g0g1_g1/p_diag($(pv[1]),$(pv[2]),$(pv[3]))/C_re")

                Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlators_xyz(fid,ens;p=pv)

                T, L  = read(fid,joinpath(ens,"lattice"))[1:2]
                m0    = only(read(fid,joinpath(ens,"quarkmasses")))
                ncfg  = first(size(Corrπ))
                title = L"{%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf{p} = %$(p0), n_{cfg}=%$ncfg, t_0 = %$(t0)"
                ylabel= L"\textrm{pseudolog}_{10} C(t)"
                xlabel= L"t"
                    
                corrs  = [Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2] 
                labels = ["Corrπ", "Corrρ", "CorrT1", "CorrT2", "CorrR1", "CorrR2", "CorrR3", "CorrR4", "CorrD1", "CorrD2"]
                markers = [:circle, :diamond, :dtriangle, :pentagon, :rect, :rtriangle, :utriangle, :star6, :xcross, :vline]
                
                mi = 1
                plt  = plot(legend=:outerright; xlabel, ylabel, title)
                for (C_tmp,l) in zip(corrs,labels)
                    C_tmp .= pseudolog10.(C_tmp)
                    C  = dropdims(mean(C_tmp,dims=(1,2)),dims=(1,2))
                    ΔC = dropdims(std(mean(C_tmp,dims=(2)),dims=1),dims=(1,2))/sqrt(ncfg)
                    scatter!(plt,1:T,C,yerr=ΔC,label=l,marker=markers[mi])
                    mi += 1
                end 
                
                savefig(plt,"temp.pdf")
                if backend_name() == :pgfplotsx
                    savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p0.tex") )
                end
                append_pdf!(joinpath(plotpath,plotname),"temp.pdf",cleanup=true)
                isinteractive() && display(plt)

                if three_by_three
                    Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2 = correlators_xyz_3x3(fid,ens;p=pv)
                    corrs  = [Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2]
                    labels = ["Corr_γ0γi_γi","Corr_γi_γ0γi", "Corr_γ0γi_γ0γi", "Corrγ0γiT1", "Corrγ0γiT2"]
                    
                    mi = 1
                    plt  = plot(legend=:outerright; xlabel, ylabel, title)
                    for (C_tmp,l) in zip(corrs,labels)
                        C_tmp .= pseudolog10.(C_tmp)
                        C  = dropdims(mean(C_tmp,dims=(1,2)),dims=(1,2))
                        ΔC = dropdims(std(mean(C_tmp,dims=(2)),dims=1),dims=(1,2))/sqrt(ncfg)
                        scatter!(plt,1:T,C,yerr=ΔC,label=l,marker=markers[mi])
                        mi += 1
                    end

                    if backend_name() == :pgfplotsx
                        ispath(texpath) || mkpath(texpath)
                        savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$(p0)_3by3.tex") )
                    end

                    savefig(plt,"temp.pdf")
                    append_pdf!(joinpath(plotpath,plotname3x3),"temp.pdf",cleanup=true)
                    isinteractive() && display(plt)    
                end
            end
        end
    end
end