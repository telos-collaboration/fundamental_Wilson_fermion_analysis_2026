function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
pseudolog10(x) = sign(x) * log10(abs(x) + 1)
function write_correlation_matrix(file_in,file_out;combined=true,plotting=true,plotpath="",only_ens=nothing)
    isfile(file_out) && rm(file_out)
    fid = h5open(file_in)
    
    if combined
        ensembles = keys(fid)
    else
        ensembles = AbstractString[]
        for e in keys(fid)
            append!(ensembles,joinpath.(e,keys(fid[e])))
        end
    end

    if plotting
        plotname = "diagrams_t0$(t0)_deriv_$deriv.pdf"
        plotname3x3 = "diagrams_t0$(t0)_deriv_$(deriv)_3x3.pdf"
        texpath  = joinpath(plotpath,"diagrams_tex")
        ispath(texpath) || mkpath(texpath)
        ispath(plotpath) || mkpath(plotpath)
        isfile(joinpath(plotpath,plotname)) && rm(joinpath(plotpath,plotname))
        isfile(joinpath(plotpath,plotname3x3)) && rm(joinpath(plotpath,plotname3x3))
    end

    # Restrict ourselves to only the specified ensembles in the optional keyword argument
    ensembles = isnothing(only_ens) ? ensembles : intersect(ensembles,only_ens)
    
    @showprogress desc="Construct correlation matrices" for ens in ensembles

        _copy_lattice_parameters(file_out,file_in;group=ens)
        T, L = fid["$ens/lattice"][1:2]
        p_external    = read(fid,"$ens/p_external")

        for p0 in p_external 
            if p0 == "p(0,0,0)"
                Corrπ0, Corrρ0 = correlatorsp000(fid,ens)
                h5write(file_out,joinpath(ens,p0,"correlator_pion"),Corrπ0)
                h5write(file_out,joinpath(ens,p0,"correlator_rho") ,Corrρ0)
                h5write(file_out,joinpath(ens,p0,"Nsrc") ,size(Corrπ0)[2])
            else 
                pv = ScatteringI1._parse_momentum(p0)
                # Check if the correlators for the 3x3 correlation matrix have been measured by 
                # checking if the parsed data contains the first relevant dataset
                three_by_three = haskey(fid,"$ens/p($(pv[1]),$(pv[2]),$(pv[3]))/rho_g0g1_g1/p_diag($(pv[1]),$(pv[2]),$(pv[3]))/C_re")
                Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlators_xyz(fid,ens;p=pv)
                
                Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
                Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
                h5write(file_out,joinpath(ens,p0,"correlation_matrix"),Corr)
                h5write(file_out,joinpath(ens,p0,"correlator_pion"),Corrπ)
                h5write(file_out,joinpath(ens,p0,"Nsrc"),size(Corr2π)[2])

                if three_by_three
                    Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2 = correlators_xyz_3x3(fid,ens;p=pv)
                    Corr_3x3_ext = ScatteringI1.pipi_rho_matrix_3x3_extension(Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2,L)
                    h5write(file_out,joinpath(ens,p0,"Nsrc_3x3"),size(Corr_3x3_ext)[2])
                    h5write(file_out,joinpath(ens,p0,"correlation_matrix_3x3_ext"),Corr_3x3_ext)
                end

                if plotting
                    T, L = read(fid,joinpath(ens,"lattice"))[1:2]
                    m0   = only(read(fid,joinpath(ens,"quarkmasses")))
                    
                    corrs  = [Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2] 
                    labels = ["Corrπ", "Corrρ", "CorrT1", "CorrT2", "CorrR1", "CorrR2", "CorrR3", "CorrR4", "CorrD1", "CorrD2"]
                    markers = [:circle, :diamond, :dtriangle, :pentagon, :rect, :rtriangle, :utriangle, :star6, :xcross, :vline]
                    
                    mi = 1
                    plt  = plot(ylabel=L"${\rm pseudolog}_{10} C_(t)$",xlabel=L"t",legend=:outerright)
                    for (C_tmp,l) in zip(corrs,labels)
                        C_tmp .= pseudolog10.(C_tmp)
                        ncfg = first(size(C_tmp))
                        C  = dropdims(mean(C_tmp,dims=(1,2)),dims=(1,2))
                        ΔC = dropdims(std(mean(C_tmp,dims=(2)),dims=1),dims=(1,2))/sqrt(ncfg)
                        
                        title = L"${%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf p = %$(p0), n_{cfg}=%$ncfg, t_0 = %$(t0)$"   
                        scatter!(plt,1:T,C,yerr=ΔC,label=l,marker=markers[mi],title=title)
                        mi += 1
                    end

                    savefig(plt,"temp.pdf")
                    savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$p0.tex") )
                    append_pdf!(joinpath(plotpath,plotname),"temp.pdf",cleanup=true)
                    isinteractive() && display(plt)

                    if three_by_three
                        corrs  = [Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2]
                        labels = ["Corr_γ0γi_γi","Corr_γi_γ0γi", "Corr_γ0γi_γ0γi", "Corrγ0γiT1", "Corrγ0γiT2"]

                        mi = 1
                        plt  = plot(ylabel=L"${\rm pseudolog}_{10} C_(t)$",xlabel=L"t",legend=:outerright)
                        for (C_tmp,l) in zip(corrs,labels)
                            C_tmp .= pseudolog10.(C_tmp)
                            ncfg = first(size(C_tmp))
                            C  = dropdims(mean(C_tmp,dims=(1,2)),dims=(1,2))
                            ΔC = dropdims(std(mean(C_tmp,dims=(2)),dims=1),dims=(1,2))/sqrt(ncfg)
                            
                            title = L"${%$T} \times {%$L}^3: am^f_0={%$m0}, \mathbf p = %$(p0), n_{cfg}=%$ncfg, t_0 = %$(t0)$"   
                            scatter!(plt,1:T,C,yerr=ΔC,label=l,marker=markers[mi],title=title)
                            mi += 1
                        end
    
                        savefig(plt,"temp.pdf")
                        savefig(plot!(plt,tex_output_standalone = true), joinpath(texpath,"$(ens)_$(p0)_3by3.tex") )
                        append_pdf!(joinpath(plotpath,plotname3x3),"temp.pdf",cleanup=true)
                        isinteractive() && display(plt)    
                    end


                end
            end
        end
    end
end
