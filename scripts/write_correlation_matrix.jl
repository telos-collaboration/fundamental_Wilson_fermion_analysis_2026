function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function write_correlation_matrix(file_in,file_out;combined=true)
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
                
            end
        end
    end
end
