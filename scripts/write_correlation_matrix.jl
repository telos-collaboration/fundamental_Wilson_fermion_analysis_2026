using Pkg; Pkg.activate(".")
using ScatteringI1
using HDF5
using Statistics
using LatticeUtils
using ProgressMeter
include("read_rhopipi_diagrams.jl")

function pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
    L3, L6 = L^3, L^6
    Corr2π = @. (CorrD1 - CorrD2)/L6 + (CorrR1 + CorrR2 - CorrR3 - CorrR4)/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(ComplexF64,(2,2,N,nhits,T))
    corr[1,1,:,:,:] =  @. Corrρ/L3 + 0*im
    corr[1,2,:,:,:] =  @. 0        + im*(CorrT1-CorrT2)/L3
    corr[2,1,:,:,:] =  @. 0        + im*(CorrT2-CorrT1)/L3
    corr[2,2,:,:,:] =  @. Corr2π   + 0*im
    return corr
end
function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"(APE|Wuppertal)") ,keys(file))
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,entries)
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function write_p001_correlation_matrix(file_in,file_out)
    isfile(file_out) && rm(file_out)
    fid       = h5open(file_in) 
    ensembles = keys(fid)

    @showprogress desc="Write correlation matrices" for ens in ensembles
        _copy_lattice_parameters(file_out,file_in;group=ens)
        T, L = fid["$ens/lattice"][1:2]
        p_external = fid["$ens/p_external"][1:2]
        
        if "p(0,0,1)" ∈ p_external
            Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp001(fid,ens;p=1)
            Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
            Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
            
            h5write(file_out,joinpath(ens,"p(0,0,1)","correlation_matrix"),Corr)
            h5write(file_out,joinpath(ens,"p(0,0,1)","correlator_pion"),Corrπ)
        end
        if "p(0,0,0)" ∈ p_external
            Corrπ0, Corrρ0 = correlatorsp000(fid,ens;p=1)
            h5write(file_out,joinpath(ens,"p(0,0,0)","correlator_pion"),Corrπ0)
            h5write(file_out,joinpath(ens,"p(0,0,0)","correlator_rho") ,Corrρ0)
        end
    end
end
file_in  = "data/isospin1.hdf5"
file_out = "data/isospin1_corr_p001_V3.hdf5"
write_p001_correlation_matrix(file_in,file_out)
