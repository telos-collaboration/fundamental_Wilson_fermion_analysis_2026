using Pkg; Pkg.activate(".")
using ScatteringI1
using HDF5
using Statistics
using LatticeUtils
using ProgressMeter

# TODO: Average over sources here, so that the code deals with differing number of sources
function pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
    L3, L6 = L^3, L^6
    Corr2π = @. (CorrD1 - CorrD2)/L6 + (CorrR1 + CorrR2 - CorrR3 - CorrR4)/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(ComplexF64,(2,2,N,nhits,T))
    @assert size(Corr2π) == size(Corrρ) == size(CorrT1) == size(CorrT2) 
    corr[1,1,:,:,:] =  @. Corrρ/L3 + 0*im
    corr[1,2,:,:,:] =  @. 0        + im*(CorrT1-CorrT2)/L3
    corr[2,1,:,:,:] =  @. 0        + im*(CorrT2-CorrT1)/L3
    corr[2,2,:,:,:] =  @. Corr2π   + 0*im
    return corr
end
function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function _parse_momentum(p0)
    rx = r"\(([0-9]),([0-9]),([0-9])\)"
    m = match(rx,p0)
    return parse.(Int,m.captures) 
end
function _are_permutations(p0,p1)
    return sort(p0) == sort(p1) 
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
    
    for ens in ensembles
        _copy_lattice_parameters(file_out,file_in;group=ens)
        T, L = fid["$ens/lattice"][1:2]
        p_external = fid["$ens/p_external"][]
        p_external_parsed = _parse_momentum.(p_external)
        @show ens
        @show p_external
        
        if "p(0,0,0)" ∈ p_external
            Corrπ0, Corrρ0 = correlatorsp000(fid,ens;p=1)
            h5write(file_out,joinpath(ens,"p(0,0,0)","correlator_pion"),Corrπ0)
            h5write(file_out,joinpath(ens,"p(0,0,0)","correlator_rho") ,Corrρ0)
        end
        if "p(0,0,1)" ∈ p_external
            Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(fid,ens;p=1)
            Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
            Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
            
            h5write(file_out,joinpath(ens,"p(0,0,1)","correlation_matrix"),Corr)
            h5write(file_out,joinpath(ens,"p(0,0,1)","correlator_pion"),Corrπ)
        end
        if "p(0,0,2)" ∈ p_external
            Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(fid,ens;p=2)
            Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
            Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
            
            h5write(file_out,joinpath(ens,"p(0,0,2)","correlation_matrix"),Corr)
            h5write(file_out,joinpath(ens,"p(0,0,2)","correlator_pion"),Corrπ)
        end
        if "p(0,1,1)" ∈ p_external
            Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp011(fid,ens;p=1)
            Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
            Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
            
            h5write(file_out,joinpath(ens,"p(0,1,1)","correlation_matrix"),Corr)
            h5write(file_out,joinpath(ens,"p(0,1,1)","correlator_pion"),Corrπ)
        end
        if "p(1,1,0)" ∈ p_external
            Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp110(fid,ens;p=1)
            Corr2π = pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
            Corr   = pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
            
            h5write(file_out,joinpath(ens,"p(1,1,0)","correlation_matrix"),Corr)
            h5write(file_out,joinpath(ens,"p(1,1,0)","correlator_pion"),Corrπ)
        end
    end
end
file_in  = "data/isospin1_merged.hdf5"
file_in2 = "data/isospin1_sorted.hdf5"
file_out = "data/isospin1_corr.hdf5"
file_out2= "data/isospin1_corr_allruns.hdf5"
write_correlation_matrix(file_in,file_out)
write_correlation_matrix(file_in2,file_out2;combined=false)
