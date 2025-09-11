using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using DelimitedFiles: readdlm
using HDF5: h5open, h5write
using LatticeUtils: log_meff
using ScatteringI1
using Statistics
include("utils_swap.jl")

function _copy_lattice_parameters_eigenvalues(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        entry == "p_external" && continue
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end

function mean_error_cov(corr)
    me = dropdims(mean(corr, dims=1);dims=1)
    sd = dropdims(std(corr, dims=1);dims=1)
    cv = cov(corr, dims=1)
    sd .= sd ./ sqrt(size(corr,1))
    cv .= cv ./ size(corr,1)
    return me, sd, cv
end

function write_meson_correlators(outfile,ens,p,irrep,Corr)
    meff, Δmeff = log_meff(Corr')
    C, ΔC, Ccov = mean_error_cov(Corr)
    exists_file = isfile(outfile)
    if exists_file
        fid = h5open(outfile)
        exists_group = haskey(fid,joinpath(ens,p,irrep,"C"))
        close(fid)
        if exists_group
            return 
        end
    end
    h5write(outfile,joinpath(ens,p,irrep,"C"),C)
    h5write(outfile,joinpath(ens,p,irrep,"Delta_C"),ΔC)
    h5write(outfile,joinpath(ens,p,irrep,"cov_C"),Ccov)
    h5write(outfile,joinpath(ens,p,irrep,"meff"),meff)
    h5write(outfile,joinpath(ens,p,irrep,"Delta_meff"),Δmeff)    
end

function write_all_eigenvalues(infile,outfile; maxhits=typemax(Int), average_equivalent_momenta=true, metadata, swap_metadata)    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    # get metadate for specific momentum
    data = readdlm(metadata,',',skipstart=1)
    ensembles = unique(data[:,1])

    for ens in ensembles
        # copy basic information for every ensemble 
        _copy_lattice_parameters_eigenvalues(outfile,infile;group=ens)    
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        h5write(outfile,"$ens/p_external",p_external)

        # get correlators at vanishing momentum for every ensemble
        p = "p(0,0,0)"
        Corrπ = read_meson_correlator(h5dset,ens,p,"correlator_pion";average_equivalent_momenta)
        write_meson_correlators(outfile,ens,p,"pi",Corrπ)        
        Corrρ = read_meson_correlator(h5dset,ens,p,"correlator_rho";average_equivalent_momenta)
        write_meson_correlators(outfile,ens,p,"T1",Corrρ)
    end

    @showprogress desc="Write eigenvalues:" enabled=true for row in eachrow(data)

        ens, p, id = row[1], row[2], row[14]
        t0, deriv, gevp, symmetrise = Int(row[8]), Bool(row[9]), Bool(row[10]), Bool(row[11])

        Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits,average_equivalent_momenta)    
        # For the 2x2 problem the eigenvalues should always be relabelled (swapped) at swap_t = t0 when using the GEVP
        eigvals, Δeigvals, eigvals_cov = ScatteringI1.variational_analysis(Corr;t0,deriv,gevp,symmetrise,swap=gevp,swap_t=t0)
        meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv,gevp,symmetrise,swap=gevp,swap_t=t0)
        eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)

        three_by_three = haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
        if three_by_three
            Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits,average_equivalent_momenta)
            Corr3x3[1:2,1:2,:,:] .= Corr
            eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3 = ScatteringI1.variational_analysis(Corr3x3;t0,deriv,gevp,symmetrise)
            eigvals_3x3, Δeigvals_3x3 = real.(eigvals_3x3), real.(Δeigvals_3x3), real.(eigvals_cov_3x3)
            if !isempty(swap_metadata)
                eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3 = swap_eigvals(eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3, swap_metadata, ens, p, "A1", id)
            end
            meff_3x3, Δmeff_3x3 = ScatteringI1.effective_masses(Corr3x3;t0,deriv,gevp,symmetrise)
        end

        Corrπ = read_meson_correlator(h5dset,ens,p,"correlator_pion";average_equivalent_momenta)
        write_meson_correlators(outfile,ens,p,"pi",Corrπ)
        if haskey(h5dset,joinpath(ens,p,"B1"))
            B1 = dropdims(mean(read(h5dset,joinpath(ens,p,"B1")),dims=2);dims=2)
            write_meson_correlators(outfile,ens,p,"B1",B1)
        end 
        if haskey(h5dset,joinpath(ens,p,"E"))
            E = dropdims(mean(read(h5dset,joinpath(ens,p,"E")),dims=2);dims=2)
            write_meson_correlators(outfile,ens,p,"E",E)
        end 
   
        h5write(outfile,joinpath(ens,p,"A1",id,"eigvals"),eigvals)
        h5write(outfile,joinpath(ens,p,"A1",id,"Delta_eigvals"),Δeigvals)
        h5write(outfile,joinpath(ens,p,"A1",id,"cov_eigvals"),eigvals_cov)
        h5write(outfile,joinpath(ens,p,"A1",id,"t0"),t0)
        h5write(outfile,joinpath(ens,p,"A1",id,"gevp"),gevp)
        h5write(outfile,joinpath(ens,p,"A1",id,"deriv"),deriv)
        h5write(outfile,joinpath(ens,p,"A1",id,"symmetrise"),symmetrise)
        h5write(outfile,joinpath(ens,p,"A1",id,"average_equivalent_momenta"),average_equivalent_momenta)
        h5write(outfile,joinpath(ens,p,"A1",id,"momenta"),ascii(momenta))
        h5write(outfile,joinpath(ens,p,"A1",id,"sources"),[sources...])
        h5write(outfile,joinpath(ens,p,"A1",id,"Corr2x2"),Corr)
        h5write(outfile,joinpath(ens,p,"A1",id,"meff"),meff)
        h5write(outfile,joinpath(ens,p,"A1",id,"Delta_meff"),Δmeff)            
        if three_by_three
            h5write(outfile,joinpath(ens,p,"A1",id,"Corr3x3"),Corr3x3)
            h5write(outfile,joinpath(ens,p,"A1",id,"momenta_3x3"),ascii(momenta3x3))
            h5write(outfile,joinpath(ens,p,"A1",id,"sources_3x3"),[sources3x3...])
            h5write(outfile,joinpath(ens,p,"A1",id,"eigvals_3x3"),eigvals_3x3)
            h5write(outfile,joinpath(ens,p,"A1",id,"Delta_eigvals_3x3"),Δeigvals_3x3)
            h5write(outfile,joinpath(ens,p,"A1",id,"cov_eigvals_3x3"),eigvals_cov_3x3)
            h5write(outfile,joinpath(ens,p,"A1",id,"meff_3x3"),meff_3x3)
            h5write(outfile,joinpath(ens,p,"A1",id,"Delta_meff_3x3"),Δmeff_3x3)            
        end
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_in"
        help = "HDF5 file containing the parsed data"
        required = true
        "--h5file_out"
        help = "HDF5 output file containing the correlation matrices"
        required = true
        "--metadata"
        help = "CSV file containing the parameters for the variational analysis"
        required = true
        "--avg"
        help = "Average over equivalent external momenta"
        required = true
        arg_type = Bool
        "--maxhits"
        help = "Maximal number of stochastic sources to include"
        default = typemax(Int)
        "--swap_metadata"
        help = "CSV containing data for eigenvalue relabelling"
        default = "" 
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    maxhits = args["maxhits"]
    metadata = args["metadata"]
    swap_metadata = args["swap_metadata"]
    average_equivalent_momenta = args["avg"]
    write_all_eigenvalues(args["h5file_in"],args["h5file_out"]; maxhits, average_equivalent_momenta, metadata, swap_metadata)
end
main()

