using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using DelimitedFiles: readdlm
using HDF5: h5open, h5write
using LatticeUtils: log_meff
using ScatteringI1
using Statistics

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

function write_all_eigenvalues(infile,outfile; maxhits=typemax(Int), average_equivalent_momenta=true, metadata)    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    ensembles = keys(h5dset)
    @showprogress desc="Write eigenvalues:" enabled=true for ens in ensembles

        _copy_lattice_parameters_eigenvalues(outfile,infile;group=ens)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        h5write(outfile,"$ens/p_external",p_external)

        for p in p_external

            p == "p(0,0,0)" && continue
            # TODO: Deal with p(0,0,0)

            # get metadate for specific momentum
            data = readdlm(metadata,',',skipstart=1)
            metadata_ind = findfirst(i -> isequal(joinpath(ens,p),joinpath(data[i,1:2]...)),1:first(size(data)))
            t0, deriv, gevp, symmetrise = data[metadata_ind,8:11]
            t0::Int, deriv::Bool, gevp::Bool, symmetrise::Bool = Int(t0), Bool(deriv), Bool(gevp), Bool(symmetrise) 

            Corrπ = read_pion_correlator(h5dset,ens,p;average_equivalent_momenta)
            meffπ, Δmeffπ = log_meff(Corrπ')
            Cπ, ΔCπ, Cπcov = mean_error_cov(Corrπ)
            h5write(outfile,joinpath(ens,p,"correlator_pion"),Corrπ)
            h5write(outfile,joinpath(ens,p,"Cpi"),Cπ)
            h5write(outfile,joinpath(ens,p,"Delta_Cpi"),ΔCπ)
            h5write(outfile,joinpath(ens,p,"cov_Cpi"),Cπcov)
            h5write(outfile,joinpath(ens,p,"meff_pi"),meffπ)
            h5write(outfile,joinpath(ens,p,"Delta_meff_pi"),Δmeffπ)


            Corr, sources, momenta = read_correlation_matrix(h5dset,ens,p,"correlation_matrix";maxhits,average_equivalent_momenta)    
            eigvals, Δeigvals, eigvals_cov = ScatteringI1.variational_analysis(Corr;t0,deriv,gevp,symmetrise)
            eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)
            meff, Δmeff = ScatteringI1.effective_masses(Corr;t0,deriv,gevp,symmetrise)

            three_by_three = haskey(h5dset[ens][p],"correlation_matrix_3x3_ext")
            if three_by_three
                Corr3x3, sources3x3, momenta3x3 = read_correlation_matrix(h5dset,ens,p,"correlation_matrix_3x3_ext";maxhits,average_equivalent_momenta)
                Corr3x3[1:2,1:2,:,:] .= Corr
                eigvals_3x3, Δeigvals_3x3, eigvals_cov_3x3 = ScatteringI1.variational_analysis(Corr3x3;t0,deriv,gevp,symmetrise)
                eigvals_3x3, Δeigvals_3x3 = real.(eigvals_3x3), real.(Δeigvals_3x3), real.(eigvals_cov_3x3)
                meff_3x3, Δmeff_3x3 = ScatteringI1.effective_masses(Corr3x3;t0,deriv,gevp,symmetrise)
            end

            if haskey(h5dset,joinpath(ens,p,"B1"))
                B1 = dropdims(mean(read(h5dset,joinpath(ens,p,"B1")),dims=2);dims=2)
                meffB1, ΔmeffB1 = log_meff(B1')
                CB1, ΔCB1, covCB1 = mean_error_cov(B1)
                h5write(outfile,joinpath(ens,p,"B1","correlator_B1"),B1)
                h5write(outfile,joinpath(ens,p,"B1","C"),CB1)
                h5write(outfile,joinpath(ens,p,"B1","Delta_C"),ΔCB1)
                h5write(outfile,joinpath(ens,p,"B1","cov_C"),covCB1)
                h5write(outfile,joinpath(ens,p,"B1","meff"),meffB1)
                h5write(outfile,joinpath(ens,p,"B1","Delta_meff"),ΔmeffB1)
            end 
            if haskey(h5dset,joinpath(ens,p,"E"))
                E = dropdims(mean(read(h5dset,joinpath(ens,p,"E")),dims=2);dims=2)
                meffE, ΔmeffE = log_meff(E')
                CE, ΔCE, covCE = mean_error_cov(E)
                h5write(outfile,joinpath(ens,p,"E","correlator_E"),E)
                h5write(outfile,joinpath(ens,p,"E","C"),CE)
                h5write(outfile,joinpath(ens,p,"E","Delta_C"),ΔCE)
                h5write(outfile,joinpath(ens,p,"E","cov_C"),covCE)
                h5write(outfile,joinpath(ens,p,"E","meff"),meffE)
                h5write(outfile,joinpath(ens,p,"E","Delta_meff"),ΔmeffE)
            end 
   
            h5write(outfile,joinpath(ens,p,"A1","eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"A1","Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"A1","cov_eigvals"),eigvals_cov)
            h5write(outfile,joinpath(ens,p,"A1","t0"),t0)
            h5write(outfile,joinpath(ens,p,"A1","gevp"),gevp)
            h5write(outfile,joinpath(ens,p,"A1","deriv"),deriv)
            h5write(outfile,joinpath(ens,p,"A1","symmetrise"),symmetrise)
            h5write(outfile,joinpath(ens,p,"A1","average_equivalent_momenta"),average_equivalent_momenta)
            h5write(outfile,joinpath(ens,p,"A1","momenta"),ascii(momenta))
            h5write(outfile,joinpath(ens,p,"A1","sources"),[sources...])
            h5write(outfile,joinpath(ens,p,"A1","Corr2x2"),Corr)
            h5write(outfile,joinpath(ens,p,"A1","meff"),meff)
            h5write(outfile,joinpath(ens,p,"A1","Delta_meff"),Δmeff)            
            if three_by_three
                h5write(outfile,joinpath(ens,p,"A1","Corr3x3"),Corr3x3)
                h5write(outfile,joinpath(ens,p,"A1","momenta_3x3"),ascii(momenta3x3))
                h5write(outfile,joinpath(ens,p,"A1","sources_3x3"),[sources3x3...])
                h5write(outfile,joinpath(ens,p,"A1","eigvals_3x3"),eigvals_3x3)
                h5write(outfile,joinpath(ens,p,"A1","Delta_eigvals_3x3"),Δeigvals_3x3)
                h5write(outfile,joinpath(ens,p,"A1","cov_eigvals_3x3"),eigvals_cov_3x3)
                h5write(outfile,joinpath(ens,p,"A1","meff_3x3"),meff_3x3)
                h5write(outfile,joinpath(ens,p,"A1","Delta_meff_3x3"),Δmeff_3x3)            
            end
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
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    maxhits = args["maxhits"]
    metadata = args["metadata"]
    average_equivalent_momenta = args["avg"]
    write_all_eigenvalues(args["h5file_in"],args["h5file_out"]; maxhits, average_equivalent_momenta, metadata)    
end
main()

