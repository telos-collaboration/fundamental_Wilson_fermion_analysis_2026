using Pkg; Pkg.activate("src/src_jl")
using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using ProgressMeter: @showprogress
using HDF5: h5open, h5write
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

function mean_and_error(corr)
    me   = dropdims(mean(corr, dims=1);dims=1)
    sd   = dropdims(std(corr, dims=1);dims=1)
    err  = sd ./ sqrt(size(corr,1))
    return me, err
end

function write_all_eigenvalues(infile,outfile; t0, deriv, maxhits=typemax(Int), average_equivalent_momenta, gevp, symmetrise)    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    ensembles = keys(h5dset)
    @showprogress desc="Write eigenvalues:" enabled=true for ens in ensembles
        _copy_lattice_parameters_eigenvalues(outfile,infile;group=ens)
        p0 = read(h5dset,"$ens/p_external")
        p_external = ifelse(average_equivalent_momenta,unique_momenta(p0),p0)
        h5write(outfile,"$ens/p_external",p_external)
        for p in p_external

            Corrπ = read_pion_correlator(h5dset,ens,p;average_equivalent_momenta)
            Cπ, ΔCπ = mean_and_error(Corrπ)
            h5write(outfile,joinpath(ens,p,"correlator_pion"),Corrπ)
            h5write(outfile,joinpath(ens,p,"Cpi"),Cπ)
            h5write(outfile,joinpath(ens,p,"Delta_Cpi"),ΔCπ)
            p == "p(0,0,0)" && continue

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
                CB1, ΔCB1 = mean_and_error(B1)
                h5write(outfile,joinpath(ens,p,"correlator_B1"),B1)
                h5write(outfile,joinpath(ens,p,"B1","C"),CB1)
                h5write(outfile,joinpath(ens,p,"B1","Delta_C"),ΔCB1)
            end 
            if haskey(h5dset,joinpath(ens,p,"E"))
                E = dropdims(mean(read(h5dset,joinpath(ens,p,"E")),dims=2);dims=2)
                CE, ΔCE = mean_and_error(E)
                h5write(outfile,joinpath(ens,p,"correlator_E"),E)
                h5write(outfile,joinpath(ens,p,"E","C"),CE)
                h5write(outfile,joinpath(ens,p,"E","Delta_C"),ΔCE)
            end 
   
            h5write(outfile,joinpath(ens,p,"A1","eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"A1","Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"A1","cov_eigvals"),eigvals_cov)
            h5write(outfile,joinpath(ens,p,"t0"),t0)
            h5write(outfile,joinpath(ens,p,"gevp"),gevp)
            h5write(outfile,joinpath(ens,p,"deriv"),deriv)
            h5write(outfile,joinpath(ens,p,"symmetrise"),symmetrise)
            h5write(outfile,joinpath(ens,p,"average_equivalent_momenta"),average_equivalent_momenta)
            h5write(outfile,joinpath(ens,p,"momenta"),ascii(momenta))
            h5write(outfile,joinpath(ens,p,"sources"),[sources...])
            h5write(outfile,joinpath(ens,p,"Corr2x2"),Corr)
            h5write(outfile,joinpath(ens,p,"meff"),meff)
            h5write(outfile,joinpath(ens,p,"Delta_meff"),Δmeff)            
            if three_by_three
                h5write(outfile,joinpath(ens,p,"Corr3x3"),Corr3x3)
                h5write(outfile,joinpath(ens,p,"momenta_3x3"),ascii(momenta3x3))
                h5write(outfile,joinpath(ens,p,"sources_3x3"),[sources3x3...])
                h5write(outfile,joinpath(ens,p,"eigvals_3x3"),eigvals_3x3)
                h5write(outfile,joinpath(ens,p,"Delta_eigvals_3x3"),Δeigvals_3x3)
                h5write(outfile,joinpath(ens,p,"cov_eigvals_3x3"),eigvals_cov_3x3)
                h5write(outfile,joinpath(ens,p,"meff_3x3"),meff_3x3)
                h5write(outfile,joinpath(ens,p,"Delta_meff_3x3"),Δmeff_3x3)            
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
        "--gevp"
        help = "Use GEVP analysis, otherwise use a simple EVP"
        required = true
        arg_type = Bool
        "--t0"
        help = "Reference time t0 used in GEVP analysis"
        default = 1
        arg_type = Int
        "--deriv"
        help = "Perform a numerical derivative first"
        required = true
        arg_type = Bool
        "--avg"
        help = "Average over equivalent external momenta"
        required = true
        arg_type = Bool
        "--maxhits"
        help = "Maximal number of stochastic sources to include"
        default = typemax(Int)
        "--symmetrise"
        help = "Symmetrise correlation matrix with resprect to T/2"
        required = true
        arg_type = Bool
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    t0 = args["t0"]
    gevp = args["gevp"]
    deriv = args["deriv"]
    maxhits = args["maxhits"]
    symmetrise = args["symmetrise"]
    average_equivalent_momenta = args["avg"]
    write_all_eigenvalues(args["h5file_in"],args["h5file_out"]; gevp, t0, deriv, maxhits, average_equivalent_momenta, symmetrise)    
end
main()

