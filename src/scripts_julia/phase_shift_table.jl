using Pkg; Pkg.activate("src/src_jl")
using HDF5
using Plots
using Statistics
using SpecialFunctions
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

function uncertainty_from_samples(samples)
    # Yannick encodes missing/failed determinations with 0
    # First, remove those samples and sort the values
    s = filter(!iszero,samples)
    sort!(s)
    # monitor how many measurements we have discarded
    # If we discard more than 10% of all samples, we do not 
    # report a result
    p = length(s)/length(samples) 
    if p < 0.9
        return NaN, NaN
    end
    # Now find the median and find the symmetric interval such that 
    # erf(1/sqrt(2)) ≈ 68% lie within this interval
    x = median(s)
    ind = last(findmin(y->abs(y-x),s))
    shift = Int(round(length(s)*erf(1/sqrt(2))/2))
    x_upper = s[ind+shift]
    x_lower = s[ind-shift]

    Δx_upper = abs(x - x_upper)
    Δx_lower = abs(x - x_lower)
    Δx = max(Δx_upper, Δx_lower)
    return x, Δx
end

function phase_shift_table_csv(infile,outfile,pattern)
    h5id = h5open(infile)
    ensembles = keys(h5id)
    io = open(outfile,"w")
    for ens in ensembles
        if endswith(ens,pattern)
            Λ = read(h5id[ens],"lattice")
            Nt, Ns = Λ[1], Λ[2]
            momenta = filter(startswith("p"),keys(h5id[ens])) 
            for mom in momenta
                irreps = filter( !isequal("pi"), keys(h5id[ens][mom]))
                for ir in irreps
                    # obtain the number of energy levels from the available data sets
                    h5channel = h5id[ens][mom][ir]
                    # read fitted energy levels
                    E = read(h5channel,"E")
                    ΔE = read(h5channel,"Delta_E")
                    levels = filter(startswith("lv"),keys(h5channel))
                    for (i,l) in enumerate(levels)
                        # We need to work out the uncertainties 
                        s_sample = read(h5channel[l]["sample"],"s")
                        δ_sample = read(h5channel[l]["sample"],"PS")
                        δ, Δδ = uncertainty_from_samples(real.(δ_sample))
                        rs, Δrs = uncertainty_from_samples(real.(sqrt.(s_sample)))
                        println(io,"$ens;$Nt;$Ns;$mom;$ir;$i;$(E[i]);$(ΔE[i]);$(rs);$(Δrs);$δ;$Δδ")
                    end
                end
            end
        end
    end
    close(io)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_in"
        help = "HDF5 file containing the parsed data"
        required = true
        "--outfile"
        help = "CSV output containing the list of all parsed runs"
        required = true
        "--pattern"
        help = "label of the ensemble in the form 'beta{coupling}m{mass}' "
        default = nothing
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    phase_shift_table_csv(args["h5file_in"],args["outfile"],args["pattern"])
end 
main()