using Pkg; Pkg.activate("scattering/src/src_jl")
using LatticeUtils
using DelimitedFiles
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

function latex_table(datafile,metadatafile,outfile)
    metadata = readdlm(metadatafile,';',skipstart=0)
    csv_data = readdlm(datafile,';',skipstart=1)

    header = raw"""
    \begin{tabular}[t]{|c|c|c|c|c|c|c|c|}
        \hline 
        $N_L$ & $|\vec{d}|^2$ & $\Lambda$ & $n$ & $aE$ & $a\sqrt{s}$ & $\delta_1$ & Incl. \\ \hline \hline"""
    footer = raw"""
        \hline
    \end{tabular}"""

    io = open(outfile,"w")
    counter = 0
    println(io,header)
    for (i,r) in enumerate(eachrow(csv_data))
        ens, Nt, Ns, mom, ir, lv, E, ΔE, rs, Δrs, δ, Δδ = r[1:end]
        d2  = length(findall('1',mom))
        key = ens*mom*ir*"lv$(lv-1)"
        # We perform no fit for the medium ensemble
        # Thus, we always denote the phase shifts as not included
        incl = startswith(ens,"beta7.05m-0.863") ? false : key ∈ metadata
        incl = incl ? "yes" : "no" 
        # see if the spatial volume changes
        # once it has changed twice, create a new tabular environment
        if i == 1 ||  csv_data[i-1,3] != csv_data[i,3]
            counter += 1 
            if counter == 3
                println(io,footer)
                println(io,"\\quad")
                println(io,header)
            elseif counter > 1
                println(io,"    \\hline")
            end
        end
        println(io,"    $Ns & $d2 & $ir & $lv & $(errorstring(E, ΔE; nsig=1)) & $(errorstring(rs, Δrs; nsig=1)) & $(errorstring(δ, Δδ; nsig=1)) & $incl \\\\")
    end
    println(io,footer)
    close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--csv_in"
        help = "HDF5 file containing the parsed data"
        required = true
        "--outfile"
        help = "CSV output containing the list of all parsed runs"
        required = true
        "--metadata"
        help = "label of the ensemble in the form 'beta{coupling}m{mass}' "
        default = nothing
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    latex_table(args["csv_in"],args["metadata"],args["outfile"])
end 
main()
