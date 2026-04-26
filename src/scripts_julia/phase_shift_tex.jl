using Pkg; Pkg.activate("src/src_jl")
using LatticeUtils
using DelimitedFiles
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

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
end 

metadata = readdlm("metadata/fit_scatter_input.csv",';',skipstart=0)
csv_data = readdlm("tmp/beta6.9m-0.92.csv",';',skipstart=1)

for r in eachrow(csv_data)
    ens, Nt, Ns, mom, ir, lv, E, ΔE, rs, Δrs, δ, Δδ = r[1:end]
    d2  = length(findall('1',mom))
    key = ens*mom*ir*"lv$(lv-1)"
    incl = startswith(ens,"beta7.05m-0.863") ? false : key ∈ metadata
    incl = incl ? "yes" : "no" 
    println("$Ns & $d2 & $ir & $lv & $(errorstring(E, ΔE; nsig=1)) & $(errorstring(rs, Δrs; nsig=1)) & $(errorstring(δ, Δδ; nsig=1)) & $incl")
end