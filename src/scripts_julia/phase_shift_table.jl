using Pkg; Pkg.activate("src/src_jl")
using HDF5
using Plots
using Statistics
using SpecialFunctions

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
    # If the upper and lower bound are very asymmetric, we will discard
    # them, since the underlying distribution appears to be non-Gaussian
    Δx_upper = abs(x - x_upper)
    Δx_lower = abs(x - x_lower)
    if (1/2) < Δx_upper/Δx_lower < 2 
        Δx = max(Δx_upper, Δx_lower)
        return x, Δx
    end
    return NaN, NaN
end

file = "data_assets/isospin1_fit_scatter.hdf5"
h5id = h5open(file)
ensembles = keys(h5id)

@show file

is_heavy(ensemble) = endswith(ensemble,"beta6.9m-0.92")
is_medium(ensemble) = endswith(ensemble,"beta7.05m-0.863")
is_light(ensemble) = endswith(ensemble,"beta7.05m-0.867")

for ens in ensembles
    if is_heavy(ens)
        Λ = read(h5id[ens],"lattice")
        Nt, Ns = Λ[1], Λ[2]
        momenta = filter(startswith("p"),keys(h5id[ens])) 
        for mom in momenta
            irreps = filter( !isequal("pi"), keys(h5id[ens][mom]))
            for ir in irreps
                # obtain the number of energy levels from the available data sets
                h5channel = h5id[ens][mom][ir]
                levels = filter(startswith("lv"),keys(h5channel))
                for l in levels
                    # We need to work out the uncertainties 
                    s_sample = read(h5channel[l]["sample"],"s")
                    δ_sample = read(h5channel[l]["sample"],"PS")
                    δ, Δδ = uncertainty_from_samples(real.(δ_sample))
                    s, Δs = uncertainty_from_samples(real.(s_sample))
                end
            end
        end
    end
end