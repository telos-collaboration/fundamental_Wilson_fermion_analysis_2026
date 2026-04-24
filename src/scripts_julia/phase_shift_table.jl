using Pkg; Pkg.activate("src/src_jl")
using HDF5

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
                end
            end
        end
        break
    end
end