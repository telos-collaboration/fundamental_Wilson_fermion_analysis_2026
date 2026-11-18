using Pkg; Pkg.activate("scattering/src/src_jl")
using HDF5: h5open, h5write
using ScatteringI1

h5file = "data_assets/scattering/isospin1_eigenvalues.hdf5"
h5id = h5open(h5file)
ens_id = 1
mom_id = 1
ens = keys(h5id)
mom = filter(startswith("p("),keys(h5id[ens[ens_id]]))
# filter out vanishing momenta
mom = filter(!isequal("p(0,0,0)"),mom)

# obtain the correlation matrices
Corr2x2 = read(h5id[ens[ens_id]][mom[mom_id]]["A1"],"Corr2x2")
Corr3x3 = read(h5id[ens[ens_id]][mom[mom_id]]["A1"],"Corr3x3")

# perform GEVP now with fixed t-t0 such that t0 > t/2
t0 = 1
swap = false
gevp = true
deriv = false
swap_t = 0
binsize = 2
symmetrise = true
gevp_fixed_t0_diff = true
eigvals_resamples = ScatteringI1.variational_analysis_samples(Corr2x2;t0,deriv,gevp,symmetrise,swap,swap_t,binsize,gevp_fixed_t0_diff)
