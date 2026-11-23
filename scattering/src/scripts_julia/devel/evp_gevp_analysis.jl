using Pkg; Pkg.activate("scattering/src/src_jl")
using HDF5: h5open, h5write
using ScatteringI1
using LatticeUtils
using Plots
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.9,size=(400,300),titlefontsize=11)

h5file = "data_assets/scattering/isospin1_fitresults_evp_gevp.hdf5"
h5id = h5open(h5file)
ensembles = keys(h5id)
ens_id = 1
fit_id = 1
mom_id = 1
fits = filter(contains("evp"),keys(h5id[ensembles[ens_id]]))
moms = filter(!isequal("p(0,0,0)"),keys(h5id[ensembles[ens_id]][fits[1]]))

E  = hcat([read(h5id[ensembles[ens_id]][f][moms[mom_id]],"A1/E") for f in fits]...)
ΔE = hcat([read(h5id[ensembles[ens_id]][f][moms[mom_id]],"A1/Delta_E") for f in fits]...)

using Plots
scatter(E',yerr=ΔE')