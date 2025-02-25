using Pkg; Pkg.activate(".")
using ScatteringI1
using HDF5
using Statistics
using LatticeUtils
using ProgressMeter

function variational_analysis(Corr;t0,maxhits=typemax(Int),deriv=true)

    nhits = size(Corr)[4]
    h     = min(nhits,maxhits)
    Corr  = dropdims(mean(Corr[:,:,:,1:h,:],dims=4),dims=4)
    Corr  = correlator_folding(Corr;t_dim=4,sign=+1)

    if deriv
        Corr = correlator_derivative(Corr;t_dim=4)
    end

    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    eigvals, Δeigvals = LatticeUtils.apply_jackknife(eigvals_resamples;dims=2)
    eigvals_cov = LatticeUtils.cov_jackknife_eigenvalues(eigvals_resamples)

    return eigvals, Δeigvals, eigvals_cov
end
function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"(APE|Wuppertal)") ,keys(file))
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,entries)
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end

outfile  = "data/isospin1_eigenvalues_t0_8.hdf5"
hdf5file = "data/isospin1_corr.hdf5"
h5dset   = h5open(hdf5file)
maxhits  = typemax(Int)
t0       = 8

isfile(outfile) && rm(outfile)

ensembles = keys(h5dset)
@showprogress for ens in ensembles
    _copy_lattice_parameters(outfile,hdf5file;group=ens)
    p_external = h5dset["$ens/p_external"][]
    for p in p_external
        p == "p(0,0,0)" && continue
        Corr = h5dset[joinpath(ens,p,"correlation_matrix")][]
        eigvals, Δeigvals, eigvals_cov = variational_analysis(Corr;t0,maxhits,deriv=true)
        eigvals, Δeigvals = real.(eigvals), real.(Δeigvals)

        h5write(outfile,joinpath(ens,p,"eigvals"),eigvals)
        h5write(outfile,joinpath(ens,p,"Delta_eigvals"),Δeigvals)
        h5write(outfile,joinpath(ens,p,"cov_eigvals"),eigvals_cov)
    end
end
