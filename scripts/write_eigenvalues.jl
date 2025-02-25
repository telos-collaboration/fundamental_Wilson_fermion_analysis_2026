using Pkg; Pkg.activate(".")
using ScatteringI1
using HDF5
using Statistics
using LatticeUtils
using ProgressMeter

function swap_eigval_numbering(old,t0,T)
    new = copy(old)
    @. new[1,:,1:t0-1] = old[2,:,1:t0-1]
    @. new[2,:,1:t0-1] = old[1,:,1:t0-1]
    @. new[1,:,T-t0+2:T] = old[2,:,T-t0+2:T]
    @. new[2,:,T-t0+2:T] = old[1,:,T-t0+2:T]
    return new
end
function variational_analysis(Corr;t0,maxhits=typemax(Int),deriv=true)

    nhits, T = size(Corr)[4:5]
    h     = min(nhits,maxhits)
    Corr  = dropdims(mean(Corr[:,:,:,1:h,:],dims=4),dims=4)
    Corr  = correlator_folding(Corr;t_dim=4,sign=+1)

    if deriv
        Corr = correlator_derivative(Corr;t_dim=4)
    end

    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0, T)
    eigvals, Δeigvals = LatticeUtils.apply_jackknife(eigvals_resamples;dims=2)
    eigvals_cov = LatticeUtils.cov_jackknife_eigenvalues(eigvals_resamples)

    return eigvals, Δeigvals, eigvals_cov
end
function _copy_lattice_parameters(outfile,infile;group="")
    file = h5open(infile)[group]
    entries = filter(!contains(r"p\([0-9],[0-9],[0-9]\)") ,keys(file))
    for entry in entries
        label = joinpath(group,entry)
        h5write(outfile,label,read(file,entry))
    end
end
function write_all_eigenvalues(infile,outfile; t0, deriv, maxhits=typemax(Int))
    
    h5dset   = h5open(infile)
    isfile(outfile) && rm(outfile)

    ensembles = keys(h5dset)
    @showprogress for ens in ensembles
        _copy_lattice_parameters(outfile,infile;group=ens)
        p_external = h5dset["$ens/p_external"][]
        for p in p_external
            p == "p(0,0,0)" && continue
            Corr = h5dset[joinpath(ens,p,"correlation_matrix")][]
            eigvals, Δeigvals, eigvals_cov = variational_analysis(Corr;t0,maxhits,deriv)
            eigvals, Δeigvals = real.(eigvals), real.(Δeigvals), real.(eigvals_cov)

            h5write(outfile,joinpath(ens,p,"eigvals"),eigvals)
            h5write(outfile,joinpath(ens,p,"Delta_eigvals"),Δeigvals)
            h5write(outfile,joinpath(ens,p,"cov_eigvals"),eigvals_cov)
        end
    end
end

outfile = "data/isospin1_eigenvalues_t0_3_deriv.hdf5"
infile  = "data/isospin1_corr.hdf5"
t0      = 3
deriv   = true
write_all_eigenvalues(infile,outfile; t0, deriv)