using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
plotlyjs()
function _average_correlator(Corr)
    N, nhits, T = size(Corr)
    CorrAvg = dropdims(mean(Corr,dims=2),dims=2)
    C  = dropdims(mean(CorrAvg,dims=1),dims=1)
    ΔC = dropdims(std(CorrAvg,dims=1),dims=1)./sqrt(N)
    return C, ΔC
end

# Lt32Ls24beta6.9m1-0.92m2-0.92
hdf5file = "isospin1_raw.hdf5"
h5dset = h5open(hdf5file)

Corr = h5dset["E1/p(0,0,1)/r1/p_diag(0,0,1)/C_re"][]
Corr = h5dset["E1/p(0,0,1)/d/p_diag(0,0,1)/C_re"][]
Corr = h5dset["E1/p(0,0,1)/pi/p_diag(0,0,1)/C_re"][]
Corr = h5dset["E1/p(0,0,1)/rho_g33/p_diag(0,0,1)/C_re"][]
Corr = h5dset["E1/p(1,1,0)/t1_g3/p_diag(1,1,0)/C_im"][]

C_re, ΔC_re = _average_correlator(Corr)

plt = plot()
scatter!(plt,C_re,yerr=ΔC_re,yscale=:log10)

