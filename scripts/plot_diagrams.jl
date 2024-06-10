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

hdf5file = "isospin1_h16_raw.hdf5"
h5dset = h5open(hdf5file)

CorrD  = h5dset["E1/p(0,0,1)/d/p_diag(0,0,1)/C_re"][]
Corrπ  = h5dset["E1/p(0,0,1)/pi/p_diag(0,0,1)/C_re"][]
Corrρ  = h5dset["E1/p(0,0,1)/rho_g33/p_diag(0,0,1)/C_re"][]
CorrT1 = -h5dset["E1/p(0,0,1)/t1_g3/p_diag(0,0,1)/C_im"][]
CorrT2 =  h5dset["E1/p(0,0,1)/t2_g3/p_diag(0,0,1)/C_im"][]

CorrR1 = h5dset["E1/p(0,0,1)/r1/p_diag(0,0,1)/C_re"][]
CorrR2 = h5dset["E1/p(0,0,1)/r2/p_diag(0,0,1)/C_re"][]
CorrR3 = h5dset["E1/p(0,0,1)/r3/p_diag(0,0,1)/C_im"][]
CorrR4 = h5dset["E1/p(0,0,1)/r4/p_diag(0,0,1)/C_im"][]

N, nhits, T = size(CorrD)

CD,  ΔCD  = _average_correlator(CorrD)
Cπ,  ΔCπ  = _average_correlator(Corrπ)
Cρ,  ΔCρ  = _average_correlator(Corrρ)
CT1, ΔCT1 = _average_correlator(CorrT1)
CT2, ΔCT2 = _average_correlator(CorrT2)

CR1, ΔCR1 = _average_correlator(CorrR1)
CR2, ΔCR2 = _average_correlator(CorrR2)
CR3, ΔCR3 = _average_correlator(CorrR3)
CR4, ΔCR4 = _average_correlator(CorrR4)

pltPiPi = plot()
pltMes  = plot()
pltCross = plot()

scatter!(pltMes,Cπ, yerr=ΔCπ,yscale=:log10, label="(π) n_src = $nhits")
scatter!(pltMes,Cρ, yerr=ΔCρ,yscale=:log10, label="(ρ) n_src = $nhits")

scatter!(pltCross,CT1,yerr=ΔCT1,yscale=:log10,label="-Im(T1) n_src = $nhits")
scatter!(pltCross,CT2,yerr=ΔCT2,yscale=:log10,label="+Im(T2) n_src = $nhits")

scatter!(pltPiPi,CD, yerr=ΔCD ,yscale=:log10,label="(D) n_src = $nhits")
# These correlators dip below zero: Study nhit-dependence
# Also, study autocorrelation, topological freezing and related
#scatter!(pltPiPi,CR1,yerr=ΔCR1,yscale=:log10,label="(R1) n_src = $nhits")
#scatter!(pltPiPi,CR2,yerr=ΔCR2,yscale=:log10,label="(R2) n_src = $nhits")
#scatter!(pltPiPi,CR3,yerr=ΔCR3,label="(R3) n_src = $nhits")
#scatter!(pltPiPi,CR4,yerr=ΔCR4,label="(R4) n_src = $nhits")

pltMes
pltCross
pltPiPi