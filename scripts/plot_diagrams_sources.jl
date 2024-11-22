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

hdf5file = "isospin1.hdf5"
h5dset = h5open(hdf5file)

CorrD_4  =  h5dset["T32_L24_h4_p1/p(0,0,1)/d/p_diag(0,0,1)/C_re"][]
Corrπ_4  =  h5dset["T32_L24_h4_p1/p(0,0,1)/pi/p_diag(0,0,1)/C_re"][]
Corrρ_4  =  h5dset["T32_L24_h4_p1/p(0,0,1)/rho_g33/p_diag(0,0,1)/C_re"][]
CorrT1_4 = -h5dset["T32_L24_h4_p1/p(0,0,1)/t1_g3/p_diag(0,0,1)/C_im"][]
CorrT2_4 =  h5dset["T32_L24_h4_p1/p(0,0,1)/t2_g3/p_diag(0,0,1)/C_im"][]
CorrR1_4 =  h5dset["T32_L24_h4_p1/p(0,0,1)/r1/p_diag(0,0,1)/C_re"][]
CorrR2_4 =  h5dset["T32_L24_h4_p1/p(0,0,1)/r2/p_diag(0,0,1)/C_re"][]
CorrR3_4 =  h5dset["T32_L24_h4_p1/p(0,0,1)/r3/p_diag(0,0,1)/C_re"][]
CorrR4_4 =  h5dset["T32_L24_h4_p1/p(0,0,1)/r4/p_diag(0,0,1)/C_re"][]

N, nhits, T = size(CorrD_4)
CD_4,  ΔCD_4   = _average_correlator(CorrD_4)
Cπ_4,  ΔCπ_4   = _average_correlator(Corrπ_4)
Cρ_4,  ΔCρ_4   = _average_correlator(Corrρ_4)
CT1_4, ΔCT1_4  = _average_correlator(CorrT1_4)
CT2_4, ΔCT2_4  = _average_correlator(CorrT2_4)
CR1_4, ΔCR1_4 = _average_correlator(CorrR1_4)
CR2_4, ΔCR2_4 = _average_correlator(CorrR2_4)
CR3_4, ΔCR3_4 = _average_correlator(CorrR3_4)
CR4_4, ΔCR4_4 = _average_correlator(CorrR4_4)

CorrD_16  =  h5dset["T32_L24_h16_p1/p(0,0,1)/d/p_diag(0,0,1)/C_re"][]
Corrπ_16  =  h5dset["T32_L24_h16_p1/p(0,0,1)/pi/p_diag(0,0,1)/C_re"][]
Corrρ_16  =  h5dset["T32_L24_h16_p1/p(0,0,1)/rho_g33/p_diag(0,0,1)/C_re"][]
CorrT1_16 = -h5dset["T32_L24_h16_p1/p(0,0,1)/t1_g3/p_diag(0,0,1)/C_im"][]
CorrT2_16 =  h5dset["T32_L24_h16_p1/p(0,0,1)/t2_g3/p_diag(0,0,1)/C_im"][]
CorrR1_16 =  h5dset["T32_L24_h16_p1/p(0,0,1)/r1/p_diag(0,0,1)/C_re"][]
CorrR2_16 =  h5dset["T32_L24_h16_p1/p(0,0,1)/r2/p_diag(0,0,1)/C_re"][]
CorrR3_16 =  h5dset["T32_L24_h16_p1/p(0,0,1)/r3/p_diag(0,0,1)/C_re"][]
CorrR4_16 =  h5dset["T32_L24_h16_p1/p(0,0,1)/r4/p_diag(0,0,1)/C_re"][]

N, nhits, T = size(CorrD_16)
CD_16,  ΔCD_16  = _average_correlator(CorrD_16)
Cπ_16,  ΔCπ_16  = _average_correlator(Corrπ_16)
Cρ_16,  ΔCρ_16  = _average_correlator(Corrρ_16)
CT1_16, ΔCT1_16 = _average_correlator(CorrT1_16)
CT2_16, ΔCT2_16 = _average_correlator(CorrT2_16)
CR1_16, ΔCR1_16 = _average_correlator(CorrR1_16)
CR2_16, ΔCR2_16 = _average_correlator(CorrR2_16)
CR3_16, ΔCR3_16 = _average_correlator(CorrR3_16)
CR4_16, ΔCR4_16 = _average_correlator(CorrR4_16)

rel_difference(a,b)  = (2(a-b)/(a+b))
Δrel_difference(a,b,Δa,Δb) = (Δa+Δb)*(1/(a+b) + (a-b)/(a+b)^2)

diff  = rel_difference.(CD_16,CD_4)
Δdiff = Δrel_difference.(CD_16,CD_4,ΔCD_16,ΔCD_4)
scatter(diff,yerr=Δdiff,label="(D)")

diff  = rel_difference.(CT1_16,CT2_16)
Δdiff = Δrel_difference.(CT1_16,CT2_16,ΔCT1_16,ΔCT2_16)
scatter!(diff,yerr=Δdiff,label="(T) [T1 vs T2]")

diff  = rel_difference.(CR3_4,CR3_16)
Δdiff = Δrel_difference.(CR3_4,CR3_16,ΔCR3_4,ΔCR3_16)
scatter!(diff,yerr=Δdiff,label="(R3)")

diff  = rel_difference.(CR3_16,CR4_16)
Δdiff = Δrel_difference.(CR3_16,CR4_16,ΔCR3_16,ΔCR4_16)
scatter!(diff,yerr=Δdiff,label="(R) [R3 vs R4]")

diff  = rel_difference.(CR1_16,CR2_16)
Δdiff = Δrel_difference.(CR1_16,CR2_16,ΔCR1_16,ΔCR2_16)
scatter!(diff,yerr=Δdiff,label="(R) [R1 vs R2]")