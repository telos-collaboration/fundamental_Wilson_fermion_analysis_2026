using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
plotlyjs(frame=:box, legend=:topright,legendfontsize=12)
pgfplotsx(frame=:box, legend=:topright,legendfontsize=12,markersize=5,tickfontsize=12)
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

p1 = "(0,0,1)"
CorrD  =  h5dset["E1/p$p1/d/p_diag$p1/C_re"][]
Corrπ  =  h5dset["E1/p$p1/pi/p_diag$p1/C_re"][]
Corrρ  =  h5dset["E1/p$p1/rho_g33/p_diag$p1/C_re"][]
CorrT1 = -h5dset["E1/p$p1/t1_g3/p_diag$p1/C_im"][]
CorrT2 =  h5dset["E1/p$p1/t2_g3/p_diag$p1/C_im"][]
CorrR1 =  h5dset["E1/p$p1/r1/p_diag$p1/C_re"][]
CorrR2 =  h5dset["E1/p$p1/r2/p_diag$p1/C_re"][]
CorrR3 =  h5dset["E1/p$p1/r3/p_diag$p1/C_re"][]
CorrR4 =  h5dset["E1/p$p1/r4/p_diag$p1/C_re"][]

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

pltCross = plot(yscale=:log10,legend=:top)
scatter!(pltCross,CT1,yerr=ΔCT1,label=L"$-$Im(T1)",marker=:circle)
scatter!(pltCross,CT2,yerr=ΔCT2,label=L"$+$Im(T2)",marker=:star)

pltMes  = plot(yscale=:log10)
scatter!(pltMes,Cπ, yerr=ΔCπ, label=L"$\pi~  \bf p =%$p$",alpha=0.5)
scatter!(pltMes,Cρ, yerr=ΔCρ, label=L"$\rho~ \bf p =%$p$",alpha=0.5)

t_R3R4 = vcat(1:7,27:32)
pltPiPi = plot(yscale=:log10, legend_columns=2, legend=:top)
scatter!(pltPiPi,CR1,yerr=ΔCR1,label="R1",marker=:rect)
scatter!(pltPiPi,CR2,yerr=ΔCR2,label="R2",marker=:cross)
scatter!(pltPiPi,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3",marker=:pent)
scatter!(pltPiPi,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)
scatter!(pltPiPi,CD, yerr=ΔCD ,label="D2",marker=:circle)

pltR1R2 = plot(legend=:top)
t_R3R4 = vcat(6:28)
scatter!(pltR1R2,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3",marker=:pent)
scatter!(pltR1R2,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)


pltPiPi
pltR1R2
pltMes
pltCross