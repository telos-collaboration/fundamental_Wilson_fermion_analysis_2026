using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
#gr(frame=:box, legend=:topright,legendfontsize=12)
plotlyjs(frame=:box)
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=16, legendfontsize=14,markersize=5,tickfontsize=12)
include("utils.jl")

hdf5file = "isospin1_L16_h2_p2.hdf5"
hdf5file = "isospin1_L24_h16_p23.hdf5"
hdf5file = "isospin1_L24_h16.hdf5"
h5dset = h5open(hdf5file)

p1  = "(0,0,1)"
p   = 1
ens = "E1"
T, L = h5dset["$ens/lattice"][1:2]

L3, L6 = L^3, L^6
title = L"%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1"
range = vcat(1:9,25:32)

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens;p)

N, nhits, T = size(CorrD1)

CD1, ΔCD1 = _average_correlator(CorrD1/L6)
CD2, ΔCD2 = _average_correlator(CorrD2/L6)
Cπ,  ΔCπ  = _average_correlator(Corrπ/L3)
Cρ,  ΔCρ  = _average_correlator(Corrρ/L3)
CT1, ΔCT1 = _average_correlator(CorrT1/L3)
CT2, ΔCT2 = _average_correlator(CorrT2/L3)
CR1, ΔCR1 = _average_correlator(CorrR1/L3)
CR2, ΔCR2 = _average_correlator(CorrR2/L3)
CR3, ΔCR3 = _average_correlator(CorrR3/L3)
CR4, ΔCR4 = _average_correlator(CorrR4/L3)

pltCross = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10,legend=:top)
pltMes   = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10)
pltPiPi  = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10, legend_columns=1, legend=:top)
pltR3R4  = plot(; title, xlabel=L"t",ylabel=L"C(t)",legend=:top)
pltPosNeg = plot(; title, xlabel=L"t",ylabel=L"C(t)",legend=:top,yscale=:log10)

scatter!(pltCross,CT1,yerr=ΔCT1,label=L"$-$Im(T1)",marker=:circle)
scatter!(pltCross,CT2,yerr=ΔCT2,label=L"$+$Im(T2)",marker=:star)
scatter!(pltMes,  Cπ ,yerr=ΔCπ ,label=L"\pi" ,marker=:rect,alpha=0.9)
scatter!(pltMes,  Cρ ,yerr=ΔCρ ,label=L"\rho",marker=:circ,alpha=0.9)


scatter!(pltPosNeg,CD1+CR1+CR2,yerr=ΔCR1+ΔCR2+ΔCD1,label="(D1+R1+R2)",marker=:rect)
scatter!(pltPosNeg,range,(CD2+CR3+CR4)[range],yerr=(ΔCD2+ΔCR3+ΔCR4)[range],label="(D2+R3+R4)",marker=:pent)
#scatter!(pltPosNeg,CD1+CR1+CR2-(CD2+CR3+CR4),yerr=ΔCR2,label="full",marker=:cross)
#scatter!(pltPosNeg,CD1+CR1+CR2+(CD2+CR3+CR4),yerr=ΔCR2,label="full (non relative sign)",marker=:cross)

cut = 6
t_R3R4 = vcat(1:cut+1,T-cut+1:T)
scatter!(pltPiPi,CD1, yerr=ΔCD1,label="D1",marker=:pent)
scatter!(pltPiPi,CD2, yerr=ΔCD2,label="D2",marker=:circ)
scatter!(pltPiPi,CR1,yerr=ΔCR1,label="R1/2",marker=:rect)
scatter!(pltPiPi,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3/4",marker=:pent)
#scatter!(pltPiPi,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)
#scatter!(pltPiPi,CR2,yerr=ΔCR2,label="R2",marker=:cross)

#scatter!(pltPiPi,CT1,yerr=ΔCT1,label=L"$-$Im(T1)",marker=:circle)
scatter!(pltPiPi,CT2,yerr=ΔCT2,label=L"$\pm$Im(T1/2)",marker=:star)


#display(pltR3R4)
#display(pltMes)
#display(pltCross)
display(pltPiPi)
#display(pltPosNeg)
savefig(pltPiPi,"all_diagrams.pdf")
savefig(pltPosNeg,"sum_of_each_sign.pdf")


save=false
if save
    isdir("plots") || mkpath("plots")
    savefig(pltPiPi,"plots/box_diagrams.pdf")
    savefig(pltR3R4,"plots/R3R4.pdf")
    savefig(pltMes,"plots/mesons.pdf")
    savefig(pltCross,"plots/triangles.pdf")
end