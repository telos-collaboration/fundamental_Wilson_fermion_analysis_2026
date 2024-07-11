using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
plotlyjs(frame=:box, legend=:topright,legendfontsize=12)
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=16, legendfontsize=14,markersize=5,tickfontsize=12)
include("utils.jl")

# Lt32Ls24beta6.9m1-0.92m2-0.92
hdf5file = "isospin1_h4_raw.hdf5"
hdf5file = "isospin1_h4_L16_raw.hdf5"
h5dset = h5open(hdf5file)

p1  = "(1,1,0)"
p1  = "(0,0,1)"
ens = "E1"
T, L = h5dset["$ens/lattice"][1:2]
L3 = L^3
cut = 4
title = L"%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1"

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens)
Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp110(h5dset,ens)

N, nhits, T = size(CorrD1)

CD1, ΔCD1 = _average_correlator(CorrD1)
CD2, ΔCD2 = _average_correlator(CorrD2)
Cπ,  ΔCπ  = _average_correlator(Corrπ)
Cρ,  ΔCρ  = _average_correlator(Corrρ)
CT1, ΔCT1 = _average_correlator(CorrT1)
CT2, ΔCT2 = _average_correlator(CorrT2)
CR1, ΔCR1 = _average_correlator(CorrR1)
CR2, ΔCR2 = _average_correlator(CorrR2)
CR3, ΔCR3 = _average_correlator(CorrR3)
CR4, ΔCR4 = _average_correlator(CorrR4)

pltCross = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10,legend=:top)
pltMes   = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10)
pltPiPi  = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10, legend_columns=2, legend=:top)
pltR3R4  = plot(; title, xlabel=L"t",ylabel=L"C(t)",legend=:top)

scatter!(pltCross,CT1,yerr=ΔCT1,label=L"$-$Im(T1)",marker=:circle)
scatter!(pltCross,CT2,yerr=ΔCT2,label=L"$+$Im(T2)",marker=:star)

scatter!(pltMes,Cπ, yerr=ΔCπ, label=L"\pi" ,marker=:rect,alpha=0.9)
scatter!(pltMes,Cρ, yerr=ΔCρ, label=L"\rho",marker=:circ,alpha=0.9)

t_R3R4 = vcat(1:cut+1,T-cut+1:T)
scatter!(pltPiPi,CR1,yerr=ΔCR1,label="R1",marker=:rect)
scatter!(pltPiPi,CR2,yerr=ΔCR2,label="R2",marker=:cross)
scatter!(pltPiPi,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3",marker=:pent)
scatter!(pltPiPi,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)
scatter!(pltPiPi,CD1, yerr=ΔCD1,label="D1",marker=:pent)
scatter!(pltPiPi,CD2, yerr=ΔCD2,label="D2",marker=:circ)

t_R3R4 = vcat(cut:T-cut+2)
scatter!(pltR3R4,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3",marker=:pent)
scatter!(pltR3R4,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)

display(pltR3R4)
display(pltMes)
display(pltCross)
display(pltPiPi)

save=false
if save
    isdir("plots") || mkpath("plots")
    savefig(pltPiPi,"plots/box_diagrams.pdf")
    savefig(pltR3R4,"plots/R3R4.pdf")
    savefig(pltMes,"plots/mesons.pdf")
    savefig(pltCross,"plots/triangles.pdf")
end