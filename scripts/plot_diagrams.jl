using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
#gr(frame=:box, legend=:topright,legendfontsize=12)
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=11, legendfontsize=8,markersize=5,tickfontsize=12)
plotlyjs(frame=:box)
include("utils.jl")

hdf5file = "isospin1_L16_PA_h4.hdf5"
hdf5file = "isospin1_L24_h16_p23.hdf5"
hdf5file = "isospin1_L24_h16.hdf5"
hdf5file = "isospin1_L16_h4.hdf5"
hdf5file = "isospin1_L24_h16.hdf5"
h5dset = h5open(hdf5file)

p1  = "(0,0,1)"
p   = 1
ens = "E1"
T, L = h5dset["$ens/lattice"][1:2]

L3, L6 = L^3, L^6
title = L"%$T \times %$L^3, \beta=6.9, m_0^f=-0.92, \mathbf p = %$p1"
range = vcat(1:9,25:32)
range = vcat(1:32)

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

function pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
    L3, L6 = L^3, L^6
    Corr2π = CorrD1/L6+2CorrR1/L3+CorrD2/L6+2CorrR3/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2)
    N, nhits, T = size(Corr2π)
    corr = zeros(2,2,N,nhits,T)
    corr[1,1,:,:,:] = Corrρ
    corr[1,2,:,:,:] = CorrT1
    corr[2,1,:,:,:] = CorrT2
    corr[2,2,:,:,:] = Corr2π
    return corr
end
Corr2π = pipi_correlator(CorrD1,CorrR1,CorrD2,CorrR3,L)
pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2)
C2π, ΔC2π  = _average_correlator(Corr2π)

pltCross = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10)
pltMes   = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10)
pltPiPi  = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10)
pltR3R4  = plot(; title, xlabel=L"t",ylabel=L"C(t)")

t_R3R4 = 1:32
scatter!(pltCross,CT1,yerr=ΔCT1,label=L"$-$Im(T1)",marker=:circle)
scatter!(pltCross,CT2,yerr=ΔCT2,label=L"$+$Im(T2)",marker=:star)
scatter!(pltMes,  Cπ ,yerr=ΔCπ ,label=L"\pi" ,marker=:rect,alpha=0.9)
scatter!(pltMes,  Cρ ,yerr=ΔCρ ,label=L"\rho",marker=:circ,alpha=0.9)
scatter!(pltPiPi,C2π,yerr=ΔC2π,label="full (non relative sign)",marker=:cross)
scatter!(pltR3R4,t_R3R4,CR3[t_R3R4],yerr=ΔCR3[t_R3R4],label="R3",marker=:pent)
scatter!(pltR3R4,t_R3R4,CR4[t_R3R4],yerr=ΔCR4[t_R3R4],label="R4",marker=:star)

save=false
if save
    isdir("plots") || mkpath("plots")
    savefig(pltPiPi,"plots/full_correlator.pdf")
    savefig(pltR3R4,"plots/R3R4.pdf")
    savefig(pltMes,"plots/mesons.pdf")
    savefig(pltCross,"plots/triangles.pdf")
else
    display(pltCross)
    display(pltMes)
    display(pltPiPi)
    display(pltR3R4)    
end