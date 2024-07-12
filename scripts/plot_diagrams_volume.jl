using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=16, legendfontsize=14,markersize=5,tickfontsize=12)
gr(frame=:box, legend=:topright,legendfontsize=12)
include("utils.jl")

# Lt32Ls24beta6.9m1-0.92m2-0.92
file1 = "isospin1_L16_h2_p2.hdf5"
file2 = "isospin1_L24_h16_p23.hdf5"

p1A = "(0,0,2)"
p1B = "(0,0,3)"
ens = "E1"
cut = 4

# T, L, Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
ret1 = read_hdf5_file(file1,ens,p1A,2)
ret2 = read_hdf5_file(file2,ens,p1B,3)

L1 = ret1[2]
L2 = ret2[2]
ind_corr = 1
Corr1 = ret1[ind_corr] 
Corr2 = ret2[ind_corr] 

C1,  ΔC1  = _average_correlator(Corr1)
C2,  ΔC2  = _average_correlator(Corr2)

plt1 = plot(;xlabel=L"t",ylabel=L"C(t)",yscale=:log10,legend=:top)
plt2 = plot(;xlabel=L"t",ylabel=L"$C(t)$ ratio",legend=:top)

ratio(C1,C2) = C1 / C2
Δratio(C1,C2,ΔC1,ΔC2) = ΔC1 / C2 + ΔC2 * C1 / C2 ^ 2
N1 = L1^(6)
N2 = L2^(6)
s  = (N1/N2)

scatter!(plt1,C1/N1, yerr=ΔC1/N1, label=L"\pi (L=%$L1)" ,marker=:rect,alpha=0.9)
scatter!(plt1,C2/N2, yerr=ΔC2/N2, label=L"\pi (L=%$L2)" ,marker=:rect,alpha=0.9)
scatter!(plt2,ratio.(C1,C2)/s, yerr=Δratio.(C1,C2,ΔC1,ΔC2)/s, label="" ,marker=:rect,alpha=0.9)

plot!(plt2,ylims=(-0,5))

display(plt2)
#display(plt1)