using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
pgfplotsx(frame=:box, legend=:topright,labelfontsize=16, titlefontsize=16, legendfontsize=14,markersize=5,tickfontsize=12)
gr(frame=:box, legend=:topright,legendfontsize=12)
plotlyjs(frame=:box, legend=:topright,legendfontsize=12)
include("utils.jl")

# Lt32Ls24beta6.9m1-0.92m2-0.92
file1 = "isospin1.hdf5"
file2 = "isospin1.hdf5"

p1A = "(0,0,2)"
p1B = "(0,0,3)"
ens = "T32_L24_h16_p23"
cut = 4

# T, L, Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old
ret1 = read_hdf5_file(file1,ens,p1A,2)
ret2 = read_hdf5_file(file2,ens,p1B,3)

L1 = ret1[2]
L2 = ret2[2]

plt1 = plot(;xlabel="t",ylabel="C(t)",yscale=:log10,legend=:top)
plt3 = plot(;xlabel="t",ylabel="D1/D2 - ratio",legend=:top)
ratio(C1,C2) = C1 / C2
Δratio(C1,C2,ΔC1,ΔC2) = ΔC1 / C2 + ΔC2 * C1 / C2 ^ 2
N1 = L1^(6)
N2 = L2^(6)
s  = (N1/N2)

for ind_corr in [11,12]
    Corr1 = ret1[ind_corr] 
    Corr2 = ret2[ind_corr] 
    C1,  ΔC1  = _average_correlator(Corr1)
    C2,  ΔC2  = _average_correlator(Corr2)
    scatter!(plt1,C1/N1, yerr=ΔC1/N1, label="π: L=$L1, p=$p1A" ,marker=:rect,alpha=0.9)
    scatter!(plt1,C2/N2, yerr=ΔC2/N2, label="π: L=$L2, p=$p1B" ,marker=:rect,alpha=0.9)
end

# compare old and new D1: They differ by factor of s0
CorrD1a = ret1[11]
CorrD1b = ret1[13] 
CD1a, ΔCD1a = _average_correlator(CorrD1a)
CD1b, ΔCD1b = _average_correlator(CorrD1b)
s0 = 1/48
scatter!(plt3,ratio.(CD1a*s0,CD1b), yerr=Δratio.(CD1a*s0,CD1b,ΔCD1a*s0,ΔCD1b), label="" ,marker=:rect,alpha=0.9)

plot!(plt2,ylims=(0,2))
plot!(plt3,ylims=(0,2))
display(plt1)
display(plt3)
#savefig(plt1,"pions_different_L.pdf")
#savefig(plt2,"pions_different_L_ratio.pdf")