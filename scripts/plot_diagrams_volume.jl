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
file1 = "isospin1_h4_raw.hdf5"
file2 = "isospin1_h4_L16_raw.hdf5"

p1  = "(0,0,0)"
ens = "E1"

cut = 4

# T, L, Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
ret1 = read_hdf5_file(file1,ens,p1)
ret2 = read_hdf5_file(file2,ens,p1)

L1 = ret1[2]
L2 = ret2[2]
Corr1 = ret1[3] 
Corr2 = ret2[3] 

title = L"\beta=6.9, m_0^f=-0.92, \mathbf p = %$p1"
C1,  ΔC1  = _average_correlator(Corr1)
C2,  ΔC2  = _average_correlator(Corr2)

plt1 = plot(; title, xlabel=L"t",ylabel=L"C(t)",yscale=:log10,legend=:top)
plt2 = plot(; title, xlabel=L"t",ylabel=L"$C(t)$ ratio",legend=:top)

scatter!(plt1,C1, yerr=ΔC1, label=L"\pi (L=%$L1)" ,marker=:rect,alpha=0.9)
scatter!(plt1,C2, yerr=ΔC2, label=L"\pi (L=%$L2)" ,marker=:rect,alpha=0.9)

scatter!(plt2,C1 ./ C2, yerr=ΔC1 ./ C2 .+ ΔC2 .* C1 ./ C2 .^ 2 , label="" ,marker=:rect,alpha=0.9)

display(plt2)