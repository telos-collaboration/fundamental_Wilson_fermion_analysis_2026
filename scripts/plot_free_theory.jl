using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using HDF5
using Statistics
using LaTeXStrings
include("utils.jl")
gr(frame=:box,legendfontsize=12)

hdf5file = expanduser("~/Downloads/free_theory_results_v2.hdf5")
h5dset = h5open(hdf5file)

p   = 1
ens = "T16L8"
T,L = h5dset["$ens/lattice"][1:2]
m = -0.6

Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp001(h5dset,ens;p)
#average over  stochastic samples
CD2 = dropdims(mean(CorrD2,dims=2),dims=(1,2))
CD1 = dropdims(mean(CorrD1,dims=2),dims=(1,2))
Cπ = dropdims(mean(Corrπ,dims=2),dims=(1,2))
Cρ = dropdims(mean(Corrρ,dims=2),dims=(1,2))

m_p(m,p)       = m + sum( 1.0 .- cos.(p))
N_p(m,p)       = m_p(m,p)^2 + sum(sin.(p).^2)
ap_dot_aq(p,q) = sum(sin.(p).*sin.(q))

function analytic_free_pion(T,L,m,P)
    mom_t = (0:T-1).*(2π/T)
    mom_s = (0:L-1).*(2π/L)
    C = zeros(ComplexF64,T)
    for t in 0:T-1, q0 in mom_t, p0 in mom_t, p1 in mom_s, p2 in mom_s, p3 in mom_s
        q = [q0, p1 + P[2], p2 + P[3], p3 + P[4]]
        p = [p0, p1       , p2       , p3       ]
        C[t+1] += exp(+im*(p0-q0)*t)*(m_p(m,p)*m_p(m,q) + ap_dot_aq(p,q))/N_p(m,p)/N_p(m,q)
    end
    return @. 4C /(L^3*T^2)
end
P = [0,0,0,1*2π/L]
C = analytic_free_pion(T,L,m,P)

plt = plot(legend=:top)
scatter!(plt,abs.(Cπ),yscale=:log10,label="π (measured)")
scatter!(plt,real.(C),yscale=:log10,label="π (analytic)")
