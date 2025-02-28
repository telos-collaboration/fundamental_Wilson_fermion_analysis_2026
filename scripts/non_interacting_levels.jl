using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using DelimitedFiles
using LatticeUtils

mπ,Δmπ = 0.38649,0.00051

for L in [12,14,16,24]  
    for (px,py,pz) in [(0,0,1),(0,1,1),(1,1,1),(0,0,2),(0,1,2),(1,1,2)]
        E, ΔE = non_interacting_energy_2P_lattice(mπ,Δmπ,px,py,pz,L)
        if E < 1.05
            println("L=$L: p=($px,$py,$pz): E_ππ(n.i.)_cont = $(errorstring(E, ΔE))")
        end
    end
end