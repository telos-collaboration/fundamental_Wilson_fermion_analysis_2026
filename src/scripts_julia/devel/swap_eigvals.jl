using Pkg; Pkg.activate("src/src_jl")
using HDF5
using LatticeUtils
using Plots
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

fid = h5open("data_assets/isospin1_eigenvalues_evp_gevp.hdf5")
r  = keys(fid)[3]
id = "gevp_t0_8"
p  = "p(0,0,1)"
irrep = "A1"

T   = read(fid[r],"lattice")[1]
ev  = read(fid[r]["$p/$irrep/$id"],"eigvals_3x3")
Δev = read(fid[r]["$p/$irrep/$id"],"Delta_eigvals_3x3")

function swap_eigvals(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = t:T-t+2
    @. new[op1,r1] = old[op2,r1]
    @. new[op2,r1] = old[op1,r1]
    return new
end
function swap_eigvals_start(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = 1:t-1
    r2 = T-t+2:T
    @. new[op1,r1] = old[op2,r1]
    @. new[op2,r1] = old[op1,r1]
    @. new[op1,r2] = old[op2,r2]
    @. new[op2,r2] = old[op1,r2]
    return new
end

op1,op2,t = 1,3,8
ev = swap_eigvals_start(ev,op1,op2,t)
Δev = swap_eigvals_start(Δev,op1,op2,t)
println("$r,\"$p\",$irrep,$id,$op1,$op2,-$t")

op1,op2,t = 2,3,13
ev = swap_eigvals(ev,op1,op2,t)
Δev = swap_eigvals(Δev,op1,op2,t)
println("$r,\"$p\",$irrep,$id,$op1,$op2,$t")

t = filter(!isequal(T÷2+1),1:T)
f = abs
plt = plot(yscale=:log10,legend=:top)
plot_correlator!(plt,1:T,f.(ev[1,t]),Δev[1,t],markersize=3,markershape=:rect,label="eigval #1 (3x3)")
plot_correlator!(plt,1:T,f.(ev[2,t]),Δev[2,t],markersize=3,markershape=:rect,label="eigval #2 (3x3)")    
plot_correlator!(plt,1:T,f.(ev[3,t]),Δev[3,t],markersize=3,markershape=:rect,label="eigval #3 (3x3)")    
