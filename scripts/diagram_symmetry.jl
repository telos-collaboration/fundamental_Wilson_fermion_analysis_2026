using Pkg; Pkg.activate(".")
using HDF5
using Statistics
using Plots
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function symmetry_ratio(Corr;t_dim)
    T      = size(Corr,t_dim)
    left   = 2:T÷2
    right  = T÷2+2:T
    Cleft  = selectdim(Corr,t_dim,left)
    Cright = selectdim(Corr,t_dim,reverse(right))
    ratio  = @. Cleft/Cright
    return left, ratio
end
function symmetrise_correlator(Corr;t_dim,sign)
    T      = size(Corr,t_dim)
    left   = 2:T÷2
    right  = T÷2+2:T
    Cleft  = selectdim(Corr,t_dim,left)
    Cright = selectdim(Corr,t_dim,reverse(right))
    C0     = selectdim(Corr,t_dim,1:1)
    Cmid   = selectdim(Corr,t_dim,T÷2+1:T÷2+1)
    Cleft_sym   = @. (Cleft + sign*Cright)/2 
    Cright_sym  = @. (sign*Cleft + Cright)/2 
    C_symmetric = cat(C0,Cleft_sym,Cmid,reverse(Cright_sym,dims=t_dim),dims=t_dim)
    return C_symmetric
end
function symmetry_ratio_avg(Corr;t_dim,conf_dim,avg_dims)
    Nconf = size(Corr,conf_dim)
    t, r  = symmetry_ratio(Corr;t_dim)
    R     = dropdims(mean(r,dims=avg_dims),dims=avg_dims)
    ΔR    = dropdims(std(r,dims=avg_dims),dims=avg_dims)/sqrt(Nconf)
    t_dim = t_dim - count(d -> isless(d,t_dim), avg_dims)
    return t, R, ΔR, t_dim
end


h5file1 = "output/data/isospin1_eigenvalues_t0_6_deriv_false.hdf5"
h5file2 = "output/data/isospin1_merged.hdf5"
fid1  = h5open(h5file1)
fid2  = h5open(h5file2)
Corr1 = read(fid1[(first(keys(fid1)))],"p(0,0,1)/Corr3x3")
Corr2 = read(fid2[(first(keys(fid2)))],"p(0,0,0)/pi/p_diag(0,0,0)/C_re")

t1, r1, Δr1, t_dim1 = symmetry_ratio_avg(Corr1;t_dim=4,conf_dim=3,avg_dims=(3))
t2, r2, Δr2, t_dim2 = symmetry_ratio_avg(Corr2;t_dim=3,conf_dim=1,avg_dims=(1,2))


Corr1_sym = symmetrise_correlator(Corr1;t_dim=4,sign=+1)

Corr1_sym

t1, r1, Δr1, t_dim1 = symmetry_ratio_avg(Corr1_sym;t_dim=4,conf_dim=3,avg_dims=(3))
scatter(t1,real.(r1[1,1,:]),yerr=real.(Δr1[1,1,:]))


