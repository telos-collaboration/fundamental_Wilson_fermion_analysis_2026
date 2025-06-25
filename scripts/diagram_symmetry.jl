using Pkg; Pkg.activate(".")
using HDF5
using Statistics
using Plots
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

function fold_ratio(Corr,T,t_dim)
    left   = 2:T÷2
    right  = T÷2+2:T
    Cleft  = selectdim(Corr,t_dim,left)
    Cright = selectdim(Corr,t_dim,reverse(right))
    ratio  = @. abs(Cleft/Cright)
    return left, ratio
end


h5file1 = "output/data/isospin1_eigenvalues_t0_6_deriv_false.hdf5"
h5file2 = "output/data/isospin1_merged.hdf5"
fid1  = h5open(h5file1)
fid2  = h5open(h5file2)
Corr1 = read(fid1[(first(keys(fid1)))],"p(0,0,1)/Corr3x3")
Corr2 = read(fid2[(first(keys(fid2)))],"p(0,0,1)/d/p_diag(0,0,1)/C_re")
T1 = read(fid1[(first(keys(fid1)))],"lattice")[1]
T2 = read(fid2[(first(keys(fid2)))],"lattice")[1]

t_dim = 4
c_dim = 3
avg_dims = (3)
plot_dims = (1,2)

left, ratio = fold_ratio(Corr1,T,t_dim)
Nconf = size(ratio,c_dim)
C     = dropdims(mean(ratio,dims=avg_dims),dims=avg_dims)
ΔC    = dropdims(std(ratio,dims=avg_dims),dims=avg_dims)/sqrt(Nconf)
t_dim = t_dim - count(d -> isless(t_dim,d), avg_dims)

plt = plot()
for i in 1:3, j in 1:3
    scatter!(plt,left,C[i,j,:],yerr = ΔC[i,j,:],ylims=(0.9,1.1),label="($i,$j)")
end
plt