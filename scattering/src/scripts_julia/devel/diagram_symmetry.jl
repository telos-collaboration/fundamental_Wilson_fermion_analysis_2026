using Pkg; Pkg.activate(".")
using HDF5
using Statistics
using Plots
using LatticeUtils
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.6,ms=6,titlefontsize=11)
SHAPES = [:circle, :pentagon, :diamond, :hexagon, :rect, :octagon, :star4, :dtriangle,:star5, :ltriangle, :star6, :rtriangle, :star7, :utriangle, :star8]

function avg_corr(corr;conf_dim,avg_dims=conf_dim)
    Nconf = size(corr,conf_dim)
    C  = dropdims(mean(corr,dims=avg_dims),dims=avg_dims)
    ΔC = dropdims(std(corr,dims=avg_dims),dims=avg_dims)/sqrt(Nconf)
    return C, ΔC
end
function symmetry_ratio(Corr;t_dim)
    T      = size(Corr,t_dim)
    left   = 2:T÷2
    right  = T÷2+2:T
    Cleft  = selectdim(Corr,t_dim,left)
    Cright = selectdim(Corr,t_dim,reverse(right))
    ratio  = @. Cleft/Cright
    return left, ratio
end
# This allocates more than LatticeUtils.correlatorfolding() but is 
# much shorter and a bit faster in exchange.
function correlator_folding_alt(Corr;t_dim,sign)
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
    t, r  = symmetry_ratio(Corr;t_dim)
    R, ΔR = avg_corr(r;conf_dim,avg_dims)
    t_dim = t_dim - count(d -> isless(d,t_dim), avg_dims)
    return t, R, ΔR, t_dim
end
pseudolog10(x,C=1E+4) = sign(x)*log10(abs(C*x)+1)
function plot_3x3_matrix_elements!(plt,C;conf_dim,avg_dims=conf_dim)
    Cpl, ΔCpl = avg_corr(pseudolog10.(C);conf_dim,avg_dims)
    si = 1
    for i in 1:3, j in i:3
        Cij  = Cpl[i,j,:]
        ΔCij = ΔCpl[i,j,:]
        Cji  = Cpl[i,j,:]
        ΔCji = ΔCpl[i,j,:]
        # Plot only real or imaginary part, depending on which one is non-zero 
        plot_real = sum(abs∘real,Cij) > sum(abs∘imag,Cij)
        if plot_real
            plot!(plt,real.(Cij),yerr=real.(ΔCij),marker=SHAPES[si],label="real part of ($i,$j)")
            si += 1
        else
            plot!(plt,imag.(Cij),yerr=imag.(ΔCij), marker=SHAPES[si],label="imag part of ($i,$j)")
            plot!(plt,imag.(Cji),yerr=imag.(ΔCji), marker=SHAPES[si+1],label="imag part of ($j,$i)")
            si += 2
        end
    end
    return plt
end  
function fold_3x3_correlator(c)
    c_folded = similar(c)
    c_folded[1:2,1:2,:,:] .= correlator_folding(c[1:2,1:2,:,:];t_dim=4,sign=+1)
    c_folded[1,3,:,:] .= correlator_folding(c[1,3,:,:];t_dim=2,sign=-1)
    c_folded[3,1,:,:] .= correlator_folding(c[3,1,:,:];t_dim=2,sign=-1)
    c_folded[2,3,:,:] .= correlator_folding(c[2,3,:,:];t_dim=2,sign=-1)
    c_folded[3,2,:,:] .= correlator_folding(c[3,2,:,:];t_dim=2,sign=-1)
    c_folded[3,3,:,:] .= correlator_folding(c[3,3,:,:];t_dim=2,sign=+1)
    return c_folded
end

h5file1  = "data_assets/data/isospin1_eigenvalues_t0_6_deriv_false.hdf5"
fid1     = h5open(h5file1)
Corr1    = read(fid1[(first(keys(fid1)))],"p(0,0,1)/Corr3x3")
t_dim    = 4
conf_dim = 3
avg_dims = 3

t1, r1, Δr1, t_dim1 = symmetry_ratio_avg(Corr1;t_dim,conf_dim,avg_dims)
Corr1_sym           = fold_3x3_correlator(Corr1)
ts, rs, Δrs, t_dim1 = symmetry_ratio_avg(Corr1_sym;t_dim,conf_dim,avg_dims)
plt1 = plot_3x3_matrix_elements!(plot(),Corr1;conf_dim,avg_dims)
plt2 = plot_3x3_matrix_elements!(plot(),Corr1_sym;conf_dim,avg_dims)
