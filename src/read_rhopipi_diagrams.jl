function correlators_xyz(h5dset,ens;p::Vector{Int})
    p_e = "($(p[1]),$(p[2]),$(p[3]))"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]

    # The following entries are always the same
    Corrπ  = read(group,"pi/$p_i/C_re")
    CorrD1 = read(group,"d/$p_i/C_re")
    CorrD2 = read(group,"d/$p_0/C_re")
    CorrR1 = read(group,"r1/$p_i/C_re")
    CorrR2 = read(group,"r2/$p_i/C_re")
    CorrR3 = read(group,"r3/$p_i/C_re")
    CorrR4 = read(group,"r4/$p_i/C_re")

    # The operators that result from inclusion of vector-meson operators depend on the external momentum which determines the required polarizations
    normρ  = dot(p,p)
    normT  = norm(p)
    Corrρ  = zero(Corrπ)
    CorrT1 = zero(Corrπ)
    CorrT2 = zero(Corrπ)

    for i in 1:3, j in 1:3
        sign  = ifelse(i == 2 , -1 , +1)
        try
            Corrρ = Corrρ + p[i]*p[j]*sign*read(group,"rho_g$(i)_g$(j)/$p_i/C_re")/normρ
        catch
            Corrρ = Corrρ + p[i]*p[j]*sign*read(group,"rho_g$(i)$(j)/$p_i/C_re")/normρ
        end
    end
    for i in 1:3
        CorrT1 = CorrT1 + p[i]*read(group,"t1_g$(i)/$p_i/C_im")/normT
        CorrT2 = CorrT2 + p[i]*read(group,"t2_g$(i)/$p_i/C_im")/normT
    end

    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end

function correlatorsp000(h5dset,ens)
    Corrπ  = read(h5dset,"$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re")
    Corrρ  = read(h5dset,"$ens/p(0,0,0)/rho_g1/p_diag(0,0,0)/C_re")
    return Corrπ, Corrρ
end