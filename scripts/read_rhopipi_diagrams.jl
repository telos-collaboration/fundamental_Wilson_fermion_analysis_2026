function _average_correlator(Corr)
    N, nhits, T = size(Corr)
    CorrAvg = dropdims(mean(Corr,dims=2),dims=2)
    C  = dropdims(mean(CorrAvg,dims=1),dims=1)
    ΔC = dropdims(std(CorrAvg,dims=1),dims=1)./sqrt(N)
    return C, ΔC
end
function D1(h5dset,p,ens)
    πp =  h5dset["$ens/p$p/pi/p_diag$p/C_re"][]
    π0 =  h5dset["$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re"][]
    return πp .* π0
end
function correlatorsp001(h5dset,ens;p=1)
    p1 = "(0,0,$p)"
    CorrD1_old = D1(h5dset,p1,ens)
    CorrD1 =  h5dset["$ens/p$p1/d/p_diag$p1/C_re"][]
    CorrD2 =  h5dset["$ens/p$p1/d/p_diag(0,0,0)/C_re"][]
    Corrπ  =  h5dset["$ens/p$p1/pi/p_diag$p1/C_re"][]
    try
        Corrρ  =  h5dset["$ens/p$p1/rho_g33/p_diag$p1/C_re"][]
    catch
        Corrρ  =  h5dset["$ens/p$p1/rho_g3_g3/p_diag$p1/C_re"][]
    end
    CorrR1 =  h5dset["$ens/p$p1/r1/p_diag$p1/C_re"][]
    CorrR2 =  h5dset["$ens/p$p1/r2/p_diag$p1/C_re"][]
    CorrR3 =  h5dset["$ens/p$p1/r3/p_diag$p1/C_re"][]
    CorrR4 =  h5dset["$ens/p$p1/r4/p_diag$p1/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old
end
function correlatorsp110(h5dset,ens;p=1)
    p1 = "($p,$p,0)"
    CorrD1_old = D1(h5dset,p1,ens)
    CorrD1 =  h5dset["$ens/p$p1/d/p_diag$p1/C_re"][]
    CorrD2 =  h5dset["$ens/p$p1/d/p_diag(0,0,0)/C_re"][]
    Corrπ  = h5dset["$ens/p$p1/pi/p_diag$p1/C_re"][]
    # THE FOLLOWING DOES NOT GIVE A SIGNAL. IS SOMETHING WRONG WITH THE PARSING?
    Corrρ  =  (h5dset["$ens/p$p1/rho_g11/p_diag$p1/C_re"][] + h5dset["$ens/p$p1/rho_g22/p_diag$p1/C_re"][] + h5dset["$ens/p$p1/rho_g21/p_diag$p1/C_re"][] + h5dset["$ens/p$p1/rho_g12/p_diag$p1/C_re"][])/2
    CorrT1 =  (h5dset["$ens/p$p1/t1_g1/p_diag$p1/C_im"][]   + h5dset["$ens/p$p1/t1_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrT2 =  (h5dset["$ens/p$p1/t2_g1/p_diag$p1/C_im"][]   + h5dset["$ens/p$p1/t2_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrR1 =   h5dset["$ens/p$p1/r1/p_diag$p1/C_re"][]
    CorrR2 =   h5dset["$ens/p$p1/r2/p_diag$p1/C_re"][]
    CorrR3 =   h5dset["$ens/p$p1/r3/p_diag$p1/C_re"][]
    CorrR4 =   h5dset["$ens/p$p1/r4/p_diag$p1/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old
end
function correlatorsp110_alt(h5dset,ens;p=1)
    p1 = "($p,$p,0)"
    CorrD1_old = D1(h5dset,p1,ens)
    CorrD1 =  h5dset["$ens/p$p1/d/p_diag$p1/C_re"][]
    CorrD2 =  h5dset["$ens/p$p1/d/p_diag(0,0,0)/C_re"][]
    Corrπ  = h5dset["$ens/p$p1/pi/p_diag$p1/C_re"][]
    # THE FOLLOWING IS ONLY TAKEN ONE POLARISATION INTO ACCOUNT - IT IS NOT FULLY SELFCONSISTENT
    Corrρ  =   h5dset["$ens/p$p1/rho_g11/p_diag$p1/C_re"][]
    CorrT1 =  (h5dset["$ens/p$p1/t1_g1/p_diag$p1/C_im"][] + h5dset["$ens/p$p1/t1_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrT2 =  (h5dset["$ens/p$p1/t2_g1/p_diag$p1/C_im"][] + h5dset["$ens/p$p1/t2_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrR1 =   h5dset["$ens/p$p1/r1/p_diag$p1/C_re"][]
    CorrR2 =   h5dset["$ens/p$p1/r2/p_diag$p1/C_re"][]
    CorrR3 =   h5dset["$ens/p$p1/r3/p_diag$p1/C_re"][]
    CorrR4 =   h5dset["$ens/p$p1/r4/p_diag$p1/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old
end
function correlatorsp000(h5dset,ens;p=1)
    Corrπ  = h5dset["$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re"][]
    Corrρ  = h5dset["$ens/p(0,0,0)/rho_g1/p_diag(0,0,0)/C_re"][]
    return Corrπ, Corrρ
end
function read_hdf5_file(file,ens,p1,p)
    h5dset = h5open(file)
    T, L = h5dset["$ens/lattice"][1:2]
    if p1 == "(0,0,0)"
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp000(h5dset,ens)
    elseif p1 ∈ [ "(0,0,1)" , "(0,0,2)" , "(0,0,3)"]
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp001(h5dset,ens;p)
    elseif p1 ∈ [ "(1,1,0)" , "(2,2,0)" , "(3,3,0)"]
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old = correlatorsp110(h5dset,ens;p)
    end
    return T, L, Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2, CorrD1_old
end