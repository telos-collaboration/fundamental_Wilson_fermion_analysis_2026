function correlatorsp001(h5dset,ens;p=1)
    p_e = "(0,0,$p)"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]
    CorrD1 = read(group,"d/$p_i/C_re")
    CorrD2 = read(group,"d/$p_0/C_re")
    Corrπ  = read(group,"pi/$p_i/C_re")
    Corrρ = try
        read(group,"rho_g3_g3/$p_i/C_re")
    catch
        read(group,"rho_g33/$p_i/C_re")
    end
    CorrT1 = read(group,"t1_g3/$p_i/C_im")
    CorrT2 = read(group,"t2_g3/$p_i/C_im")
    CorrR1 = read(group,"r1/$p_i/C_re")
    CorrR2 = read(group,"r2/$p_i/C_re")
    CorrR3 = read(group,"r3/$p_i/C_re")
    CorrR4 = read(group,"r4/$p_i/C_re")
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function correlatorsp110(h5dset,ens;p=1)
    p_e = "($p,$p,0)"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]
    CorrD1 =  read(group,"d/$p_i/C_re")
    CorrD2 =  read(group,"d/$p_0/C_re")
    Corrπ  =  read(group,"pi/$p_i/C_re")
    Corrρ = try
         (read(group,"rho_g1_g1/$p_i/C_re") + read(group,"rho_g1_g2/$p_i/C_re") - read(group,"rho_g2_g1/$p_i/C_re") - read(group,"rho_g2_g2/$p_i/C_re"))/2
    catch
        (read(group,"rho_g11/$p_i/C_re") + read(group,"rho_g12/$p_i/C_re") - read(group,"rho_g21/$p_i/C_re") - read(group,"rho_g22/$p_i/C_re"))/2
    end
    CorrT1 = (read(group,"t1_g1/$p_i/C_im") + read(group,"t1_g2/$p_i/C_im"))/sqrt(2)
    CorrT2 = (read(group,"t2_g1/$p_i/C_im") + read(group,"t2_g2/$p_i/C_im"))/sqrt(2)
    CorrR1 =  read(group,"r1/$p_i/C_re")
    CorrR2 =  read(group,"r2/$p_i/C_re")
    CorrR3 =  read(group,"r3/$p_i/C_re")
    CorrR4 =  read(group,"r4/$p_i/C_re")
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function correlatorsp011(h5dset,ens;p=1)
    p_e = "(0,$p,$p)"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]
    CorrD1 =  read(group,"d/$p_i/C_re")
    CorrD2 =  read(group,"d/$p_0/C_re")
    Corrπ  =  read(group,"pi/$p_i/C_re")
    # TODO: Visually inspect sign of the individual components
    Corrρ = try
        (read(group,"rho_g3_g3/$p_i/C_re") + read(group,"rho_g3_g2/$p_i/C_re") - read(group,"rho_g2_g3/$p_i/C_re") - read(group,"rho_g2_g2/$p_i/C_re"))/2
    catch
        (read(group,"rho_g33/$p_i/C_re") + read(group,"rho_g32/$p_i/C_re") - read(group,"rho_g23/$p_i/C_re") - read(group,"rho_g22/$p_i/C_re"))/2
    end
    CorrT1 = (read(group,"t1_g3/$p_i/C_im") + read(group,"t1_g2/$p_i/C_im"))/sqrt(2)
    CorrT2 = (read(group,"t2_g3/$p_i/C_im") + read(group,"t2_g2/$p_i/C_im"))/sqrt(2)
    CorrR1 =  read(group,"r1/$p_i/C_re")
    CorrR2 =  read(group,"r2/$p_i/C_re")
    CorrR3 =  read(group,"r3/$p_i/C_re")
    CorrR4 =  read(group,"r4/$p_i/C_re")
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function correlatorsp000(h5dset,ens;p=1)
    Corrπ  = h5dset["$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re"][]
    Corrρ  = h5dset["$ens/p(0,0,0)/rho_g1/p_diag(0,0,0)/C_re"][]
    return Corrπ, Corrρ
end