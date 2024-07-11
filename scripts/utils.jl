function errorstring(x,Δx;nsig=2)
    @assert Δx > 0
    sgn = x < 0 ? "-" : ""
    x = abs(x)  
    # round error part to desired number of signficant digits
    # convert to integer if no fractional part exists
    Δx_rounded = round(Δx,sigdigits=nsig) 
    # get number of decimal digits for x  
    floor_log_10 = floor(Int,log10(Δx))
    dec_digits   = (nsig - 1) - floor_log_10
    # round x, to desired number of decimal digits 
    # (standard julia function deals with negative dec_digits) 
    x_rounded = round(x,digits=dec_digits)
    # get decimal and integer part if there is a decimal part
    if dec_digits > 0
        digits_val = Int(round(x_rounded*10.0^(dec_digits)))
        digits_unc = Int(round(Δx_rounded*10.0^(dec_digits)))
        str_val = _insert_decimal(digits_val,dec_digits) 
        str_unc = _insert_decimal(digits_unc,dec_digits)
        str_unc = nsig > dec_digits ? str_unc : string(digits_unc)
        return sgn*"$str_val($str_unc)"
    else
        return sgn*"$(Int(x_rounded))($(Int(Δx_rounded)))"
    end
end
function _insert_decimal(val::Int,digits)
    str = lpad(string(val),digits,"0")
    pos = length(str) - digits
    int = rpad(str[1:pos],1,"0")
    dec = str[pos+1:end]
    inserted = int*"."*dec
    return inserted
end
function _average_correlator(Corr)
    N, nhits, T = size(Corr)
    CorrAvg = dropdims(mean(Corr,dims=2),dims=2)
    C  = dropdims(mean(CorrAvg,dims=1),dims=1)
    ΔC = dropdims(std(CorrAvg,dims=1),dims=1)./sqrt(N)
    return C, ΔC
end
function D1(h5dset,p)
    πp =  h5dset["E1/p$p/pi/p_diag$p/C_re"][]
    π0 =  h5dset["E1/p(0,0,0)/pi/p_diag(0,0,0)/C_re"][]
    return πp .* π0
end
function correlatorsp001(h5dset,ens;p=1)
    p1 = "(0,0,$p)"
    CorrD1 = D1(h5dset,p1)
    CorrD2 =  h5dset["$ens/p$p1/d/p_diag$p1/C_re"][]
    Corrπ  =  h5dset["$ens/p$p1/pi/p_diag$p1/C_re"][]
    Corrρ  =  h5dset["$ens/p$p1/rho_g33/p_diag$p1/C_re"][]
    CorrT1 = -h5dset["$ens/p$p1/t1_g3/p_diag$p1/C_im"][]
    CorrT2 =  h5dset["$ens/p$p1/t2_g3/p_diag$p1/C_im"][]
    CorrR1 =  h5dset["$ens/p$p1/r1/p_diag$p1/C_re"][]
    CorrR2 =  h5dset["$ens/p$p1/r2/p_diag$p1/C_re"][]
    CorrR3 =  h5dset["$ens/p$p1/r3/p_diag$p1/C_re"][]
    CorrR4 =  h5dset["$ens/p$p1/r4/p_diag$p1/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function correlatorsp110(h5dset,ens;p=1)
    p1 = "($p,$p,0)"
    CorrD1 = D1(h5dset,p1)
    CorrD2 = h5dset["$ens/p$p1/d/p_diag$p1/C_re"][]
    Corrπ  = h5dset["$ens/p$p1/pi/p_diag$p1/C_re"][]
    # THE FOLLOWING DOES NOT GIVE A SIGNAL. IS SOMETHING WRONG WITH THE PARSING?
    Corrρ  =  (h5dset["$ens/p$p1/rho_g11/p_diag$p1/C_re"][])
    #Corrρ  =  (h5dset["$ens/p$p1/rho_g22/p_diag$p1/C_re"][])
    CorrT1 = -(h5dset["$ens/p$p1/t1_g1/p_diag$p1/C_im"][]   + h5dset["E1/p$p1/t1_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrT2 =  (h5dset["$ens/p$p1/t2_g1/p_diag$p1/C_im"][]   + h5dset["E1/p$p1/t2_g2/p_diag$p1/C_im"][])/sqrt(2)
    CorrR1 =   h5dset["$ens/p$p1/r1/p_diag$p1/C_re"][]
    CorrR2 =   h5dset["$ens/p$p1/r2/p_diag$p1/C_re"][]
    CorrR3 =   h5dset["$ens/p$p1/r3/p_diag$p1/C_re"][]
    CorrR4 =   h5dset["$ens/p$p1/r4/p_diag$p1/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function correlatorsp000(h5dset,ens)
    Corrπ  = h5dset["$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re"][]
    Corrρ  = h5dset["$ens/p(0,0,0)/rho_g1/p_diag(0,0,0)/C_re"][]
    CorrD1 = D1(h5dset,"(0,0,0)")
    # THE FOLLOWING IS NOT FULLY SELFCONSISTENT
    # It is used to compare the correlator normalization
    CorrD2 =  h5dset["$ens/p(0,0,1)/d/p_diag(0,0,0)/C_re"][]
    CorrT1 = -h5dset["$ens/p(0,0,1)/t1_g1/p_diag(0,0,0)/C_im"][]  
    CorrT2 =  h5dset["$ens/p(0,0,1)/t2_g1/p_diag(0,0,0)/C_im"][] 
    CorrR1 =  h5dset["$ens/p(0,0,1)/r1/p_diag(0,0,0)/C_re"][]
    CorrR2 =  h5dset["$ens/p(0,0,1)/r2/p_diag(0,0,0)/C_re"][]
    CorrR3 =  h5dset["$ens/p(0,0,1)/r3/p_diag(0,0,0)/C_re"][]
    CorrR4 =  h5dset["$ens/p(0,0,1)/r4/p_diag(0,0,0)/C_re"][]
    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end
function read_hdf5_file(file,ens,p1,p)
    h5dset = h5open(file)
    T, L = h5dset["$ens/lattice"][1:2]
    if p1 == "(0,0,0)"
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp000(h5dset,ens)
    elseif p1 == "(0,0,1)" || "(0,0,2)" || "(0,0,3)"
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorsp00n(h5dset,ens;p)
    elseif p1 == "(1,1,0)" || "(2,2,0)" || "(3,3,0)"
        Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2 = correlatorspnn0(h5dset,ens;p)
    end
    return T, L, Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end