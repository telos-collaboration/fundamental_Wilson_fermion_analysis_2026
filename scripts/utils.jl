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
