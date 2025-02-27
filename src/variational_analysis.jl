function pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
    L3, L6 = L^3, L^6
    Corr2π = @. (CorrD1 - CorrD2)/L6 + (CorrR1 + CorrR2 - CorrR3 - CorrR4)/L3
    return Corr2π 
end
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(ComplexF64,(2,2,N,nhits,T))
    @assert size(Corr2π) == size(Corrρ) == size(CorrT1) == size(CorrT2) 
    corr[1,1,:,:,:] =  @. Corrρ/L3 + 0*im
    corr[1,2,:,:,:] =  @. 0        + im*(CorrT1-CorrT2)/L3
    corr[2,1,:,:,:] =  @. 0        + im*(CorrT2-CorrT1)/L3
    corr[2,2,:,:,:] =  @. Corr2π   + 0*im
    return corr
end
function swap_eigval_numbering(old,t0)
    T   = size(old)[3]
    new = copy(old)
    @. new[1,:,1:t0-1] = old[2,:,1:t0-1]
    @. new[2,:,1:t0-1] = old[1,:,1:t0-1]
    @. new[1,:,T-t0+2:T] = old[2,:,T-t0+2:T]
    @. new[2,:,T-t0+2:T] = old[1,:,T-t0+2:T]
    return new
end
function _preprocess_correlator(Corr;deriv,maxhits,sign = +1)
    Corr  = correlator_folding(Corr;t_dim=4,sign=+1)
    if deriv
        Corr = correlator_derivative(Corr;t_dim=4)
        sign = -1
    end
    return Corr, sign
end
function variational_analysis(Corr;t0,maxhits=typemax(Int),deriv=true)
    Corr, sign = _preprocess_correlator(Corr;deriv,maxhits)
    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0)
    eigvals, Δeigvals = LatticeUtils.apply_jackknife(eigvals_resamples;dims=2)
    eigvals_cov = LatticeUtils.cov_jackknife_eigenvalues(eigvals_resamples)
    return eigvals, Δeigvals, eigvals_cov
end
function effective_masses(Corr;t0,deriv,maxhits=typemax(Int))
    Corr, sign = _preprocess_correlator(Corr;deriv,maxhits)
    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0)
    eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0)
    meff, Δmeff = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))
    return meff, Δmeff
end