function pipi_correlator(CorrD1,CorrD2,CorrR1,CorrR2,CorrR3,CorrR4,L)
    L3, L6 = L^3, L^6
    Corr2π = @. (CorrD1 - CorrD2)/L6 + (CorrR1 + CorrR2 - CorrR3 - CorrR4)/L3
    return Corr2π 
end
# Note, I only include here the non-vanishing imaginary and real parts, respectively. 
# The reason is, that I do not want to include noisy matrix elements if I know that the diagram 
# must vanish analytically.  
function pipi_rho_matrix(Corr2π,Corrρ,CorrT1,CorrT2,L)
    N, nhits, T = size(Corr2π)
    L3, L6 = L^3, L^6
    corr = zeros(ComplexF64,(2,2,N,nhits,T))
    @assert size(Corr2π) == size(Corrρ) == size(CorrT1) == size(CorrT2) 
    corr[2,2,:,:,:] =  @. Corrρ/L3 + 0*im
    corr[2,1,:,:,:] =  @. 0        + im*(CorrT1-CorrT2)/L3
    corr[1,2,:,:,:] =  @. 0        + im*(CorrT2-CorrT1)/L3
    corr[1,1,:,:,:] =  @. Corr2π   + 0*im
    return corr
end
function pipi_rho_matrix_3x3_extension(Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2,L)
    @assert size(Corr_γ0γi_γi) == size(Corr_γi_γ0γi) == size(Corr_γ0γi_γ0γi) == size(Corrγ0γiT1) == size(Corrγ0γiT2)
    N, nhits, T = size(Corr_γ0γi_γi)
    L3, L6 = L^3, L^6
    corr_ext = zeros(ComplexF64,(3,3,N,nhits,T))
    corr_ext[1,3,:,:,:] = @. -im*(Corrγ0γiT1-Corrγ0γiT2)/L3 
    corr_ext[3,1,:,:,:] = @. -im*(Corrγ0γiT2-Corrγ0γiT1)/L3 
    corr_ext[2,3,:,:,:] = @. Corr_γi_γ0γi/L3
    corr_ext[3,2,:,:,:] = @. Corr_γ0γi_γi/L3
    corr_ext[3,3,:,:,:] = @. Corr_γ0γi_γ0γi/L3
    return corr_ext
end
function swap_eigval_numbering(old,t0)
    Nops, T   = size(old)[1], size(old)[3]
    new = copy(old)
    r1 = 1:t0-1
    r2 = T-t0+2:T
    if Nops == 2
        @. new[1,:,r1] = old[2,:,r1]
        @. new[2,:,r1] = old[1,:,r1]
        @. new[1,:,r2] = old[2,:,r2]
        @. new[2,:,r2] = old[1,:,r2]
    elseif Nops == 3

    end
    return new
end
function _preprocess_correlator(Corr;deriv)
    if deriv
        Corr = correlator_derivative(Corr;t_dim=4)
    end
    return Corr
end
function variational_analysis(Corr;t0,deriv=true,gevp)
    Corr = _preprocess_correlator(Corr;deriv)
    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0,gevp)
    # eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0)
    eigvals, Δeigvals = LatticeUtils.apply_jackknife(eigvals_resamples;dims=2)
    eigvals_cov = LatticeUtils.cov_jackknife_eigenvalues(eigvals_resamples)
    return eigvals, Δeigvals, eigvals_cov
end
function effective_masses(Corr;t0,deriv,gevp)
    Corr = _preprocess_correlator(Corr;deriv)
    eigvals_resamples = eigenvalues_jackknife_samples(Corr;t0,gevp)
    eigvals_resamples = swap_eigval_numbering(eigvals_resamples, t0)
    meff, Δmeff = LatticeUtils.log_meff_jackknife(real.(eigvals_resamples))
    return meff, Δmeff
end