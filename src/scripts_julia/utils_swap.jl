function swap_eigvals(old,op1,op2,t)
    T = size(old,3)
    new = copy(old)
    r1 = t:T-t+2
    @. new[op1,:,r1] = old[op2,:,r1]
    @. new[op2,:,r1] = old[op1,:,r1]
    return new
end
function swap_eigvals_start(old,op1,op2,t)
    T = size(old,3)
    new = copy(old)
    r1 = 1:t-1
    r2 = T-t+2:T
    @. new[op1,:,r1] = old[op2,:,r1]
    @. new[op2,:,r1] = old[op1,:,r1]
    @. new[op1,:,r2] = old[op2,:,r2]
    @. new[op2,:,r2] = old[op1,:,r2]
    return new
end
function _swap_eigenvalues(eigenvalues_resamples, op1, op2, t_swap)
    # load eigenvalues
    if t_swap > 0
        eigenvalues_resamples = swap_eigvals(eigenvalues_resamples,op1,op2,t_swap)
    end
    if t_swap < 0
        eigenvalues_resamples = swap_eigvals_start(eigenvalues_resamples,op1,op2,abs(t_swap))
    end
    return eigenvalues_resamples
end
function swap_eigvals(eigenvalues_resamples,csvdat,r,p,irrep,id)
    data = readdlm(csvdat,',',skipstart=1)
    # determine which swaps to perform
    for row in eachrow(data)
        if row[1:4] == [r, p, irrep, id]
            op1, op2, t_swap = row[5:7]
            eigenvalues_resamples = _swap_eigenvalues(eigenvalues_resamples, op1, op2, t_swap)
        end
    end
    # replace old eigenvalues with the relabelled ones
    return eigenvalues_resamples
end