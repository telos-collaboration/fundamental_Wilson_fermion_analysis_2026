function _swap_eigenvalues(old,op1,op2,t)
    T = size(old,3)
    new = copy(old)
    # We use the sign of t to encode which parts of the correlator to relabel
    if t > 0
        r1 = t:T-t+2
    else
        r1 = vcat(collect(1:abs(t)-1), collect(T-abs(t)+2:T))
    end
    @. new[op1,:,r1] = old[op2,:,r1]
    @. new[op2,:,r1] = old[op1,:,r1]
    return new
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