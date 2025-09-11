function swap_eigvals(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = t:T-t+2
    @. new[op1,r1] = old[op2,r1]
    @. new[op2,r1] = old[op1,r1]
    return new
end
function swap_eigvals_start(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = 1:t-1
    r2 = T-t+2:T
    #@. new[op1,r1] = old[op2,r1]
    @. new[op2,r1] = old[op1,r1]
    @. new[op1,r2] = old[op2,r2]
    @. new[op2,r2] = old[op1,r2]
    return new
end
function swap_eigvals_cov(old,op1,op2,t)
    Nops, T = size(old,1), size(old,3)
    new = copy(old)
    r1 = t:T-t+2
    #@. new[op1,:,r1] = old[op2,:,r1]
    @. new[op2,:,r1] = old[op1,:,r1]
    @. new[:,op1,r1] = old[:,op2,r1]
    @. new[:,op2,r1] = old[:,op1,r1]
    return new
end
function swap_eigvals_cov(old,op1,op2,t)
    Nops, T = size(old,1), size(old,3)
    new = copy(old)
    r1 = 1:t-1
    r2 = T-t+2:T
    @. new[op1,:,r1] = old[op2,:,r1]
    @. new[op2,:,r1] = old[op1,:,r1]
    @. new[op1,:,r2] = old[op2,:,r2]
    @. new[op2,:,r2] = old[op1,:,r2]
    @. new[:,op1,r1] = old[:,op2,r1]
    @. new[:,op2,r1] = old[:,op1,r1]
    @. new[:,op1,r2] = old[:,op2,r2]
    @. new[:,op2,r2] = old[:,op1,r2]
    return new
end
function _swap_eigenvalues(ev, Δev, ev_cov, op1, op2, t_swap)
    # load eigenvalues
    if t_swap > 0
        ev = swap_eigvals(ev,op1,op2,t_swap)
        Δev = swap_eigvals(Δev,op1,op2,t_swap)
        ev_cov = swap_eigvals_cov(ev_cov,op1,op2,t_swap)
    end
    if t_swap < 0
        ev = swap_eigvals_start(ev,op1,op2,abs(t_swap))
        Δev = swap_eigvals_start(Δev,op1,op2,abs(t_swap))
        ev_cov = swap_eigvals_cov(ev_cov,op1,op2,abs(t_swap))
    end
    return ev, Δev, ev_cov
end
function test_swaps(fid,csvdat,r,p,irrep,id)
    data = readdlm(csvdat,',',skipstart=1)
    # load data from hdf5 file
    ev  = read(fid[r]["$p/$irrep/$id"],"eigvals_3x3")
    Δev = read(fid[r]["$p/$irrep/$id"],"Delta_eigvals_3x3")
    ev_cov = read(fid[r]["$p/$irrep/$id"],"cov_eigvals_3x3")    
    # determine which swaps to perform
    for row in eachrow(data)
        if row[1:4] == [r, p, irrep, id]
            op1, op2, t_swap = row[5:7]
            ev, Δev, ev_cov = _swap_eigenvalues(ev, Δev, ev_cov, op1, op2, t_swap)
        end
    end
    # replace old eigenvalues with the relabelled ones
    fid[r]["$p/$irrep/$id/eigvals_3x3"][:,:] = ev[:,:]
    fid[r]["$p/$irrep/$id/Delta_eigvals_3x3"][:,:] = Δev[:,:]
    fid[r]["$p/$irrep/$id/cov_eigvals_3x3"][:,:,:] = ev_cov[:,:,:]
end
function swap_eigvals(ev, Δev, ev_cov,csvdat,r,p,irrep,id)
    data = readdlm(csvdat,',',skipstart=1)
    # determine which swaps to perform
    for row in eachrow(data)
        if row[1:4] == [r, p, irrep, id]
            op1, op2, t_swap = row[5:7]
            ev, Δev, ev_cov = _swap_eigenvalues(ev, Δev, ev_cov, op1, op2, t_swap)
        end
    end
    # replace old eigenvalues with the relabelled ones
    return ev, Δev, ev_cov
end