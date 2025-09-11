using Pkg; Pkg.activate("src/src_jl")
using HDF5
using LatticeUtils
using DelimitedFiles
using ArgParse: ArgParseSettings, parse_args, @add_arg_table

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
    @. new[op1,r1] = old[op2,r1]
    @. new[op2,r1] = old[op1,r1]
    @. new[op1,r2] = old[op2,r2]
    @. new[op2,r2] = old[op1,r2]
    return new
end

function swap_eigvals_cov(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = t:T-t+2
    @. new[op1,op1,r1] = old[op2,op2,r1]
    @. new[op2,op2,r1] = old[op1,op1,r1]
    @. new[op1,op2,r1] = old[op2,op1,r1]
    @. new[op2,op1,r1] = old[op1,op2,r1]
    return new
end
function swap_eigvals_cov(old,op1,op2,t)
    T = size(old,2)
    new = copy(old)
    r1 = 1:t-1
    r2 = T-t+2:T
    # TODO: This is probably not correct for elements
    # involving a potentially unswapped operator
    @. new[op1,op1,r1] = old[op2,op2,r1]
    @. new[op2,op2,r1] = old[op1,op1,r1]
    @. new[op1,op2,r1] = old[op2,op1,r1]
    @. new[op2,op1,r1] = old[op1,op2,r1]
    @. new[op1,op1,r2] = old[op2,op2,r2]
    @. new[op2,op2,r2] = old[op1,op1,r2]
    @. new[op1,op2,r2] = old[op2,op1,r2]
    @. new[op2,op1,r2] = old[op1,op2,r2]
    return new
end

function test_swaps(h5file,csvdat)
    fid  = h5open(h5file,"cw")
    data = readdlm(csvdat,',',skipstart=1)

    for row in eachrow(data)

        # determine which swaps to perform
        r, p, irrep, id = row[1:4]
        op1,op2,t_swap = row[5:7]

        # load eigenvalues
        T   = read(fid[r],"lattice")[1]
        ev  = read(fid[r]["$p/$irrep/$id"],"eigvals_3x3")
        Δev = read(fid[r]["$p/$irrep/$id"],"Delta_eigvals_3x3")
        ev_cov = read(fid[r]["$p/$irrep/$id"],"cov_eigvals_3x3")
        
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

        # replace old eigenvalues with the relabelled ones
        fid[r]["$p/$irrep/$id/eigvals_3x3"][:,:] = ev[:,:]
        fid[r]["$p/$irrep/$id/Delta_eigvals_3x3"][:,:] = Δev[:,:]
        fid[r]["$p/$irrep/$id/cov_eigvals_3x3"][:,:,:] = ev_cov[:,:,:]
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the parsed data"
        required = true
        "--swapinfo"
        help = "CSVs containing the relabelling information"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    h5file = args["h5file"]
    csvdat = args["swapinfo"]
    test_swaps(h5file,csvdat)
end
main()