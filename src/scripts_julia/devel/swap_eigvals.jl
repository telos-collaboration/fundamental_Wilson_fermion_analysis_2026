using Pkg; Pkg.activate("src/src_jl")
using HDF5
using LatticeUtils
using Plots
using DelimitedFiles
gr(fontfamily="Computer Modern",frame=:box,markeralpha=0.7,titlefontsize=11)

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

function test_swaps(h5file,csvdat,irrep,p)
    fid  = h5open(h5file)
    runs = keys(fid)
    data = readdlm(csvdat,',',skipstart=1)

    for r in runs 
        ids = keys(fid[r][p][irrep])
        gevp_ids = filter(contains("gevp"),ids)
        for id in gevp_ids
            
            # load eigenvalues
            T   = read(fid[r],"lattice")[1]
            ev  = read(fid[r]["$p/$irrep/$id"],"eigvals_3x3")
            Δev = read(fid[r]["$p/$irrep/$id"],"Delta_eigvals_3x3")

            # determine which swaps to perform
            for row in eachrow(data)
                if row[1:4] == [r, p, irrep, id]
                    op1,op2,t_swap = row[5:7]
                    @show id,op1,op2,t_swap
                    if t_swap > 0
                        ev = swap_eigvals(ev,op1,op2,t_swap)
                        Δev = swap_eigvals(Δev,op1,op2,t_swap)
                    end
                    if t_swap < 0
                        ev = swap_eigvals_start(ev,op1,op2,abs(t_swap))
                        Δev = swap_eigvals_start(Δev,op1,op2,abs(t_swap))
                    end
                end
            end

            # plot for visual inspection
            t   = filter(!isequal(T÷2+1),1:T)
            plt = plot(yscale=:log10,legend=:top)
            plot_correlator!(plt,1:T,abs.(ev[1,t]),Δev[1,t],markersize=3,markershape=:rect,label="eigval #1 (3x3)")
            plot_correlator!(plt,1:T,abs.(ev[2,t]),Δev[2,t],markersize=3,markershape=:rect,label="eigval #2 (3x3)")
            plot_correlator!(plt,1:T,abs.(ev[3,t]),Δev[3,t],markersize=3,markershape=:rect,label="eigval #3 (3x3)")
            display(plt)
        end
    end
end

irrep  = "A1"
p      = "p(0,0,1)"
h5file = "data_assets/isospin1_eigenvalues_evp_gevp.hdf5"
csvdat = "metadata/swaps.csv"
test_swaps(h5file,csvdat,irrep,p)
