function _parse_momentum(p0)
    rx = r"\(([0-9]),([0-9]),([0-9])\)"
    m = match(rx,p0)
    return parse.(Int,m.captures) 
end
function unique_momenta(p0)
    p = unique(sort,_parse_momentum.(p0))
    return ["p($(p[i][1]),$(p[i][2]),$(p[i][3]))" for i in eachindex(p)]
end
function non_interacting_energy_2P(mπ,Δmπ,p2,L)
    E1, ΔE1 = non_interacting_energy_1P(mπ,Δmπ,p2,L)
    E  = E1 + mπ
    ΔE = sqrt(Δmπ^2 + ΔE1^2)
    return E, ΔE
end
function non_interacting_energy_1P(mπ,Δmπ,p2,L)
    E1  = sqrt(mπ^2 + p2*(2*pi/L)^2)
    ΔE1 = Δmπ*mπ/E1
    return E1, ΔE1
end
function non_interacting_energy_2P_lattice(mπ,Δmπ,px,py,pz,L)
    E1, ΔE1 = non_interacting_energy_1P_lattice(mπ,Δmπ,px,py,pz,L)
    E  = E1 + mπ
    ΔE = sqrt(Δmπ^2 + ΔE1^2)
    return E, ΔE
end
function non_interacting_energy_1P_lattice(mπ,Δmπ,px,py,pz,L)
    xyz = 2*(sin(px*pi/L)^2 + sin(py*pi/L)^2 + sin(pz*pi/L)^2)
    E1  = acosh(cosh(mπ) + xyz )
    ΔE1 = Δmπ*sinh(mπ)/(sqrt(-1 + xyz + cosh(mπ))*sqrt(1 + xyz + cosh(mπ)))
    return E1, ΔE1
end
function _avg_sources(Corr;maxhits=typemax(Int))
    nhits = size(Corr)[4]
    h     = min(nhits,maxhits)
    Corr_avg = dropdims(mean(Corr[:,:,:,1:h,:],dims=4),dims=4)
    return Corr_avg, h
end
function read_correlation_matrix(h5dset,ens,p,label="correlation_matrix";maxhits=typemax(Int),average_equivalent_momenta)
    p0   = _parse_momentum(p)
    Corr = read(h5dset,joinpath(ens,p,label))
    
    Corr_avg, h = _avg_sources(Corr;maxhits)
    sources = [h]
    momenta = "($(p0[1])$(p0[2])$(p0[3])"
    
    if average_equivalent_momenta
        norm  = 1
        perms = unique(permutations(p0))
        # remove the do-nothing permutation because we already covered it 
        for perm in filter!(!isequal(p0),perms)
            key = "p($(perm[1]),$(perm[2]),$(perm[3]))"
            if haskey(h5dset[ens],joinpath(key,label))
                Corr = read(h5dset,joinpath(ens,key,label))
                Corr_tmp, h = _avg_sources(Corr;maxhits)

                @. Corr_avg = Corr_avg + Corr_tmp
                append!(sources, h)
                momenta = momenta * ",$(perm[1])$(perm[2])$(perm[3])"
                norm += 1
            end
        end
        @. Corr_avg = Corr_avg/norm 
    end
    momenta = momenta * ")"

    return Corr_avg, Tuple(sources), momenta
end
function read_pion_correlator(h5dset,ens,p;average_equivalent_momenta)
    # for the pion correlator: Average over equvivalent external momenta
    Corrπ = read(h5dset,joinpath(ens,p,"correlator_pion"))
    if average_equivalent_momenta
        label = "correlator_pion"
        perms = unique(permutations(p))
        norm = 0
        # remove the do-nothing permutation
        for perm in filter!(!isequal(p),perms)
            key = "p($(perm[1]),$(perm[2]),$(perm[3]))"
            if haskey(h5dset[ens],joinpath(key,label))
                Corr_perm = read(h5dset,joinpath(ens,key,label))
                @. Corrπ = Corrπ + Corr_perm
                norm += 1
            end
        end
        @. Corrπ = Corrπ/norm
    end
    return Corrπ
end