function correlators_xyz(h5dset,ens;p::Vector{Int})
    p_e = "($(p[1]),$(p[2]),$(p[3]))"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]

    # The following entries are always the same
    Corrπ  = read(group,"pi/$p_i/C_re")
    CorrD1 = read(group,"d/$p_i/C_re")
    CorrD2 = read(group,"d/$p_0/C_re")
    CorrR1 = read(group,"r1/$p_i/C_re")
    CorrR2 = read(group,"r2/$p_i/C_re")
    CorrR3 = read(group,"r3/$p_i/C_re")
    CorrR4 = read(group,"r4/$p_i/C_re")

    # The operators that result from inclusion of vector-meson operators depend on the external momentum which determines the required polarizations
    normρ  = dot(p,p)
    normT  = norm(p)
    Corrρ  = zero(Corrπ)
    CorrT1 = zero(Corrπ)
    CorrT2 = zero(Corrπ)

    # temporary array for reading data from hdf5 with fewer allocation
    tmp = zero(Corrπ)

    for i in 1:3, j in 1:3
        sign = ifelse(i == 2 , -1 , +1)
        pref = p[i]*p[j]*sign/normρ
        try
            Corrρ .= Corrρ .+ pref .* copyto!(tmp,group["rho_g$(i)_g$(j)/$p_i/C_re"])
        catch
            Corrρ .= Corrρ .+ pref .* copyto!(tmp,group["rho_g$(i)$(j)/$p_i/C_re"])
        end
    end
    for i in 1:3
        pref = p[i]/normT
        CorrT1 .= CorrT1 .+ pref .* copyto!(tmp,group["t1_g$(i)/$p_i/C_im"])
        CorrT2 .= CorrT2 .+ pref .* copyto!(tmp,group["t2_g$(i)/$p_i/C_im"])
    end

    return Corrπ, Corrρ, CorrT1, CorrT2, CorrR1, CorrR2, CorrR3, CorrR4, CorrD1, CorrD2
end

function correlators_xyz_3x3(h5dset,ens;p::Vector{Int})
    p_e = "($(p[1]),$(p[2]),$(p[3]))"
    p_i = "p_diag$p_e"
    p_0 = "p_diag(0,0,0)"
    group  = h5dset["$ens/p$p_e"]

    # The operators that result from inclusion of vector-meson operators depend on the external momentum which determines the required polarizations
    normρ  = dot(p,p)
    normT  = norm(p)
    shape  = size(read(group,"rho_g0g1_g1/$p_i/C_re"))

    Corr_γ0γi_γi = zeros(shape)
    Corr_γi_γ0γi = zeros(shape)
    Corr_γ0γi_γ0γi = zeros(shape)
    Corrγ0γiT1 = zeros(shape)
    Corrγ0γiT2 = zeros(shape)

    # temporary array for reading data from hdf5 with fewer allocation
    tmp = zeros(shape)
    
    for i in 1:3, j in 1:3
        sign = ifelse(i == 2 , -1 , +1)
        pref = p[i]*p[j]*sign/normρ
        Corr_γ0γi_γi   .= Corr_γ0γi_γi   .+ pref .* copyto!(tmp,group["rho_g0g$(i)_g$(j)/$p_i/C_re"])
        Corr_γi_γ0γi   .= Corr_γi_γ0γi   .+ pref .* copyto!(tmp,group["rho_g$(i)_g0g$(j)/$p_i/C_re"])
        Corr_γ0γi_γ0γi .= Corr_γ0γi_γ0γi .+ pref .* copyto!(tmp,group["rho_g0g$(i)_g0g$(j)/$p_i/C_re"])
    end
    for i in 1:3
        pref = p[i]/normT
        Corrγ0γiT1 .= Corrγ0γiT1 .+ pref .* copyto!(tmp,group["t1_g0g$(i)/$p_i/C_im"])
        Corrγ0γiT2 .= Corrγ0γiT2 .+ pref .* copyto!(tmp,group["t2_g0g$(i)/$p_i/C_im"])
    end

    return Corr_γ0γi_γi, Corr_γi_γ0γi, Corr_γ0γi_γ0γi, Corrγ0γiT1, Corrγ0γiT2
end

function correlatorsp000(h5dset,ens)
    Corrπ  = read(h5dset,"$ens/p(0,0,0)/pi/p_diag(0,0,0)/C_re")
    Corrρ  = read(h5dset,"$ens/p(0,0,0)/rho_g1/p_diag(0,0,0)/C_re")
    return Corrπ, Corrρ
end

# check momentum configuration
function momentum_is_nnn(p)
    @assert length(p) == 3
    return (allequal(p) && p[1] > 0) 
end
function momentum_is_00n(p)
    @assert length(p) == 3
    vanishing_indices = findall(iszero,p)
    return length(vanishing_indices) == 2
end
function momentum_is_0nn(p)
    @assert length(p) == 3
    vanishing_indices = findall(iszero,p)
    nonvanishing_indices = filter( i -> i ∉ vanishing_indices, eachindex(p))
    if length(vanishing_indices) == 1 # have exactly one vanishing component
        if p[nonvanishing_indices[1]] == p[nonvanishing_indices[2]] 
            return true
        end
    end
    return false
end

# Additional lattice irreps for vector meson for Yannick
function correlatorE(h5dset,ens;p::Vector{Int})
    # check that we have either a p=(nnn) momentum configuration 
    # or a p=(00n) momentum configuration (up to a permutation thereof)
    if momentum_is_nnn(p)
        Corrρ = correlatorE_nnn(h5dset,ens;p)
        return Corrρ
    elseif momentum_is_00n(p)
        Corrρ = correlatorE_00n(h5dset,ens;p)
        return Corrρ
    end
    return nothing
end
function correlatorE_00n(h5dset,ens;p::Vector{Int})
    @assert momentum_is_00n(p)
    z     = only(findall(x->!iszero(x),p))
    zp    = (z+1) % 3
    p_e   = "($(p[1]),$(p[2]),$(p[3]))"
    try
        Corrρ = read(h5dset,"$ens/p$p_e/rho_g$(zp)_g$(zp)/p_diag$p_e/C_re")
        return Corrρ
    catch
        Corrρ = read(h5dset,"$ens/p$p_e/rho_g$(zp)$(zp)/p_diag$p_e/C_re")
        return Corrρ
    end
end
function correlatorB1(h5dset,ens;p::Vector{Int})
    if momentum_is_0nn(p)
        z     = only(findall(x->iszero(x),p))
        p_e   = "($(p[1]),$(p[2]),$(p[3]))"
        try 
            Corrρ = read(h5dset,"$ens/p$p_e/rho_g$(z)$(z)/p_diag$p_e/C_re")
            return Corrρ
        catch
            Corrρ = read(h5dset,"$ens/p$p_e/rho_g$(z)_g$(z)/p_diag$p_e/C_re")
            return Corrρ
        end
    else
        return nothing
    end
end
function correlatorE_nnn(h5dset,ens;p::Vector{Int})
    @assert momentum_is_nnn(p)
    p_e   = "($(p[1]),$(p[2]),$(p[3]))"
    # TODO: Fix normalization, check signs of i=2
    try
        Corrρ = 4read(h5dset,"$ens/p$p_e/rho_g11/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g13/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g12/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g31/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g21/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g23/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g32/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g22/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g33/p_diag$p_e/C_re")
        return Corrρ
    catch
        Corrρ = 4read(h5dset,"$ens/p$p_e/rho_g1_g1/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g1_g3/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g1_g2/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g3_g1/p_diag$p_e/C_re")
        Corrρ = Corrρ .- 2read(h5dset,"$ens/p$p_e/rho_g2_g1/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g2_g3/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g3_g2/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g2_g2/p_diag$p_e/C_re")
        Corrρ = Corrρ .+  read(h5dset,"$ens/p$p_e/rho_g3_g3/p_diag$p_e/C_re")
        return Corrρ
    end
end
# function correlatorB2(h5dset,ens;p::Vector{Int})              # not needed for now
#     if momentum_is_0nn(p)
#         x, y = findall(x->!iszero(x),p)
#         p_e  = "($(p[1]),$(p[2]),$(p[3]))"
#         # TODO: Fix normalization, check signs of i=2
#         try 
#             Corrρ =          read(h5dset,"$ens/p$p_e/rho_g$(x)$(x)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .+ read(h5dset,"$ens/p$p_e/rho_g$(y)$(y)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .- read(h5dset,"$ens/p$p_e/rho_g$(x)$(y)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .- read(h5dset,"$ens/p$p_e/rho_g$(y)$(x)/p_diag$p_e/C_re")
#             return Corrρ
#         catch
#             Corrρ =          read(h5dset,"$ens/p$p_e/rho_g$(x)_g$(x)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .+ read(h5dset,"$ens/p$p_e/rho_g$(y)_g$(y)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .- read(h5dset,"$ens/p$p_e/rho_g$(x)_g$(y)/p_diag$p_e/C_re")
#             Corrρ = Corrρ .- read(h5dset,"$ens/p$p_e/rho_g$(y)_g$(x)/p_diag$p_e/C_re")
#             return Corrρ
#         end
#     else
#         return nothing
#     end
# end