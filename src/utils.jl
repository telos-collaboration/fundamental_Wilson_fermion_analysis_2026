function _parse_momentum(p0)
    rx = r"\(([0-9]),([0-9]),([0-9])\)"
    m = match(rx,p0)
    return parse.(Int,m.captures) 
end
function _are_permutations(p0,p1)
    return sort(p0) == sort(p1) 
end
function non_interacting_energy_2P(mπ,Δmπ,p2,L)
    E1, ΔE1 = non_interacting_energy_1P(mπ,Δmπ,p2,L)
    E   = E1 + mπ
    ΔE  = sqrt(Δmπ^2 + ΔE1^2)
    return E, ΔE
end
function non_interacting_energy_1P(mπ,Δmπ,p2,L)
    E1  = sqrt(mπ^2 + p2*(2*pi/L)^2)
    ΔE1 = Δmπ*mπ/E1
    return E1, ΔE1
end
