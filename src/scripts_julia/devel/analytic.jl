using LinearAlgebra
m_p(m,p)       = m + sum( 1.0 .- cos.(p))
N_p(m,p)       = m_p(m,p)^2 + sum(sin.(p).^2)
ap_dot_aq(p,q) = dot(sin.(p),sin.(q))
unitP_dot_ap(P,q) = dot(P,sin.(q))/norm(P)

function analytic_free_pion(T,L,m,P)
    mom_t = (0:T-1).*(2π/T)
    mom_s = (0:L-1).*(2π/L)
    C = zeros(ComplexF64,T)
    for t in 0:T-1, q0 in mom_t, p0 in mom_t, p1 in mom_s, p2 in mom_s, p3 in mom_s
        q = [q0, p1 + P[2], p2 + P[3], p3 + P[4]]
        p = [p0, p1       , p2       , p3       ]
        C[t+1] += exp(+im*(p0-q0)*t)*(m_p(m,p)*m_p(m,q) + ap_dot_aq(p,q))/N_p(m,p)/N_p(m,q)
    end
    return @. 4C /(L^3*T^2)
end
function analytic_free_vector(T,L,m,P)
    mom_t = (0:T-1).*(2π/T)
    mom_s = (0:L-1).*(2π/L)
    C = zeros(ComplexF64,T)
    for t in 0:T-1, q0 in mom_t, p0 in mom_t, p1 in mom_s, p2 in mom_s, p3 in mom_s
        q = [q0, p1 + P[2], p2 + P[3], p3 + P[4]]
        p = [p0, p1       , p2       , p3       ]
        C[t+1] += exp(+im*(p0-q0)*t)*(m_p(m,p)*m_p(m,q) -2unitP_dot_ap(P,q)*unitP_dot_ap(P,p) + ap_dot_aq(p,q))/N_p(m,p)/N_p(m,q)
    end
    return @. 4C /(L^3*T^2)
end
function analytic_free_tensor_vector(T,L,m,P)
    mom_t = (0:T-1).*(2π/T)
    mom_s = (0:L-1).*(2π/L)
    C = zeros(ComplexF64,T)
    for t in 0:T-1, q0 in mom_t, p0 in mom_t, p1 in mom_s, p2 in mom_s, p3 in mom_s
        q = [q0, p1 + P[2], p2 + P[3], p3 + P[4]]
        p = [p0, p1       , p2       , p3       ]
        C[t+1] += exp(+im*(p0-q0)*t)*(m_p(m,q)*sin(p[1]) - m_p(m,p)*sin(q[1]))/N_p(m,p)/N_p(m,q)
    end
    return @. 4im*C /(L^3*T^2)
end
