using Pkg; Pkg.activate(".")
using ScatteringI1
using Plots
using LaTeXStrings
include("analytic.jl")
gr(frame=:box,legendfontsize=12)

T, L = 16, 16
t = 1:T
P = [0,0,0,1*2π/L]
m = -0.6

C_PS = analytic_free_pion(T,L,m,P)
C_V  = analytic_free_vector(T,L,m,P)
C_TV = analytic_free_tensor_vector(T,L,m,P)

scatter(1:T,(real.(C_PS)),label="real(C_PS)")
scatter!(1:T ,(real.(C_V)) ,label="real(C_V)")
scatter!(1:T ,(imag.(C_V)) ,label="imag(C_V)")
scatter!(1:T ,(imag.(C_PS)),label="imag(C_PS)")
scatter!(1:T ,(real.(C_TV)),label="real(C_T_V)")
scatter!(1:T ,(imag.(C_TV)),label="imag(C_T_V)")
