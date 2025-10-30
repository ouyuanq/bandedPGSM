include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE 0.1*u^(4) + x^10 * u'' + sin(40*x) * u = x^10*exp(sin(pi*x)) - sin(40*x)*(pi^2*exp(sin(pi*x))*sin(pi*x) - pi^2*exp(sin(pi*x))*cos(pi*x)^2) + (3*pi^4*exp(sin(pi*x))*sin(pi*x)^2)/10 + (pi^4*exp(sin(pi*x))*sin(pi*x))/10 - (2*pi^4*exp(sin(pi*x))*cos(pi*x)^2)/5 + (pi^4*exp(sin(pi*x))*cos(pi*x)^4)/10 - (3*pi^4*exp(sin(pi*x))*cos(pi*x)^2*sin(pi*x))/5, u(-1) = u (1) = 1, u'(-1) = u'(1) = -pi, s.t. u = exp(sin(pi*x))

# set up for different methods
T = Float64
oscillatory_coeffs, oscillatory_v = oscillatory(T)
oscillatory_coeffs_MPG, oscillatory_v_MPG = oscillatory_MPG(T)
oscillatory_coeffs_funcs, oscillatory_v_MPG_NI = oscillatory_MPG_NI(T)
N = length(oscillatory_coeffs) - 1
ue = x -> exp(sin(pi*x))
f = x -> x^10*exp(sin(pi*x)) - sin(40*x)*(pi^2*exp(sin(pi*x))*sin(pi*x) - pi^2*exp(sin(pi*x))*cos(pi*x)^2) + (3*pi^4*exp(sin(pi*x))*sin(pi*x)^2)/10 + (pi^4*exp(sin(pi*x))*sin(pi*x))/10 - (2*pi^4*exp(sin(pi*x))*cos(pi*x)^2)/5 + (pi^4*exp(sin(pi*x))*cos(pi*x)^4)/10 - (3*pi^4*exp(sin(pi*x))*cos(pi*x)^2*sin(pi*x))/5

# accuracy
nvec = [15; 40:13:105; 2 .^ (8:10); 2044; 4088; 8192]
accuracy = zeros(length(nvec), 3)
@printf "solution accuracy of a second order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    oscillatory_R, oscillatory_Q, oscillatory_Omega = oscillatory(T, n)

    # banded PG method (test and trial space are the dual)
    u = bandedPGsolve(oscillatory_coeffs, oscillatory_v, oscillatory_R, oscillatory_Q, oscillatory_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # MPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n, n, N)
    u = MPG_Chebyshev_solve(T, oscillatory_coeffs_MPG, oscillatory_R, oscillatory_v_MPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # MPG (numerical integration)
    u = MPG_Chebyshev_NI_solve(T, oscillatory_coeffs_funcs, oscillatory_R, oscillatory_v_MPG_NI, fc)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])

    @printf "%4i   %.2e   %.2e   %.2e\n" n accuracy[i, 1] accuracy[i, 2] accuracy[i, 3]
end
open("examples/oscillatory_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end

# construction speed
nvec = 2 .^ (3:13)
time_construct = zeros(length(nvec), 3)
@printf "construction cost of a second order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    oscillatory_R, oscillatory_Q, oscillatory_Omega = oscillatory(T, n)

    # banded PG method
    ben = @benchmark bandedPGmatrix_Chebyshev($(oscillatory_coeffs), $(oscillatory_R), $(oscillatory_Q), $(oscillatory_Omega))
    time_construct[i, 1] = minimum(ben).time / 1e9

    # MPG (recurrence)
    fc = zeros(T, n-2)
    ben = @benchmark MPG_Chebyshev($(T), $(oscillatory_coeffs_MPG), $(oscillatory_R), $(oscillatory_v_MPG), $(fc))
    time_construct[i, 2] = minimum(ben).time / 1e9

    # MPG (numerical integration)
    ben = @benchmark MPG_Chebyshev_NI($(T), $(oscillatory_coeffs_funcs), $(oscillatory_R), $(oscillatory_v_MPG_NI), $(fc))
    time_construct[i, 3] = minimum(ben).time / 1e9

    @printf "%4i   %.2e   %.2e   %.2e\n" n time_construct[i, 1] time_construct[i, 2] time_construct[i, 3]
end
open("examples/oscillatory_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end