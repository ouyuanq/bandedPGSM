include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE u''' - (cos(x) + 2) * u' + 10 * exp(x) * u = x^3*exp(x^2/2 - 1/2) + 10*exp(x)*exp(x^2/2 - 1/2) + 3*x*exp(x^2/2 - 1/2) - x*exp(x^2/2 - 1/2)*(cos(x) + 2), u(-1) = u (1) = 1, u'(1) = 1, s.t. u = exp((x^2 - 1) / 2)

# set up for different methods
T = Float64
taylor_coeffs, taylor_v = taylor(T)
taylor_coeffs_GSBSPG, taylor_v_GSBSPG = taylor_GSBSPG(T)
taylor_coeffs_funcs, taylor_v_GSBSPG_NI = taylor_GSBSPG_NI(T)
N = length(taylor_coeffs) - 1
ue = x -> exp((x^2 - 1) / 2)
f = x -> x^3*exp(x^2/2 - 1/2) + 10*exp(x)*exp(x^2/2 - 1/2) + 3*x*exp(x^2/2 - 1/2) - x*exp(x^2/2 - 1/2)*(cos(x) + 2)

# accuracy
nvec = [8:2:24; 2 .^ (5:13)]
accuracy = zeros(length(nvec), 4)
@printf "solution accuracy of a third order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    taylor_R, taylor_Q, taylor_Omega = taylor(T, n)

    # banded PG method (test and trial space are the dual)
    u = bandedPGsolve(taylor_coeffs, taylor_v, taylor_R, taylor_Q, taylor_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n, n, N)
    u = GSBSPG_Chebyshev_solve(T, taylor_coeffs_GSBSPG, taylor_R, taylor_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    u = GSBSPG_Chebyshev_NI_solve(T, taylor_coeffs_funcs, taylor_R, taylor_v_GSBSPG_NI, fc)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])

    # # banded PG method (test and trial space are the same)
    # taylor_R, taylor_Q, taylor_Omega = thirdrightNeumann_Chebyshev_nondual(T, n)
    # u = bandedPGsolve(taylor_coeffs, taylor_v, taylor_R, taylor_Q, taylor_Omega, f)
    # accuracy[i, 4] = Chebyshev_L2error(u, ue, nvec[end])

    @printf "%4i   %.2e   %.2e   %.2e\n" n accuracy[i, 1] accuracy[i, 2] accuracy[i, 3]
end
open("examples/taylor_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end

# construction speed
nvec = 2 .^ (3:13)
time_construct = zeros(length(nvec), 3)
@printf "construction cost of a third order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    taylor_R, taylor_Q, taylor_Omega = taylor(T, n)

    # banded PG method
    ben = @benchmark bandedPGmatrix_Chebyshev($(taylor_coeffs), $(taylor_R), $(taylor_Q), $(taylor_Omega))
    time_construct[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    fc = zeros(T, n)
    ben = @benchmark GSBSPG_Chebyshev($(T), $(taylor_coeffs_GSBSPG), $(taylor_R), $(taylor_v_GSBSPG), $(fc))
    time_construct[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    ben = @benchmark GSBSPG_Chebyshev_NI($(T), $(taylor_coeffs_funcs), $(taylor_R), $(taylor_v_GSBSPG_NI), $(fc))
    time_construct[i, 3] = minimum(ben).time / 1e9

    @printf "%4i   %.2e   %.2e   %.2e\n" n time_construct[i, 1] time_construct[i, 2] time_construct[i, 3]
end
open("examples/taylor_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end