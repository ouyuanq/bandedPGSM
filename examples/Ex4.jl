include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE u'' - cos(4 * sin(3*pi*x^2 + 1)) * u = - cos(4 * sin(3*pi*x^2 + 1))*sin(cos(20*pi*x) + 1) - 400*pi^2*sin(cos(20*pi*x) + 1)*sin(20*pi*x)^2 - 400*pi^2*cos(cos(20*pi*x) + 1)*cos(20*pi*x), u(-1) = u (1) = sin(2), s.t. u = sin(cos(20*pi*x) + 1)

# set up for different methods
T = Float64
composite_coeffs, composite_v = composite(T)
composite_coeffs_GSBSPG, composite_v_GSBSPG = composite_GSBSPG(T)
composite_coeffs_funcs, composite_v_GSBSPG_NI = composite_GSBSPG_NI(T)
N = length(composite_coeffs) - 1
ue = x -> sin(cos(20*pi*x) + 1)
f = x -> -cos(4 * sin(3*pi*x^2 + 1))*sin(cos(20*pi*x) + 1) - 400*pi^2*sin(cos(20*pi*x) + 1)*sin(20*pi*x)^2 - 400*pi^2*cos(cos(20*pi*x) + 1)*cos(20*pi*x)

# accuracy
nvec = [2 .^ (3:7); 200; 300; 400; 600; 700; 900; 2 .^(10:13)]
accuracy = zeros(length(nvec), 3)
@printf "solution accuracy of a second order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    composite_Wtrial, composite_Wtest, composite_Omega = composite(T, n)

    # banded PG method (test and trial space are the dual)
    u = bandedPGsolve(composite_coeffs, composite_v, composite_Wtrial, composite_Wtest, composite_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n, n, N)
    u = GSBSPG_Chebyshev_solve(T, composite_coeffs_GSBSPG, composite_Wtrial, composite_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    u = GSBSPG_Chebyshev_NI_solve(T, composite_coeffs_funcs, composite_Wtrial, composite_v_GSBSPG_NI, fc)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])

    @printf "%4i   %.2e   %.2e   %.2e\n" n accuracy[i, 1] accuracy[i, 2] accuracy[i, 3]
end
open("examples/ex4_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end

# construction speed
nvec = 2 .^ (3:13)
time_construct = zeros(length(nvec), 3)
@printf "construction cost of a second order equation:\n"
@printf "   n   bandedPG    MPG(R)     MPG(NI)\n"
for i in eachindex(nvec)
    n = nvec[i]
    composite_Wtrial, composite_Wtest, composite_Omega = composite(T, n)

    # banded PG method
    ben = @benchmark bandedPGmatrix_Chebyshev($(composite_coeffs), $(composite_Wtrial), $(composite_Wtest), $(composite_Omega))
    time_construct[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    fc = zeros(T, n-2)
    ben = @benchmark GSBSPG_Chebyshev($(T), $(composite_coeffs_GSBSPG), $(composite_Wtrial), $(composite_v_GSBSPG), $(fc))
    time_construct[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    ben = @benchmark GSBSPG_Chebyshev_NI($(T), $(composite_coeffs_funcs), $(composite_Wtrial), $(composite_v_GSBSPG_NI), $(fc))
    time_construct[i, 3] = minimum(ben).time / 1e9

    @printf "%4i   %.2e   %.2e   %.2e\n" n time_construct[i, 1] time_construct[i, 2] time_construct[i, 3]
end
open("examples/ex4_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end