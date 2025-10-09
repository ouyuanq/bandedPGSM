include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE u''' - cos(x) * u'' + 10 * exp(x) * u = exp((x^2 - 1)/2) * (3*x - cos(x) + 10*exp(x) - x^2*cos(x) + x^3), u(-1) = u (1) = 1, u'(1) = 1, s.t. u = exp((x^2 - 1) / 2)

T = Float64
expcos_coeffs, expcos_v = expcos(T)
expcos_coeffs_GSBSPG, expcos_v_GSBSPG = expcos_GSBSPG(T)
N = length(expcos_coeffs) - 1

# construction speed
nvec = 2 .^ (3:13)
time_construct = zeros(length(nvec), 3)
for i in eachindex(nvec)
    @printf "Time No.%i\n" i
    
    n = nvec[i]
    expcos_Wtrial, expcos_Wtest, expcos_Omega = expcos(T, n)

    # banded PG method
    ben = @benchmark bandedPGmatrix_Chebyshev($(expcos_coeffs), $(expcos_Wtrial), $(expcos_Wtest), $(expcos_Omega))
    time_construct[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    f = zeros(T, n-2)
    ben = @benchmark GSBSPG_Chebyshev($(T), $(expcos_coeffs_GSBSPG), $(expcos_Wtrial), $(expcos_v_GSBSPG), $(f))
    time_construct[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    ben = @benchmark GSBSPG_Chebyshev_NI($(T), $(expcos_coeffs_GSBSPG), $(expcos_Wtrial), $(expcos_v_GSBSPG), $(f))
    time_construct[i, 3] = minimum(ben).time / 1e9
end
open("examples/ex1_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end

# accuracy
nvec = [8:2:30; 2 .^ (5:13)]
accuracy = zeros(length(nvec), 4)
f = x -> exp((x^2 - 1)/2) * (3*x - cos(x) + 10*exp(x) - x^2*cos(x) + x^3)
ue = x -> exp((x^2 - 1) / 2)
for i in eachindex(nvec)
    @printf "Accuracy No.%i\n" i

    n = nvec[i]
    expcos_Wtrial, expcos_Wtest, expcos_Omega = expcos(T, n)

    # banded PG method (test and trial space are the dual)
    u = bandedPGsolve(expcos_coeffs, expcos_v, expcos_Wtrial, expcos_Wtest, expcos_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-N, n, N)
    u = GSBSPG_Chebyshev_solve(T, expcos_coeffs_GSBSPG, expcos_Wtrial, expcos_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    u = GSBSPG_Chebyshev_NI_solve(T, expcos_coeffs_GSBSPG, expcos_Wtrial, expcos_v_GSBSPG, fc)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])

    # # banded PG method (test and trial space are the same)
    # expcos_Wtrial, expcos_Wtest, expcos_Omega = thirdrightNeumann_ChebyshevW_nondual(T, n)
    # u = bandedPGsolve(expcos_coeffs, expcos_v, expcos_Wtrial, expcos_Wtest, expcos_Omega, f)
    # accuracy[i, 4] = Chebyshev_L2error(u, ue, nvec[end])
end
open("examples/ex1_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end