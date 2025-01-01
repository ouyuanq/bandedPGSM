include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODEs u'' - (1 + sin(x)) * u' + exp(x) * u = f, u(-1) = 1, u(1) = 1, where f = (exp((x^2 - 1) / 2)) * (1 + x^2 + exp(x) - x*(sin(x) + 1)) s.t. u = exp((x^2 - 1) / 2)

T = Float64
expx2_coeffs, expx2_v = expx2(T)
expx2_coeffs_GSBSPG, expx2_v_GSBSPG = expx2_GSBSPG(T)

# construction speed
nvec = 2 .^ (3:12)
time_construct = zeros(length(nvec), 3)
for i in eachindex(nvec)
    @printf "Time No.%i\n" i
    
    n = nvec[i]
    expx2_Wtrial, expx2_Wtest, expx2_Omega = expx2(T, n)

    # banded PG method
    ben = @benchmark bandedPGmatrix_Chebyshev($(expx2_coeffs), $(expx2_Wtrial), $(expx2_Wtest), $(expx2_Omega))
    time_construct[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    bc = zeros(T, 2)
    f = zeros(T, n)
    ben = @benchmark GSBSPG_Chebyshev($(T), $(expx2_coeffs_GSBSPG), $(expx2_Wtrial), $(expx2_v_GSBSPG), $(f))
    time_construct[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    ben = @benchmark GSBSPG_Chebyshev_NI($(T), $(expx2_coeffs_GSBSPG), $(expx2_Wtrial), $(expx2_v_GSBSPG), $(f))
    time_construct[i, 3] = minimum(ben).time / 1e9
end
open("examples/ex1_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end

# accuracy
nvec = [8:2:30; 2 .^ (5:12)]
accuracy = zeros(length(nvec), 3)
f = x -> (exp((x^2 - 1) / 2)) * (1 + x^2 + exp(x) - x*(sin(x) + 1))
ue = x -> (exp((x^2 - 1) / 2))
for i in eachindex(nvec)
    @printf "Accuracy No.%i\n" i

    n = nvec[i]
    expx2_Wtrial, expx2_Wtest, expx2_Omega = expx2(T, n)

    # banded PG method
    u = bandedPGsolve(expx2_coeffs, expx2_v, expx2_Wtrial, expx2_Wtest, expx2_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-2, n, 2)
    L1, fc1 = GSBSPG_Chebyshev(T, expx2_coeffs_GSBSPG, expx2_Wtrial, expx2_v_GSBSPG, copy(fc))
    u1 = expx2_Wtrial * (L1 \ fc1)
    u1[1:length(expx2_v_GSBSPG)] .+= expx2_v_GSBSPG
    accuracy[i, 2] = Chebyshev_L2error(u1, ue, nvec[end])

    # GSBSPG (numerical integration)
    L2, fc2 = GSBSPG_Chebyshev_NI(T, expx2_coeffs_GSBSPG, expx2_Wtrial, expx2_v_GSBSPG, copy(fc))
    u2 = expx2_Wtrial * (L2 \ fc2)
    u2[1:length(expx2_v_GSBSPG)] .+= expx2_v_GSBSPG
    accuracy[i, 3] = Chebyshev_L2error(u2, ue, nvec[end])
end
open("examples/ex1_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end