include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE u'' - cos(sin(2*pi*x + 1)) * u = - cos(sin(2*pi*x + 1))*sin(cos(20*pi*x) + 1) - 400*pi^2*sin(cos(20*pi*x) + 1)*sin(20*pi*x)^2 - 400*pi^2*cos(cos(20*pi*x) + 1)*cos(20*pi*x), u(-1) = u (1) = sin(2), s.t. u = sin(cos(20*pi*x) + 1)

T = Float64
composite_coeffs, composite_v = composite(T)
composite_coeffs_GSBSPG, composite_v_GSBSPG = composite_GSBSPG(T)
composite_funcs = Vector{Function}(undef, 3)
composite_funcs[1] = x -> -cos(sin(2*pi*x + 1))
composite_funcs[3] = x -> 1
N = length(composite_funcs) - 1
# parameter for determining bandwidths in final matrix derived from numerical integration
ql = maximum(length.(composite_coeffs_GSBSPG) .+ (N-1:-1:-1)) 

# construction speed
nvec = 2 .^ (4:13)
time_construct = zeros(length(nvec), 3)
for i in eachindex(nvec)
    @printf "Time No.%i\n" i

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
    ben = @benchmark GSBSPG_Chebyshev_NI($(T), $(composite_funcs), $(composite_Wtrial), $(composite_v_GSBSPG), $(fc), $(ql))
    time_construct[i, 3] = minimum(ben).time / 1e9
end
open("examples/ex4_time.txt", "w") do io
    writedlm(io, [nvec time_construct])
end

# accuracy
ue = x -> sin(cos(20*pi*x) + 1)
f = x -> - cos(sin(2*pi*x + 1))*sin(cos(20*pi*x) + 1) - 400*pi^2*sin(cos(20*pi*x) + 1)*sin(20*pi*x)^2 - 400*pi^2*cos(cos(20*pi*x) + 1)*cos(20*pi*x)
nvec = [2 .^ (4:8); 300:100:800; 2 .^(10:13)]
accuracy = zeros(length(nvec), 3)

for i in eachindex(nvec)
    @printf "Accuracy No.%i\n" i

    n = nvec[i]
    composite_Wtrial, composite_Wtest, composite_Omega = composite(T, n)

    # banded PG method (test and trial space are the dual)
    u = bandedPGsolve(composite_coeffs, composite_v, composite_Wtrial, composite_Wtest, composite_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-N, n, N)
    u = GSBSPG_Chebyshev_solve(T, composite_coeffs_GSBSPG, composite_Wtrial, composite_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    u = GSBSPG_Chebyshev_NI_solve(T, composite_funcs, composite_Wtrial, composite_v_GSBSPG, fc, ql)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])
end
open("examples/ex4_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end