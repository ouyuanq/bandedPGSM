include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE epsilon * u'' - x * u = 0, u(-1) = Ai(-(1/epsilon)^(1/3)), u(1) = Ai((1/epsilon)^(1/3)), where the exact solution u = Ai((1/epsilon)^(1/3) * x)

T = Float64
airy_coeffs, airy_v = airy(T)
airy_coeffs_GSBSPG, airy_v_GSBSPG = airy_GSBSPG(T)
f = x -> 0*x

# solving speed
nvec = 2 .^ (4:16)
time_solve = zeros(length(nvec), 4)
for i in eachindex(nvec)
    @printf "Time No.%i\n" i

    n = nvec[i]
    airy_Wtrial, airy_Wtest, airy_Omega = airy(T, n)

    # banded PG method
    ben = @benchmark bandedPGsolve($(airy_coeffs), $(airy_v), $(airy_Wtrial), $(airy_Wtest), $(airy_Omega), $(f))
    time_solve[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-2, n, 2)
    ben = @benchmark GSBSPG_Chebyshev_solve($(T), $(airy_coeffs_GSBSPG), $(airy_Wtrial), $(airy_v_GSBSPG), $(fc))
    time_solve[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    ben = @benchmark GSBSPG_Chebyshev_NI_solve($(T), $(airy_coeffs_GSBSPG), $(airy_Wtrial), $(airy_v_GSBSPG), $(fc))
    time_solve[i, 3] = minimum(ben).time / 1e9

    # ultraspherical spectral method
    airy_coeffs_ultraS, airy_bc, airy_bcvals = airy_ultraS(T, n)
    ben = @benchmark ultraSsolve($(airy_coeffs_ultraS), $(airy_bc), $(airy_bcvals), $(f))
    time_solve[i, 4] = minimum(ben).time / 1e9
end
open("examples/ex2_time.txt", "w") do io
    writedlm(io, [nvec time_solve])
end

# accuracy
nvec = [2 .^ (4:14); 19500:50:20200; 2 .^ (15:16)]
accuracy = zeros(length(nvec), 4)
ue = x -> airyai(1e3 * x)
for i in eachindex(nvec)
    @printf "Accuracy No.%i\n" i

    n = nvec[i]
    airy_Wtrial, airy_Wtest, airy_Omega = airy(T, n)

    # banded PG method
    u = bandedPGsolve(airy_coeffs, airy_v, airy_Wtrial, airy_Wtest, airy_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-2, n, 2)
    u = GSBSPG_Chebyshev_solve(T, airy_coeffs_GSBSPG, airy_Wtrial, airy_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    u = GSBSPG_Chebyshev_NI_solve(T, airy_coeffs_GSBSPG, airy_Wtrial, airy_v_GSBSPG, fc)
    accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])

    # ultraspherical spectral method
    airy_coeffs_ultraS, airy_bc, airy_bcvals = airy_ultraS(T, n)
    u = ultraSsolve(airy_coeffs_ultraS, airy_bc, airy_bcvals, f)
    accuracy[i, 4] = Chebyshev_L2error(u, ue, nvec[end])
end
open("examples/ex2_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end