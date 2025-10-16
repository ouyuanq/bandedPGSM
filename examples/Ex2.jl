include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving ODE epsilon * u'' - x * u = 0, u(-1) = Ai(-(1/epsilon)^(1/3)), u(1) = Ai((1/epsilon)^(1/3)), where the exact solution u = Ai((1/epsilon)^(1/3) * x)

# set up for different methods
T = Float64
airy_coeffs, airy_v = airy(T)
airy_coeffs_GSBSPG, airy_v_GSBSPG = airy_GSBSPG(T)
airy_coeffs_funcs, airy_v_GSBSPG_NI = airy_GSBSPG_NI(T)
N = length(airy_coeffs) - 1
ue = x -> airyai(1e3 * x)
f = x -> 0*x

# accuracy
nvec = [2 .^ (4:14); 19600:100:20100; 2 .^ (15:16)]
accuracy = zeros(length(nvec), 5)
@printf "solution accuracy of Airy equation:\n"
@printf "    n   bandedPG    MPG(R)     MPG(NI)    ultraS\n"
for i in eachindex(nvec)
    n = nvec[i]
    airy_Wtrial, airy_Wtest, airy_Omega = airy(T, n)

    # banded PG method
    u = bandedPGsolve(airy_coeffs, airy_v, airy_Wtrial, airy_Wtest, airy_Omega, f)
    accuracy[i, 1] = Chebyshev_L2error(u, ue, nvec[end])

    airy_Wtest = BandedMatrix((0 => Ones(T, n-N),), (n, n-N), (N, 0))
    airy_Omega = Diagonal(Ones(T, n))

    u = bandedPGsolve(airy_coeffs, airy_v, airy_Wtrial, airy_Wtest, airy_Omega, f)
    accuracy[i, 5] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n, n, N)
    u = GSBSPG_Chebyshev_solve(T, airy_coeffs_GSBSPG, airy_Wtrial, airy_v_GSBSPG, fc)
    accuracy[i, 2] = Chebyshev_L2error(u, ue, nvec[end])

    # GSBSPG (numerical integration)
    if n <= 8192
        u = GSBSPG_Chebyshev_NI_solve(T, airy_coeffs_funcs, airy_Wtrial, airy_v_GSBSPG_NI, fc)
        accuracy[i, 3] = Chebyshev_L2error(u, ue, nvec[end])
    end

    # ultraspherical spectral method
    airy_coeffs_ultraS, airy_bc, airy_bcvals = airy_ultraS(T, n)
    u = ultraSsolve(airy_coeffs_ultraS, airy_bc, airy_bcvals, f)
    accuracy[i, 4] = Chebyshev_L2error(u, ue, nvec[end])

    @printf "%5i   %.2e   %.2e   %.2e   %.2e\n" n accuracy[i, 1] accuracy[i, 2] accuracy[i, 3] accuracy[i, 4]
end
open("examples/ex2_accuracy.txt", "w") do io
    writedlm(io, [nvec accuracy])
end

# solving speed
nvec = 2 .^ (4:16)
time_solve = zeros(length(nvec), 4)
@printf "construction and solution cost of Airy equation:\n"
@printf "    n   bandedPG    MPG(R)     MPG(NI)    ultraS\n"
for i in eachindex(nvec)
    n = nvec[i]
    airy_Wtrial, airy_Wtest, airy_Omega = airy(T, n)

    # banded PG method
    ben = @benchmark bandedPGsolve($(airy_coeffs), $(airy_v), $(airy_Wtrial), $(airy_Wtest), $(airy_Omega), $(f))
    time_solve[i, 1] = minimum(ben).time / 1e9

    # GSBSPG (recurrence)
    fc = Chebyshev_rhs_NI(T, f, n-N, n, N)
    ben = @benchmark GSBSPG_Chebyshev_solve($(T), $(airy_coeffs_GSBSPG), $(airy_Wtrial), $(airy_v_GSBSPG), $(fc))
    time_solve[i, 2] = minimum(ben).time / 1e9

    # GSBSPG (numerical integration)
    if n <= 8192
        ben = @benchmark GSBSPG_Chebyshev_NI_solve($(T), $(airy_coeffs_funcs), $(airy_Wtrial), $(airy_v_GSBSPG_NI), $(fc))
        time_solve[i, 3] = minimum(ben).time / 1e9
    end

    # ultraspherical spectral method
    airy_coeffs_ultraS, airy_bc, airy_bcvals = airy_ultraS(T, n)
    ben = @benchmark ultraSsolve($(airy_coeffs_ultraS), $(airy_bc), $(airy_bcvals), $(f))
    time_solve[i, 4] = minimum(ben).time / 1e9

    @printf "%5i   %.2e   %.2e   %.2e   %.2e\n" n time_solve[i, 1] time_solve[i, 2] time_solve[i, 3] time_solve[i, 4]
end
open("examples/ex2_time.txt", "w") do io
    writedlm(io, [nvec time_solve])
end