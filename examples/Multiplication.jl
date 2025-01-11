## test two ways for construction of multiplication operator, namely recursion on coefficients of function and order reduction
include("../src/BandedPG.jl")
using BenchmarkTools, DelimitedFiles, Printf

# first example: multiplication of variable degree polynomials in C^{(10)} basis
N = 10000
nvec = 2 .^ (5:11)
time_C10 = zeros(length(nvec), 2)
T = Float64
for i in eachindex(nvec)
    @printf "No.%d iteration \n" i
    fn = nvec[i]
    f_coeffs = rand(T, fn)
    f = Fun(Chebyshev(), f_coeffs)

    M = view(Multiplication(f, Ultraspherical(10)), 1:N, 1:N)

    ben = @benchmark BandedMatrix($(M))
    time_C10[i, 1] = minimum(ben).time / 1e9

    ben = @benchmark multiplication($(f_coeffs), 10, $(N), $(N))
    time_C10[i, 2] = minimum(ben).time / 1e9
end
open("examples/Multiplication_variable_C10.txt", "w") do io
    writedlm(io, [nvec time_C10])
end

# second example: multiplication of low-degree polynomials in C^{(2)} basis
nvec = 2 .^ (5:21)
time_C2 = zeros(length(nvec), 2)
T = Float64
f_coeffs = rand(T, 20)
f = Fun(Chebyshev(), f_coeffs)
for i in eachindex(nvec)
    @printf "No.%d iteration \n" i
    n = nvec[i]

    M = view(Multiplication(f, Ultraspherical(2)), 1:n, 1:n)

    ben = @benchmark BandedMatrix($(M))
    time_C2[i, 1] = minimum(ben).time / 1e9

    ben = @benchmark multiplication($(f_coeffs), 2, $(n), $(n))
    time_C2[i, 2] = minimum(ben).time / 1e9
end
open("examples/Multiplication_low_C2.txt", "w") do io
    writedlm(io, [nvec time_C2])
end