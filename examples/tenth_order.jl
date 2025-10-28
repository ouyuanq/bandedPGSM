include("../src/BandedPG.jl")
include("Examples.jl")
using BenchmarkTools, DelimitedFiles, Printf

# solving tenth order ODE u^(10) + cosh(x)*u^(8) + x^2*u^(6) + x^4*u^(4) + cos(x)*u^(2) + x^2*u = 0
# u(-1)=u(1)=0, u'(-1)=u'(1)=0, u^(k)(-1)=u^(k)(1)=0, 2<=k<=4

nvec = 2 .^ (5:20)
time_con = zeros(length(nvec), 2)
time_sol = zeros(length(nvec), 2)
T = Float64
tenth_coeffs_br, tenth_v = tenth(T)
N = 10
@printf "construction and solution cost of a tenth order equation:\n"
@printf "    n    US(con)    US(sol)   accUS(con)   accUS(sol)\n"
for i in eachindex(nvec)
    n = nvec[i]

    # original US method
    tenth_coeffs, tenth_bc, tenth_bcvals = tenth_ultraS(T, n)
    L = matrix_ultraS(tenth_coeffs, tenth_bc)
    ben = @benchmark matrix_ultraS($(tenth_coeffs), $(tenth_bc))
    time_con[i, 1] = minimum(ben).time / 1e9

    f_US = [tenth_bcvals; zeros(T, n-10)]
    x1 = L \ f_US
    ben = @benchmark $(L) \ $(f_US)
    time_sol[i, 1] = minimum(ben).time / 1e9

    # accelerated US method
    ben = @benchmark begin
        W = tenth($(T), $(n))
        f_br = zeros($(T), $(n) - $(N))
        matrix_ultraS_br($(tenth_coeffs_br), f_br, W, $(tenth_v))
    end
    time_con[i, 2] = minimum(ben).time / 1e9

    W = tenth(T, n)
    f_br = zeros(T, n - N)
    L_br, f_br = matrix_ultraS_br(tenth_coeffs_br, f_br, W, tenth_v)
    x2 = Wlmul!(W, L_br \ f_br)
    axpy!(true, tenth_v, view(x2, 1:length(tenth_v)))

    ben = @benchmark begin
        x2 = Wlmul!($(W), $(L_br) \ $(f_br))
        axpy!(true, $(tenth_v), view(x2, 1:length($(tenth_v))))
    end
    time_sol[i, 2] = minimum(ben).time / 1e9

    @printf "%5i   %.2e   %.2e    %.2e     %.2e\n" n time_con[i, 1] time_sol[i, 1] time_con[i, 2] time_sol[i, 2]
end
open("examples/tenth_time.txt", "w") do io
    writedlm(io, [nvec time_con time_sol])
end