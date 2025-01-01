## examples of ODEs
import SpecialFunctions: airy

function airy(::Type{T}; epsilon=T(1e-9)) where {T}
    # epsilon * u'' - x * u = 0, u(-1) = Ai(-(1/epsilon)^(1/3)), u(1) = Ai((1/epsilon)^(1/3))
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = [0; -one(T)]
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = [epsilon]

    v = [1 -1; 1 1] \ [airyai(-(1 / epsilon)^(1 / 3)); airyai((1 / epsilon)^(1 / 3))]

    lincoeffs, v
end

airy(::Type{T}, n::Integer) where {T} = secondDirichlet_ChebyshevW(T, n)

function airy_GSBSPG(::Type{T}; epsilon=T(1e-9)) where {T}
    # epsilon * u'' - x * u = 0, u(-1) = Ai(-(1/epsilon)^(1/3)), u(1) = Ai((1/epsilon)^(1/3))
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = [0; -one(T)]
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = [epsilon]

    v = [1 -1; 1 1] \ [airyai(-(1 / epsilon)^(1 / 3)); airyai((1 / epsilon)^(1 / 3))]

    lincoeffs, v
end

function airy_ultraS(::Type{T}, n::Integer; epsilon=T(1e-9)) where {T}
    # epsilon * u'' - x * u = 0, u(-1) = Ai(-(1/epsilon)^(1/3)), u(1) = Ai((1/epsilon)^(1/3))
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = [0; -one(T)]
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = [epsilon]

    bc = ones(T, 2, n)
    bc[1, 2:2:end] .= -1
    bcvals = [airyai(-(1 / epsilon)^(1 / 3)); airyai((1 / epsilon)^(1 / 3))]

    lincoeffs, bc, bcvals
end

function expx2(::Type{T}) where {T}
    # u'' - (1 + sin(x)) * u' + exp(x) * u = f, u(-1) = 1, u(1) = 1, where f is chosen s.t. u = exp((x^2 - 1) / 2)
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = coeffs(exp, T)
    lincoeffs[2] = coeffs(x -> - 1 - sin(x), T)
    lincoeffs[3] = ones(T, 1)

    v = [1 -1; 1 1] \ [1; 1]

    lincoeffs, v
end

expx2(::Type{T}, n::Integer) where {T} = secondDirichlet_ChebyshevW(T, n)

function expx2_GSBSPG(::Type{T}) where {T}
    # u'' - (1 + sin(x)) * u' + exp(x) * u = f, u(-1) = 1, u(1) = 1, where f is chosen s.t. u = exp((x^2 - 1) / 2)
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = exptaylor(T, 17)
    lincoeffs[2] = sintaylor(T, 17)
    lincoeffs[2][1] += 1
    lmul!(-1, lincoeffs[2])
    lincoeffs[3] = ones(T, 1)

    v = [1 -1; 1 1] \ [1; 1]

    lincoeffs, v
end

function exptaylor(::Type{T}, n::Integer) where {T}
    # compute the Taylor expansion of exp(x) at 0 till degree n monomials (x^n)
    a = Vector{T}(undef, n + 1)
    a[1] = 1
    for i = 1:n
        a[i+1] = 1 / factorial(i)
    end

    a
end

function sintaylor(::Type{T}, n::Integer) where {T}
    # compute the Taylor expansion of sin(x) at 0 till degree n monomials (x^n)
    a = zeros(T, n + 1)
    for i = 1:2:n
        a[i+1] = 1 / factorial(i)
    end
    lmul!(-1, view(a, 4:4:n+1))

    a
end

function exptan(::Type{T}; a = T(5e4)) where {T}
    # (a*x^2 + 1) * u' + u = 0, u(-1) = 1
    lincoeffs = Vector{Vector{T}}(undef, 2)

    lincoeffs[1] = ones(T, 1)
    lincoeffs[2] = [a/2+1; zero(T); a/2]

    v = ones(T, 1)

    lincoeffs, v
end

exptan(::Type{T}, n::Integer) where {T} = leftDirichlet_ChebyshevW(T, n)

function exptan_GSBSPG(::Type{T}; a = T(5e4)) where {T}
    # (a*x^2 + 1) * u' + u = 0, u(-1) = 1
    lincoeffs = Vector{Vector{T}}(undef, 2)

    lincoeffs[1] = ones(T, 1)
    lincoeffs[2] = [one(T); zero(T); a]

    v = ones(T, 1)

    lincoeffs, v
end

function tenth(::Type{T}) where {T}
    # u^(10) + cosh(x)*u^(8) + x^2*u^(6) + x^4*u^(4) + cos(x)*u^(2) + x^2*u = 0
    # u(-1)=u(1)=0, u'(-1)=u'(1)=1, u^(k)(-1)=u^(k)(1)=0, 2<=k<=4

    lincoeffs = Vector{Vector{T}}(undef, 11)
    for j = 2:2:10
        lincoeffs[j] = zeros(T, 0)
    end
    lincoeffs[1] = [T(0.5); 0; T(0.5)]
    lincoeffs[3] = coeffs(cos, T)
    lincoeffs[5] = [T(0.375); 0; T(0.5); 0; T(0.125)]
    lincoeffs[7] = [T(0.5); 0; T(0.5)]
    lincoeffs[9] = coeffs(cosh, T)
    lincoeffs[11] = [T(1)]

    C = ones(T, 10, 10)
    C[1, 2:2:end] .= -1
    C[4, :] .= (0:9).^2
    C[3, 1:2:end] .= -C[4, 1:2:end]
    C[3, 2:2:end] .= C[4, 2:2:end]
    C[6, :] .= ldiv!(3, broadcast(*, view(C, 4, :), (0:9).^2 .- 1))
    C[5, 1:2:end] .= C[6, 1:2:end]
    C[5, 2:2:end] .= -C[6, 2:2:end]
    C[8, :] .= ldiv!(5, broadcast(*, view(C, 6, :), (0:9).^2 .- 4))
    C[7, 1:2:end] .= -C[8, 1:2:end]
    C[7, 2:2:end] .= C[8, 2:2:end]
    C[10, :] .= ldiv!(7, broadcast(*, view(C, 8, :), (0:9).^2 .- 9))
    C[9, 1:2:end] .= C[10, 1:2:end]
    C[9, 2:2:end] .= -C[10, 2:2:end]

    v = C \ [zeros(T, 2); ones(T, 2); zeros(T, 6)]

    lincoeffs, v
end

tenth(::Type{T}, n::Integer) where {T} = tenthorder(T, n)

function tenth_ultraS(::Type{T}, n::Integer) where {T}
    # u^(10) + cosh(x)*u^(8) + x^2*u^(6) + x^4*u^(4) + cos(x)*u^(2) + x^2*u = 0
    # u(-1)=u(1)=0, u'(-1)=u'(1)=1, u^(k)(-1)=u^(k)(1)=0, 2<=k<=4
    lincoeffs = Vector{Vector{T}}(undef, 11)
    for j = 2:2:10
        lincoeffs[j] = zeros(T, 0)
    end
    lincoeffs[1] = [T(0.5); 0; T(0.5)]
    lincoeffs[3] = coeffs(cos, T)
    lincoeffs[5] = [T(0.375); 0; T(0.5); 0; T(0.125)]
    lincoeffs[7] = [T(0.5); 0; T(0.5)]
    lincoeffs[9] = coeffs(cosh, T)
    lincoeffs[11] = [T(1)]

    bc = ones(T, 10, n)
    bc[1, 2:2:end] .= -1
    broadcast!(*, view(bc, 4, :), 0:n-1, 0:n-1)
    copyto!(view(bc, 3, :), view(bc, 4, :))
    lmul!(-1, view(bc, 3, 1:2:n))
    ldiv!(3, broadcast!(*, view(bc, 6, :), view(bc, 4, :), (0:n-1).^2 .- 1))
    copyto!(view(bc, 5, :), view(bc, 6, :))
    lmul!(-1, view(bc, 5, 2:2:n))
    ldiv!(5, broadcast!(*, view(bc, 8, :), view(bc, 6, :), (0:n-1).^2 .- 4))
    copyto!(view(bc, 7, :), view(bc, 8, :))
    lmul!(-1, view(bc, 7, 1:2:n))
    ldiv!(7, broadcast!(*, view(bc, 10, :), view(bc, 8, :), (0:n-1).^2 .- 9))
    copyto!(view(bc, 9, :), view(bc, 10, :))
    lmul!(-1, view(bc, 9, 2:2:n))

    bcvals = zeros(T, 10)
    bcvals[3:4] .= 1

    lincoeffs, bc, bcvals
end