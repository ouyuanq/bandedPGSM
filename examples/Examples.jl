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

function expcos(::Type{T}) where {T}
    # u''' - cos(x) * u'' + 10 * exp(x) * u = f, u(-1) = u (1) = 1, u'(1) = 1
    lincoeffs = Vector{Vector{T}}(undef, 4)

    lincoeffs[1] = lmul!(10, coeffs(exp, T))
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = lmul!(-1, coeffs(cos, T))
    lincoeffs[4] = ones(T, 1)

    v = [3; 0; 1] ./ 4

    lincoeffs, v
end

expcos(::Type{T}, n::Integer) where {T} = thirdrightNeumann_ChebyshevW(T, n)

function expcos_GSBSPG(::Type{T}) where {T}
    # u''' - cos(x) * u'' + 10 * exp(x) * u = f, u(-1) = u (1) = 1, u'(1) = 1
    lincoeffs = Vector{Vector{T}}(undef, 4)

    lincoeffs[1] = lmul!(10, exptaylor(T, 17))
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = lmul!(-1, costaylor(T, 16))
    lincoeffs[4] = ones(T, 1)

    v = [3; 0; 1] ./ 4

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

## a second order differential equation with composite coefficients
function composite(::Type{T}) where {T}
    # u'' - cos(sin(2*pi*x + 1)) * u = f, u(-1) = u (1) = sin(2)
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = coeffs(x -> -cos(sin(2*pi*x + 1)), T)
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = ones(T, 1)

    v = [sin(T(2)); 0]

    lincoeffs, v
end

composite(::Type{T}, n::Integer) where {T} = secondDirichlet_ChebyshevW(T, n)

function composite_GSBSPG(::Type{T}) where {T}
    # u'' - cos(sin(2*pi*x + 1)) * u = f, u(-1) = u (1) = sin(2)
    lincoeffs = Vector{Vector{T}}(undef, 3)

    lincoeffs[1] = poly_coeffs(T, x -> -cos(sin(2*pi*x + 1)), 92)
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = ones(T, 1)

    v = [sin(T(2)); 0]

    lincoeffs, v
end

## a third order differential equation with integral condition
function expcos(::Type{T}) where {T}
    # u''' - exp(x) * u'' + cos(x) * u = f, u(-1) = cos(sin(1)), u'(-1) = -10*pi*cos(1)*sin(sin(1)), u(0) = cos(sin(1))
    lincoeffs = Vector{Vector{T}}(undef, 4)

    lincoeffs[1] = coeffs(cos, T)
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = lmul!(-1, coeffs(exp, T))
    lincoeffs[4] = ones(T, 1)

    v = [1 -1 1; 0 1 -4; 1 0 -1] \ [cos(sin(1)); -10*pi*cos(1)*sin(sin(1)); cos(sin(1))]

    lincoeffs, T.(v)
end

function expcos_GSBSPG(::Type{T}) where {T}
    # u''' - exp(x) * u'' + cos(x) * u = f, u(-1) = cos(sin(1)), u'(-1) = -10*pi*cos(1)*sin(sin(1)), u(0) = cos(sin(1))
    lincoeffs = Vector{Vector{T}}(undef, 4)

    lincoeffs[1] = costaylor(T, 16)
    lincoeffs[2] = zeros(T, 0)
    lincoeffs[3] = lmul!(-1, exptaylor(T, 16))
    lincoeffs[4] = ones(T, 1)

    v = [1 -1 1; 0 1 -4; 1 0 -1] \ [cos(sin(1)); -10*pi*cos(1)*sin(sin(1)); cos(sin(1))]

    lincoeffs, T.(v)
end

## Taylor expansions
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

function costaylor(::Type{T}, n::Integer) where {T}
    # compute the Taylor expansion of cos(x) at 0 till degree n monomials (x^n)
    a = zeros(T, n + 1)
    for i = 1:2:n+1
        a[i] = 1 / factorial(i-1)
    end
    lmul!(-1, view(a, 3:4:n+1))

    a
end