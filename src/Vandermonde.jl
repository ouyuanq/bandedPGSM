# computing the coefficients of monomial expansion of a function
import Polynomials.coeffs as pcoeffs
import Polynomials.fit as pfit

function poly_coeffs(::Type{T}, f::Function, n::Integer) where T
    # compute the first n+1 monomial coefficients of function f through solving a Vandermonde matrix
    
    # x = chebpts(10n+1, T)
    x = range(-1, 1, length = 10n+1)

    # # Vandermonde matrix
    # V = Matrix{T}(undef, length(x), n+1)
    # V[:, 1] .= 1
    # for i = 2:n+1
    #     broadcast!(*, view(V, :, i), view(V, :, i-1), x)
    # end

    # rhs = f.(x)

    # return ldiv!(Vector{T}(undef, n+1), qr!(V), rhs)

    return pcoeffs(pfit(x, f.(x), n))
end