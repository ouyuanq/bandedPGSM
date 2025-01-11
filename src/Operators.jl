# basic operators for ultraspherical spectral method and boundary treatment

# differential operator
function diffmat(::Type{T}, m::Integer, n::Integer, lambda::Integer) where T
    # banded differential matrix
    if lambda == 0
        D = BandedMatrix{T}(Eye(m, n))
    else
        D = BandedMatrix{T}((lambda => 2^(lambda - 1)*factorial(lambda - 1) .* (lambda:min(m+lambda-1, n-1)),), (m, n), (-lambda, lambda))
    end

    D
end

diffmat(n::Integer, lambda::Integer) = diffmat(Float64, n, n, lambda)  # default type
diffmat(m::Integer, n::Integer, lambda::Integer) = diffmat(Float64, m, n, lambda)

# conversion operator
function convertmat(::Type{T}, m::Integer, n::Integer, lambda::Integer, mu::Integer) where T
    # banded conversion matrix to transfer a C^{λ} series to a C^{μ} series
    # S = S_{mu-1}S_{mu-2}...S_{lambda}, where S_{mu} = (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2)) * diagonal([mu/mu mu/(mu+1) mu/(mu+2) ...], 0)
    if lambda < mu
        # S_{mu-1} with total bandwidth 2*(mu-lambda)
        Su = 2*(mu-lambda)
        S = BandedMatrix{T}((0 => Fill(1, n), 2 => Fill(-1, n-2)), (m, n), (0, Su))
        if mu == 1
            Sl = view(S.data, Su-1:2:Su+1, 2:n)
            broadcast!(/, Sl, Sl, 2)
        else
            Sl = view(S.data, Su-1:2:Su+1, :)
            lmul!(mu-1, Sl)
            broadcast!(/, Sl, Sl, (mu-1:mu+n-2)')
        end
        Sdata = S.data

        # multiplication from left to right
        @inbounds for k = mu-2:-1:lambda
            k2 = 2*(mu-k-1)  # bandwidth before multiplication
            for j = n:-1:3
                j2 = j-2
                for i = j2:-2:j2-k2
                    Sij = inbands_getindex(Sdata, Su, i, j) - inbands_getindex(Sdata, Su, i, j2)
                    inbands_setindex!(Sdata, Su, Sij, i, j)
                end
            end
            # scaling matrix
            if k == 0
                Sl = view(S.data, Su-k2-1:2:Su+1, 2:n)
                broadcast!(/, Sl, Sl, 2)
            else
                Sl = view(S.data, Su-k2-1:2:Su+1, :)
                lmul!(k, Sl)
                broadcast!(/, Sl, Sl, (k:k+n-1)')
            end
        end
    elseif lambda == mu
        S = BandedMatrix{T}(Eye(m, n), (0, 0))
    else
        @error "μ should be no less than λ"
    end
    
    S
end

convertmat(n::Integer, lambda::Integer, mu::Integer) = convertmat(Float64, n, n, lambda, mu)  # default type
convertmat(m::Integer, n::Integer, lambda::Integer, mu::Integer) = convertmat(Float64, m, n, lambda, mu)

## fast conversion operator (banded with only diagonal and one superdiagonal)
function fastconv!(a::AbstractVector{T}, mu::Real, lambda::Real) where {T<:Number}
    # convert C^{mu} series a to C^{lambda} basis
    # The result is stored in the initial vector
    @assert isinteger(lambda - mu) "lambda - mu must be an integer."

    for i = mu:lambda-1
        fastconv1!(a, i)
    end

    a
end

function fastconv1!(a::AbstractVector{T}, lambda::Real) where {T<:Number}
    # apply the conversion matrix S^{lambda} to vector a. The result is stored in the initial vector
    # Note that S_{lambda} = (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2)) * diagonal([lambda/lambda lambda/(lambda+1) lambda/(lambda+2) ...], 0)

    n = length(a)

    # apply the scaling matrix
    if lambda == 0
        ldiv!(2, a)
        a[1] *= 2
    else
        lmul!(lambda, a)
        for i in eachindex(a)
            a[i] = a[i] / (lambda + i - 1)
        end
    end

    # apply the bidigonal matrix
    @inbounds for i = 1:n-2
        a[i] -= a[i+2]
    end

    a
end

function multmat(a::AbstractVector, lambda::Integer, m::Integer, n::Integer)
    fa = Fun(Chebyshev(), a)
    if lambda == 0
        Mf = Multiplication(fa, Chebyshev())
    elseif lambda > 0
        Mf = Multiplication(fa, Ultraspherical(lambda))
    else
        error("Cannot construct a multiplication matrix with negative lambda")
    end

    M = BandedMatrix(view(Mf, 1:m, 1:n))
end