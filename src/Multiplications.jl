## fast construction for multiplication operator in ultraspherical spectral method

## in-place construction of multiplication operator
function multiplication!(M::BandedMatrix{T}, a::AbstractVector, lambda::Integer) where {T}
    # calculate M_{lambda}[a] by order reduction, i.e., M_{lambda}[a] = S_{lambda-1}...S_1 M_{1}[a] S_1^{-1}...S_{lambda-1}^{-1}, where a is the Chebyshev coefficients of the function
    # for lambda = 0 or 1, the multiplication operator is structured and can be constructed directly
    # the output is a BandedMatrix, the same as that of ApproxFun
    # the elements of M are assumed to be zeros

    Mdata, Mu = M.data, M.u
    la = length(a)
    @assert la <= M.l+1 && la <= M.u+1 "The bandwidths are too small"
    if la == 1
        # constant term
        Mdata[Mu+1, :] .= a[1]
        return M
    end

    a2 = zeros(T, 2*la-1)  # auxiliary vector for fast assignment of Toeplitz part
    axpy!(T(0.5), a, view(a2, la:2*la-1))
    axpy!(T(0.5), a, view(a2, la:-1:1))

    # Toeplitz part
    for i = axes(Mdata, 2)
        Mdata[Mu+2-la:Mu+la, i] .= a2
    end

    if lambda == 0
        # Hankel part + rank-1 part
        for j = 1:min(la - 1, size(M, 2))
            for i = 2:min(la - j + 1, size(M, 1))
                M[i, j] += a[i+j-1] / 2
            end
        end
    elseif lambda >= 1
        # Hankel part
        for j = 1:min(la - 2, size(M, 2))
            for i = 1:min(la - j - 1, size(M, 1))
                M[i, j] -= a[i+j+1] / 2
            end
        end

        if lambda > 1
            # order reduction

            # Note that S_{mu} = (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2)) * diagonal([mu/mu mu/(mu+1) mu/(mu+2) ...], 0) and the scaling diagonal matrix should be applied from left and right at the same first
            for i = 1:lambda-1
                # scaling
                Sscale!(M, i, la - 1)
                # right division of (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2))
                Srdiv!(M, la - 1)
                # left multiplication of (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2))
                Slmul!(M, la - 1)
            end
        end
    else
        error("lambda must be a positive number")
    end

    M
end

# construct a multiplication matrix with fixed dimension
function multiplication(a::AbstractVector, lambda::Integer, m::Integer, n::Integer)
    la = length(a)
    T = eltype(a)
    M = BandedMatrix{T}(undef, m, n, la, la)
    multiplication!(M, a, lambda)
end