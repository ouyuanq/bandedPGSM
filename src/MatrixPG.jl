## solving ODEs with banded Petrov-Galerkin methods

function bandedPGsolve(lincoeffs::Vector{Vector{T}}, v::AbstractVector{T}, R::BandedMatrix{T}, Q::BandedMatrix{T}, Omega::Diagonal{T}, f) where {T}
    # solve the linear system derived from banded Petrov-Galerkin methods
    N = length(lincoeffs)-1
    @assert size(R, 1) == size(Q, 1) && R.l == Q.l == N "Incompatible of trial and test spaces"
    @assert size(Omega, 1) == size(Q, 1) "Incompatible of inner product matrix and test spaces"
    n = size(R, 1)

    # construction for linear system corresponding to main equation
    A, Au = Lmatrix_Chebyshev(lincoeffs, n)

    # construction for rhs (assume that f is expressed in Chebyshev basis)
    if isa(f, Function)
        b = fastconv!(coeffs(f, T), 0, N)
    elseif isa(f, AbstractVector)
        b = fastconv!(f, 0, N)
    end
    if length(b) > n
        deleteat!(b, n+1:length(b))
    else
        append!(b, Zeros(T, n-length(b)))
    end
    # subtract the interpolation
    mul!(view(b, 1:colrange(A, length(v))[end]), view(A, 1:colrange(A, length(v))[end], 1:length(v)), v, -1, true)
    OQ = Dlmul!(Omega, copy(Q))  # merge Q^T * Omega = (Omega * Q)^T
    Wtlmul!(OQ, b)

    # full PG matrix, i.e., Q^T * Omega * A * R
    AW = Wrmul!(A, R, Au)
    WtAW = Wtlmul!(OQ, AW, Au)
    AB, ipiv = LinearAlgebra.LAPACK.gbtrf!(A.l, A.u-A.l, n-N, WtAW.data)
    x = LinearAlgebra.LAPACK.gbtrs!('N', A.l, A.u-A.l, n-N, AB, ipiv, b)

    # convert back to Chebyshev series x = R * x + v
    Wlmul!(R, x)
    if length(x) > length(v)
        axpy!(true, v, view(x, 1:length(v)))
    else
        axpy!(true, view(v, 1:length(x)), x)
    end

    x
end

function bandedPGmatrix_Chebyshev(lincoeffs::Vector{Vector{T}}, R::BandedMatrix{T}, Q::BandedMatrix{T}, Omega::Diagonal{T}) where {T}
    # construct the matrix of banded Petrov-Galerkin methods related to Chebyshev basis
    # lincoeffs are coefficients of linear ODEs, R and Q are matrices related to recombination of trial and test basis and Omega is the inner product matrix
    @assert size(R, 1) == size(Q, 1) "Incompatible of trial and test spaces"
    @assert size(Omega, 1) == size(Q, 1) "Incompatible of inner product matrix and test spaces"
    n = size(R, 1)

    # linear operator
    A, Au = Lmatrix_Chebyshev(lincoeffs, n)
    # full PG matrix, i.e., Q^T * Omega * A * R
    AW = Wrmul!(A, R, Au)
    OQ = Dlmul!(Omega, copy(Q))  # merge Q^T * Omega = (Omega * Q)^T
    WtAW = Wtlmul!(OQ, AW, Au)

    WtAW
end

function Lmatrix_Chebyshev(lincoeffs::Vector{Vector{T}}, n::Integer) where T
    # A = matrix(lincoeffs, n) returns a matrix A of the US representation of the differential operator
    # A*u = (lincoeffs[N+1]*D^{N} + ... + lincoeffs[2]*D^{1} + lincoeffs[1]*D^{0})*u where N is the differential order

    N = length(lincoeffs) - 1
    lencoeffs = length.(lincoeffs)

    # the bandwidths of A
    Al = maximum(lencoeffs .- (0:N)) - 1
    Au = maximum(lencoeffs .+ (2*N:-1:N)) - 1

    A = BandedMatrix{T}(undef, (n, n), (Al+N, Au+Al+2*N))  # extra lower bandwidth for right multiplication of transformation operator and extra upper bandwidth for factorization
    fill!(A.data, zero(T))
    Mu = maximum(lencoeffs)-1  # bandwidths for multiplication operator
    M = BandedMatrix{T}(undef, (n+2*N, n), (Mu, Mu))

    # construction for each term
    Alnow, Aunow = -1, -1
    for i = 0:N
        if !isempty(lincoeffs[i+1])
            # nonzero terms
            fill!(M.data, zero(T))  # empty the storage
            multiplication!(M, lincoeffs[i+1], i)
            # effective elements of multiplication operator
            Mdata_i = view(M.data, Mu+2-lencoeffs[i+1]:Mu+lencoeffs[i+1], 1:n-i)
            if i > 0
                # differential operator application
                broadcast!(*, Mdata_i, Mdata_i, ((2^(i-1)*factorial(i-1)) .* (i:n-1))')
            end
            # add M * D to A
            axpy!(true, Mdata_i, view(A.data, A.u+2-i-lencoeffs[i+1]:A.u-i+lencoeffs[i+1], i+1:n))
            Alnow, Aunow = max(Alnow, lencoeffs[i+1]-1-i), max(Aunow, lencoeffs[i+1]-1+i)
        end
        if i < N
            Slmul!(A, i, i+1, Alnow, Aunow)
            Aunow += 2
        end
    end

    A, Au
end