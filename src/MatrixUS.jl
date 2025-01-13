## solving ODEs with ultraspherical spectral methods

function ultraSsolve(lincoeffs::Vector{Vector{T}}, bc::AbstractMatrix{T}, bcvals::AbstractVector{T}, f) where {T}
    # solve the linear system derived from ultraspherical spectral methods
    N = length(lincoeffs)-1
    @assert size(bc, 1) == size(bcvals, 1) == N "Incompatible of boundary conditions"

    n = size(bc, 2)

    # construction for linear system
    A = matrix_ultraS(lincoeffs, bc)

    # construction for rhs (assume that f is expressed in Chebyshev basis)
    if isa(f, Function)
        b = fastconv!(coeffs(f, T), 0, N)
    elseif isa(f, AbstractVector)
        b = fastconv!(f, 0, N)
    end
    prepend!(b, bcvals)
    if length(b) > n
        deleteat!(b, n+1:length(b))
    else
        append!(b, Zeros(T, n-length(b)))
    end

    x = A \ b

    x
end

# original construction
function matrix_ultraS(lincoeffs::Vector{Vector{T}}, bc::AbstractMatrix{T}) where T
    # A = matrix(lincoeffs, n) returns a matrix A of the US representation of the differential operator
    # A*u = (lincoeffs[N+1]*D^{N} + ... + lincoeffs[2]*D^{1} + lincoeffs[1]*D^{0})*u where N is the differential order

    N = length(lincoeffs) - 1
    n = size(bc, 2)
    lencoeffs = length.(lincoeffs)

    # the bandwidths of A
    Al = maximum(lencoeffs .- (0:N)) - 1
    Au = maximum(lencoeffs .+ (2*N:-1:N)) - 1

    L = BandedMatrix(Zeros(T, n, n), (Al+N, Au-N))  # banded part of final matrix 
    A = _BandedMatrix(L.data, Base.OneTo(n-N), Al, Au)

    # construction for each term
    Alnow, Aunow = -2*N, -2*N
    for i = 0:N
        if !isempty(lincoeffs[i+1])
            # multiplication operator
            M = multmat(lincoeffs[i+1], i, n+N-2*i, n-i)

            # left multiplication of differentiation operator
            if i > 0
                broadcast!(*, M.data, M.data, ((2^(i-1)*factorial(i-1)) .* (i:n-1))')
            end

            # right multiplication of conversion operator
            Ci = convertmat(T, n-N, n+N-2*i, i, N)
            Ai = _BandedMatrix(view(A.data, :, i+1:n), Base.OneTo(n-N), Al+i, Au-i)
            mybandedmul!(Ai, Ci, M, true, true)
        end
    end

    # assign low-rank part
    for i in axes(bc, 1)
        rowi = rowrange(L, i)
        L[i, rowi] = view(bc, i, rowi)
    end

    Lbc = AlmostBandedMatrix(L, ApplyMatrix(*, Matrix{T}(I, n, size(bc, 1)), bc))

    Lbc
end

## US + basis recombination
function matrix_ultraS_br(lincoeffs::Vector{Vector{T}}, n::Integer, br::Function, v::Vector{T}) where T
    # the same as matrix_ultraS except that boundary conditions are treated by basis recombination through function br

    N = length(lincoeffs) - 1
    lencoeffs = length.(lincoeffs)

    # the bandwidths of A
    Al = maximum(lencoeffs .- (0:N)) - 1
    Au = maximum(lencoeffs .+ (2*N:-1:N)) - 1

    A = BandedMatrix(Zeros(T, n-N, n), (Al+N, Au))  # banded part of final matrix 
    Mu = maximum(lencoeffs)-1  # bandwidths for multiplication operator
    M = BandedMatrix(Zeros(T, n+N, n), (Mu, Mu))

    # construction for each term
    Alnow, Aunow = -2*N, -2*N
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

    # construct the rhs
    b = zeros(T, n-N)
    mul!(b, view(A, 1:length(b), 1:length(v)), v, -1, true)

    # left multiplication of matrix related to basis recombination
    W = br(T, n)
    AW = Wrmul!(A, W, Au)

    AW, W, b
end