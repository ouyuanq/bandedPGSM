# basic functions and matrices used in Shenfun
# see Mortensen, Mikael. "A Generic and Strictly Banded Spectral Petrov–Galerkin Method for Differential Equations with Polynomial Coefficients." SIAM Journal on Scientific Computing 45.1 (2023): A123-A146.

################################  Chebyshev T case ########################################
# the scale of Chebyshev T polynomials from diagonal Jacobi polynomials, i.e., T_n = g_n^{-1/2, -1/2} * P_n^{-1/2, -1/2}
function Chebyshev_scale(::Type{T}, n::Integer) where {T}
    # the scale factor g used for Chebyshev T polynomials
    g = Vector{T}(undef, n)
    gi = one(T)
    for i in eachindex(g)
        g[i] = gi
        gi *= 2 * i / (2 * i - 1)
    end

    g
end

# the normalization factor of Chebyshev T polynomials, i.e., (d(T_n)^k / dx^k, d(T_n)^k / dx^k)_{w^{-1/2, -1/2}} = h_n^k
function Chebyshev_nfactor(::Type{T}, n::Integer, k::Integer) where {T}
    # the normalization factor of Chebyshev T polynomials
    @assert k >= 0 "The differential order cannot be negative"

    hk = Vector{T}(undef, n)
    hk .= pi / 2
    if k == 0
        hk[1] *= 2
    else
        for j = -k+1:k-1
            broadcast!(*, hk, hk, j:j+n-1)
        end
        broadcast!(*, hk, hk, 0:n-1)
    end

    hk
end

# matrices related to recurrence relation of multiplying x, i.e., x^q Q = (A^q)^T Q
function Chebyshev_coreA(::Type{T}, m::Integer, n::Integer, k::Integer) where {T}
    # the innermost A for x*(d(T_n)/dx^k) where k is the differential order
    if k == 0
        # the normal recurrence relation for g_{n}^{-1/2, -1/2} \ P_n^{-1/2, -1/2}
        A = BandedMatrix{T}((1 => Fill(1, min(m, n - 1)), -1 => Fill(1, min(m - 1, n))), (m, n), (1, 1))
        A[2, 1] = 2
        ldiv!(2, A.data)
    else
        # the recurrence relation corresponding to g_{n}^{-1/2, -1/2} \ P_n^{k-1/2, k-1/2}
        A = Ultra_unscaledA(T, m, n, k - 1 / 2)
        gA = Chebyshev_scale(T, max(m, n) + k)
        ldiv!(Diagonal(view(gA, k+1:k+m)), A)
        rmul!(A.data, Diagonal(view(gA, k+1:k+n)))
    end

    A
end

function Chebyshev_A(::Type{T}, m::Integer, n::Integer, k::Integer, q::Integer) where {T}
    # the shifted matrix for x^q*(d(T_n)/dx^k), which is abbreviated as A^(k, q)_(k, k)
    # the dimension of matrix is m-by-n
    @assert q >= 0 "The index q of x^q must be nonnegative"

    if q == 0
        # x^0 = 1, i.e., the shifted identity matrix
        A = BandedMatrix{T}(undef, m, n, 0, 0)
        A.data .= 1
    elseif q == 1
        # no matrix power needed
        A = Chebyshev_coreA(T, m, n, k)
        # divided by normalization factors if necessary
        if k > 0
            nfA = normalfactors(T, max(m, n) + k, k, -1 / 2, -1 / 2)
            ldiv!(Diagonal(view(nfA, k+1:m+k)), A)
            rmul!(A.data, Diagonal(view(nfA, k+1:n+k)))
        end
    else
        mn = max(m, n)
        # first construct the matrix related to multiplying x
        A = Chebyshev_coreA(T, mn, mn, k)
        # divided by normalization factors if necessary
        if k > 0
            nfA = normalfactors(T, mn + k, k, -1 / 2, -1 / 2)
            ldiv!(Diagonal(view(nfA, k+1:mn+k)), A)
            rmul!(A.data, Diagonal(view(nfA, k+1:mn+k)))
        end
        # Note that the above A is already a shifted matrix, i.e., the zeros in first k rows and k columns of original A are omitted

        # q-th matrix power (q >= 1)
        Aq = A * A
        for j = 3:q
            Aq = Aq * A
        end
        A = Aq

        # truncation
        A = _BandedMatrix(view(A.data, :, 1:n), Base.OneTo(m), q, q)
    end

    A
end

# matrices related to recurrence relation of derivative, i.e., d(Q)^(k-l)/dx^(k-l) = (B^l)^T d(Q)^k/dx^k
function Chebyshev_coreB(::Type{T}, m::Integer, n::Integer) where {T}
    # the innermost B for Q = B^T dQ/dx
    B = BandedMatrix{T}(undef, m, n, 1, 1)
    Bdata = B.data
    Bdata[2, :] .= 0
    Bdata[1, :] .= -1
    Bdata[3, :] .= 1
    B[2, 1] = 2
    B[1, 2] = 0
    Bdata1 = view(Bdata, 1, 3:n)
    broadcast!(/, Bdata1, Bdata1, 2:2:2*n-4)
    Bdata3 = view(Bdata, 3, 1:n)
    broadcast!(/, Bdata3, Bdata3, 2:2:2*n)

    B
end

function Chebyshev_B(::Type{T}, m::Integer, n::Integer, k::Integer, l::Integer) where {T}
    # the shifted matrix for B^l, which is abbreviated as B^(l)_(k)
    # the dimension of matrix is m-by-n
    # Note that only rows are shifted upward
    @assert l >= 0 "The index l of B^l must be nonnegative"

    if l == 0
        # the shifted identity matrix
        B = BandedMatrix{T}((k => Fill(1, min(m, n - k)),), (m, n), (-k, k))
    elseif l == 1
        # no matrix power needed
        B = Chebyshev_coreB(T, m + k, n)
        B = _BandedMatrix(view(B.data, :, 1:n), Base.OneTo(m), l - k, k + l)
    else
        # first construct the matrix related to derivative
        mkn = max(m + k, n)
        B = Chebyshev_coreB(T, mkn, mkn)

        # l-th matrix power (l >= 1)
        Bl = B * B
        for j = 3:l
            Bl = Bl * B
        end
        B = Bl

        # shifted matrix
        B = _BandedMatrix(view(B.data, :, 1:n), Base.OneTo(m), l - k, k + l)
    end

    B
end

# building blocks for Chebyshev T polynomials
function Chebyshev_L(::Type{T}, m::Integer, n::Integer, k::Integer, q::Integer, l::Integer) where {T}
    # the matrix L^{(k,q,l)} = A_{k, k}^{k, q} B_{k}^{l} where A_{k, k}^{k, q} is a m-by-n matrix and B_{k}^{l} is a n-by-n matrix
    L = Chebyshev_A(T, m, n, k, q) * Chebyshev_B(T, n, n, k, l)

    L
end

# GSBSPG discretization of differential operators by Chebyshev T polynomials
function GSBSPG_Chebyshev(::Type{T}, lincoeffs::Vector{Vector{T}}, K::BandedMatrix{T}, bc::AbstractVector{T}, fc::AbstractVector{T}) where {T}
    # construct the GSBSPG discretization of operator lincoeffs[1]*u + lincoeffs[2]*u' + ... + lincoeffs[N+1]*u^{N} by Chebyshev T polynomials
    # the coefficients are store as monomials, i.e., lincoeffs[1] = [0; 0; 1] respresenting x^2
    # bc is a low degree Chebyshev polynomial that satisfies boundary conditions and fc is a vector containing inner products of rhs and test functions

    # differential order
    N = length(lincoeffs) - 1

    # dimension
    n, m = size(K)

    # final matrix
    ql = maximum(length.(lincoeffs) .+ (N-1:-1:-1))
    # L should be a (ql-N, ql+N) banded matrix, but we construct it with extra lower diagonals and upper diagonals for multiplication with K and factorization respectively
    L = BandedMatrix{T}(undef, m, n, ql, 2 * ql + N)
    L.data .= 0

    for l in eachindex(lincoeffs)
        # B matrix, i.e., B_{N}^{N-l+1}, whose bandwidths are (N-l+1-N, N-l+1+N)
        B = Chebyshev_B(T, n - l + 1, n, N, N - l + 1)

        # A matrix for each monomials
        if !iszero(lincoeffs[l])
            dcoeffs = length(lincoeffs[l]) - 1 # degree of coefficient
            # A_{k, k}^{k, dcoeffs} whose bandwidths are (dcoeffs, dcoeffs)
            sizeA = max(m, n - l + 1) + dcoeffs
            Aq = Chebyshev_A(T, sizeA, sizeA, N, 0)
            A = Chebyshev_A(T, sizeA, sizeA, N, 1)

            for q in eachindex(lincoeffs[l])
                if q > 1
                    # matrix power
                    Aq = Aq * A
                end

                if lincoeffs[l][q] != 0
                    # nonzero monomials (A_{k, k}^{k, q})
                    Akq = _BandedMatrix(view(Aq.data, :, 1:n-l+1), Base.OneTo(m), q - 1, q - 1)
                    mybandedmul!(L, Akq, B, lincoeffs[l][q], true)
                end
            end
        end
    end

    # add bc to rhs
    bcrange = colrange(L, length(bc))
    mul!(view(fc, bcrange), view(L, bcrange, 1:length(bc)), bc, -1, true)

    # equip homogeneous boundary conditions to L (with ql+N nonzero upper diagonals)
    L = Wrmul!(L, K, ql + N)

    L, fc
end

# construction and solution of GSBSPG
function GSBSPG_Chebyshev_solve(::Type{T}, lincoeffs::Vector{Vector{T}}, K::BandedMatrix{T}, bc::AbstractVector{T}, fc::AbstractVector{T}) where {T}
    # construct the GSBSPG discretization (recursion) of operator lincoeffs[1]*u + lincoeffs[2]*u' + ... + lincoeffs[N+1]*u^{N} by Chebyshev T polynomials and solve by banded solver

    n = size(K, 1)
    if length(fc) >= n
        f = fc[1:n]
    else
        f = copy(fc)
        append!(f, Zeros(T, n - length(fc)))
    end

    # matrix
    L, u = GSBSPG_Chebyshev(T, lincoeffs, K, bc, f)

    # standard solver (note that L is constructed with extra upper diagonals for factorization)
    ldiv!(lu!(L), view(u, 1:size(K, 2)))
    
    # convert back to Chebyshev series x = K * x + bc
    Wlmul!(K, u)
    if length(u) > length(bc)
        axpy!(true, bc, view(u, 1:length(bc)))
    else
        axpy!(true, view(bc, 1:length(u)), u)
    end

    u
end

# rhs for Chebyshev T polynomials
function Chebyshev_rhs(::Type{T}, f, m::Integer, n::Integer, k::Integer) where {T}
    # compute the rhs in GSBSPG, i.e., f = B_{k}^{k} fhat where fhat = (f, Qtilde) and Qtilde = h\Q is scaled test function
    # Note that m, n are dimensions of coefficient matrix

    # compute the nomal Chebyshev coefficients
    if isa(f, Function)
        # function handle
        fc = coeffs(f, n + k - 1, T)
    else
        # Chebyshev coefficients
        fc = zeros(T, n + k)
        if length(f) < n + k
            copyto!(fc, f)
        else
            copyto!(fc, view(f, 1:n+k))
        end
    end
    # multiply factors <T_n, T_n>
    broadcast!(*, fc, fc, pi / 2)
    fc[1] *= 2

    # divide the normalization factors h
    broadcast!(/, fc, fc, Chebyshev_nfactor(T, n + k, 0))

    # multiplying B_{k}^{k}
    fc = Chebyshev_B(T, m, n + k, k, k) * fc

    fc
end


################################  Legendre case ########################################
# the scale of Legendre polynomials from diagonal Jacobi polynomials, i.e., L_n = g_n^{0, 0} * P_n^{0, 0}
function Legendre_scale(::Type{T}, n::Integer) where {T}
    # the scale factor g used for Legendre polynomials

    g = Ones(T, n)

    g
end

# the normalization factor of Legendre polynomials, i.e., (d(L_n)^k / dx^k, d(L_n)^k / dx^k)_{w^{0, 0}} = h_n^k
function Legendre_nfactor(::Type{T}, n::Integer, k::Integer) where {T}
    # the normalization factor of Legendre polynomials
    @assert k >= 0 "The differential order cannot be negative"

    hk = Vector{T}(undef, n)
    hk .= 2
    broadcast!(/, hk, hk, 1:2:2*n-1)
    for j = 1-k:k
        broadcast!(*, hk, hk, j:j+n-1)
    end

    hk
end

# matrices related to recurrence relation of multiplying x, i.e., x^q Q = (A^q)^T Q
function Legendre_coreA(::Type{T}, m::Integer, n::Integer, k::Integer) where {T}
    # the innermost A for x*(d(L_n)/dx^k) where k is the differential order
    if k == 0
        # the normal recurrence relation for g_{n}^{0, 0} \ P_n^{0, 0}
        A = BandedMatrix{T}(undef, m, n, 1, 1)
        A.data[2, :] .= 0
        A.data[1, :] .= 0:n-1
        A.data[3, :] .= 1:n
        rdiv!(A.data, Diagonal(1:2:2*n-1))
    else
        # the recurrence relation corresponding to g_{n}^{0, 0} \ P_n^{k, k}
        A = Ultra_unscaledA(T, m, n, k)

        ## NOTE: gA is constant 1 and we omit the scaling here
        # gA = Legendre_scale(T, max(m, n) + k)
        # ldiv!(Diagonal(view(gA, k+1:k+m)), A)
        # rmul!(A.data, Diagonal(view(gA, k+1:k+n)))
    end

    A
end

function Legendre_A(::Type{T}, m::Integer, n::Integer, k::Integer, q::Integer) where {T}
    # the shifted matrix for x^q*(d(L_n)/dx^k), which is abbreviated as A^(k, q)_(k, k)
    # the dimension of matrix is m-by-n
    @assert q >= 0 "The index q of x^q must be nonnegative"

    if q == 0
        # x^0 = 1, i.e., the shifted identity matrix
        A = BandedMatrix{T}(undef, m, n, 0, 0)
        A.data .= 1
    elseif q == 1
        # no matrix power needed
        A = Legendre_coreA(T, m, n, k)
        # divided by normalization factors if necessary
        if k > 0
            nfA = normalfactors(T, max(m, n) + k, k, 0, 0)
            ldiv!(Diagonal(view(nfA, k+1:m+k)), A)
            rmul!(A.data, Diagonal(view(nfA, k+1:n+k)))
        end
    else
        mn = max(m, n)
        # first construct the matrix related to multiplying x
        A = Legendre_coreA(T, mn, mn, k)
        # divided by normalization factors if necessary
        if k > 0
            nfA = normalfactors(T, mn + k, k, 0, 0)
            ldiv!(Diagonal(view(nfA, k+1:mn+k)), A)
            rmul!(A.data, Diagonal(view(nfA, k+1:mn+k)))
        end
        # Note that the above A is already a shifted matrix, i.e., the zeros in first k rows and k columns of original A are omitted

        # q-th matrix power (q >= 1)
        Aq = A * A
        for j = 3:q
            Aq = Aq * A
        end
        A = Aq

        # truncation
        A = _BandedMatrix(view(A.data, :, 1:n), Base.OneTo(m), q, q)
    end

    A
end

# matrices related to recurrence relation of derivative, i.e., d(Q)^(k-l)/dx^(k-l) = (B^l)^T d(Q)^k/dx^k
function Legendre_coreB(::Type{T}, m::Integer, n::Integer) where {T}
    # the innermost B for Q = B^T dQ/dx
    B = BandedMatrix{T}(undef, m, n, 1, 1)
    Bdata = B.data
    Bdata[2, :] .= 0
    Bdata[1, :] .= -1
    Bdata[3, :] .= 1
    B[1, 2] = 0
    Bdata1 = view(Bdata, 1, 3:n)
    broadcast!(/, Bdata1, Bdata1, 5:2:2*n-1)
    Bdata3 = view(Bdata, 3, 1:n)
    broadcast!(/, Bdata3, Bdata3, 1:2:2*n-1)

    B
end

function Legendre_B(::Type{T}, m::Integer, n::Integer, k::Integer, l::Integer) where {T}
    # the shifted matrix for B^l, which is abbreviated as B^(l)_(k)
    # the dimension of matrix is m-by-n
    # Note that only rows are shifted upward
    @assert l >= 0 "The index l of B^l must be nonnegative"

    if l == 0
        # the shifted identity matrix
        B = BandedMatrix{T}((k => Fill(1, min(m, n - k)),), (m, n), (-k, k))
    elseif l == 1
        # no matrix power needed
        B = Legendre_coreB(T, m + k, n)
        B = _BandedMatrix(view(B.data, :, 1:n), Base.OneTo(m), l - k, k + l)
        # B = B[k+1:m+k, 1:n]
    else
        # first construct the matrix related to derivative
        mkn = max(m + k, n)
        B = Legendre_coreB(T, mkn, mkn)

        # l-th matrix power (l >= 2)
        Bl = B * B
        for j = 3:l
            Bl = Bl * B
        end
        B = Bl

        # shifted matrix
        B = _BandedMatrix(view(B.data, :, 1:n), Base.OneTo(m), l - k, k + l)
    end

    B
end

# building blocks for Legendre polynomials
function Legendre_L(::Type{T}, m::Integer, n::Integer, k::Integer, q::Integer, l::Integer) where {T}
    # the matrix L^{(k,q,l)} = A_{k, k}^{k, q} B_{k}^{l} where A_{k, k}^{k, q} is a m-by-n matrix and B_{k}^{l} is a n-by-n matrix
    L = Legendre_A(T, m, n, k, q) * Legendre_B(T, n, n, k, l)

    L
end

# GSBSPG discretization of differential operators by Legendre polynomials
function GSBSPG_Legendre(::Type{T}, lincoeffs::Vector{Vector{T}}, K::BandedMatrix{T}, bc::AbstractVector{T}, fc::AbstractVector{T}) where {T}
    # construct the GSBSPG discretization of operator lincoeffs[1]*u + lincoeffs[2]*u' + ... + lincoeffs[N+1]*u^{N} by Legendre polynomials
    # the coefficients are store as monomials, i.e., lincoeffs[1] = [0; 0; 1] respresenting x^2
    # bc is a low degree Legendre polynomial that satisfies boundary conditions and fc is a vector containing inner products of rhs and test functions

    # differential order
    N = length(lincoeffs) - 1

    # dimension
    n, m = size(K)
    mn = max(m, n)

    # final matrix
    ql = maximum(length.(lincoeffs) .+ (N-1:-1:-1))
    L = BandedMatrix{T}(undef, m, n, ql - N, ql + N)
    L.data .= 0

    for l in eachindex(lincoeffs)
        # B matrix, i.e., B_{k}^{l}
        B = Legendre_B(T, n - l + 1, n, N, N - l + 1)

        # A matrix for each monomials
        if !iszero(lincoeffs[l])
            dcoeffs = length(lincoeffs[l]) - 1 # degree of coefficient
            # A_{k, k}^{k, dcoeffs} whose bandwidths are (dcoeffs, dcoeffs)
            sizeA = max(m, n - l + 1) + dcoeffs
            Aq = Legendre_A(T, sizeA, sizeA, N, 0)
            A = Legendre_A(T, sizeA, sizeA, N, 1)

            for q in eachindex(lincoeffs[l])
                if q > 1
                    # matrix power
                    Aq = Aq * A
                end

                if lincoeffs[l][q] != 0
                    # nonzero monomials (A_{k, k}^{k, q})
                    Akq = _BandedMatrix(view(Aq.data, :, 1:n-l+1), Base.OneTo(m), q - 1, q - 1)
                    mybandedmul!(L, Akq, B, lincoeffs[l][q], true)
                end
            end
        end
    end

    # add bc to rhs
    bcrange = colrange(L, length(bc))
    mul!(view(fc, bcrange), view(L, bcrange, 1:length(bc)), bc, -1, true)

    # equip homogeneous boundary conditions to L
    L = L * K

    L, fc
end


################################  general case ########################################
# normalization factors for k-th derivative of P_n^(alpha, beta)
function normalfactors(::Type{T}, n::Integer, k::Integer, alpha::Number, beta::Number) where {T}
    # compute the normalization factors (n+alpha+beta+1)_k / 2^k from k to k+n-1 where (a)_k = Gamma(a+k)/Gamma(a) is teh Pochhammer symbol

    nf = ones(T, n)
    ab1 = T(alpha + beta + 1)
    for i in eachindex(nf)
        faci = T(i - 1 + ab1)
        for j = 0:k-1
            nf[i] *= faci + j
        end
    end
    ldiv!(2^k, nf)

    nf
end

# recurrence relation for unscaled diagonal Jacobi polynomials, i.e., P_n^{alpha, alpha}
function Ultra_unscaledA(::Type{T}, m::Integer, n::Integer, alpha::Number) where {T}
    # matrix related to recurrence relation for unscaled diagonal Jacobi polynomials
    A = BandedMatrix{T}(undef, m, n, 1, 1)
    Adata = A.data
    Adata[2, :] .= 0
    p, q = T(alpha - 1), T(2 * alpha - 1)
    for i in axes(Adata, 2)
        p, q = p + 1, q + 2
        Adata[1, i] = p / q
        Adata[3, i] = (2 * i * (p + alpha + 1)) / ((q + 1) * q)
    end

    A
end