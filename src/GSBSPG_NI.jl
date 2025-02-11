# basic functions and matrices used in Shenfun
# see Mortensen, Mikael. "A Generic and Strictly Banded Spectral Petrov–Galerkin Method for Differential Equations with Polynomial Coefficients." SIAM Journal on Scientific Computing 45.1 (2023): A123-A146.
# all the matrices for the inner product forms are computed by numerical integration

# GSBSPG discretization of differential operators by Chebyshev T polynomials with numerical integration
function GSBSPG_Chebyshev_NI(::Type{T}, lincoeffs::Vector{Vector{T}}, K::BandedMatrix{T}, bc::AbstractVector{T}, fc::AbstractVector{T}) where {T}
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

    # determine the number of points for quadrature
    nq = cld(maximum(length.(lincoeffs) .+ (n+N-2:-1:n-2)) + (m + 2 * N - 1), 2)
    GCx, GCw = gausschebyshevt(nq)
    if eltype(GCx) != T
        GCx = convert(Vector{T}, GCx)
        GCw = convert(Vector{T}, GCw)
    end

    # get the values of linear coefficients on Gauss-Chebyshev points
    xCoeffs = similar(lincoeffs)
    for l in eachindex(lincoeffs)
        if iszero(lincoeffs[l])
            # zero coefficient
            xCoeffs[l] = zeros(T, 0)
        else
            xCoeffs[l] = Vector{T}(undef, nq)
            xCoeffs[l] .= lincoeffs[l][1]
        end
    end
    GCxq = copy(GCx)
    lencoeffs = length.(lincoeffs)
    for q = 2:maximum(lencoeffs)
        # x^q trem
        if q > 2
            broadcast!(*, GCxq, GCxq, GCx)
        end

        # for each coefficient
        for l in eachindex(lincoeffs)
            if q <= lencoeffs[l]
                axpy!(lincoeffs[l][q], GCxq, xCoeffs[l])
            end
        end
    end

    # get values of basis of trial space and its derivatives on Gauss-Chebyshev points
    xTrial = Vector{Vector{Vector{T}}}(undef, N + 1)
    xTrial[1] = Vector{Vector{T}}(undef, 2)
    # the first two Chebyshev T polynomials, 1 and x
    xTrial[1][1], xTrial[1][2] = ones(T, nq), copyto!(GCxq, GCx)
    gChebyshev = Chebyshev_scale(T, max(m, n) + N)
    Trialfactors = Vector{Vector{T}}(undef, N)
    for k = 1:N
        # the k-th derivative of Chebyshev polynomial g_{n+k} \ psi_{n+k}^{k, -1/2, -1/2} * P_{n}^{k-1/2, k-1/2}
        xTrial[k+1] = Vector{Vector{T}}(undef, 2)
        # the first two diagonal Jacobi polynomials, 1 and (k+1/2)*x
        xTrial[k+1][1], xTrial[k+1][2] = ones(T, nq), lmul!(k + 1 / 2, copy(GCx))
        Trialfactors[k] = normalfactors(T, n, k, -1 / 2, -1 / 2)
        broadcast!(*, Trialfactors[k], Trialfactors[k], view(gChebyshev, 1:n))
        Trialfactors[k][1:k] .= 0
        # Note that the constants are not multiplied in order for three-term recurrence relation used later
    end

    # get values of basis of test space on Gauss-Chebyshev points (a window of L.l(ql-N) + L.u(ql+N) + 1 bases at the same time)
    xTest = Vector{Vector{T}}(undef, 2 * ql + 1)
    xTest[1], xTest[2] = ones(T, nq), lmul!(N + 1 / 2, copy(GCx))
    # the factor (1-x^2)^N in test basis
    for t in eachindex(xTest[1])
        x2N = (1 - GCx[t]^2)^N
        xTest[1][t] *= x2N
        xTest[2][t] *= x2N
    end
    # an extra vector for computing recurrence relation and others
    rrtemp = Vector{T}(undef, nq)
    # compute the first 2 * (ql + N) + 1 test basis using three-term recurrence
    for i = 3:length(xTest)
        xTest[i] = copy(xTest[i-2])
        # Note that n = i - 2 and alpha = beta = N - 1/2 are the indexes for recurrence now
        ttr!(xTest[i], xTest[i-1], GCx, i - 2, N - 1/2, rrtemp)
    end
    # normalization factors for test function (not multiplied for recurrence used later)
    Testfactors = view(normalfactors(T, m + N, N, -1 / 2, -1 / 2), N+1:m+N)
    broadcast!(*, Testfactors, Testfactors, view(gChebyshev, N+1:m+N))
    broadcast!(/, Testfactors, Testfactors, view(Chebyshev_nfactor(T, m + N, N), N+1:m+N))

    # Now we construct the matrix of inner product column by column
    xInnerPorduct = Vector{T}(undef, nq)
    xTrialj = Vector{T}(undef, nq)
    for j in axes(L, 2)
        # the values of trial basis in j-th column
        xTrialj .= 0
        for l in eachindex(lincoeffs)
            if !iszero(lincoeffs[l]) && l <= j
                # the (l-1)-th differential term is lincoeffs[l] * d(T_n)^{l-1}/dx^{l-1}
                if l == 1
                    broadcast!(*, rrtemp, xTrial[l][1], xCoeffs[l])
                    axpy!(true, rrtemp, xTrialj) # no scaling
                else
                    broadcast!(*, rrtemp, xTrial[l][1], xCoeffs[l])
                    axpy!(Trialfactors[l-1][j], rrtemp, xTrialj) # with scaling
                end
            end
        end

        # core part: computing the inner product of j-th column
        iTest = 1
        # note that L is a (ql-N, ql+N) banded matrix
        for i = max(1, j - ql - N):min(m, j + ql - N)
            # values of product of functions in the inner product
            # multiply the trial function and the test function
            broadcast!(*, xInnerPorduct, xTest[iTest], xTrialj)

            # numerical integration (Note that factors for test functions are multiplied here)
            L[i, j] = dot(GCw, xInnerPorduct) * Testfactors[i]

            iTest += 1
        end

        # post processing for generating new terms through three-term recurrence relation
        # the trial bases
        if j < n
            # the Chebyshev basis
            xTrial[1][1], xTrial[1][2] = xTrial[1][2], xTrial[1][1]
            broadcast!(*, rrtemp, xTrial[1][1], GCx)
            axpby!(2, rrtemp, -1, xTrial[1][2])

            # the derivatives (diagonal Jacobi basis)
            for k = 1:N
                if k < j
                    xTrial[k+1][1], xTrial[k+1][2] = xTrial[k+1][2], xTrial[k+1][1]
                    # Note that n = j - k and alpha = beta = k - 1/2 are the indexes for recurrence now
                    ttr!(xTrial[k+1][2], xTrial[k+1][1], GCx, j - k, k - 1/2, rrtemp)
                end
            end
        end

        # the test bases
        if j > ql + N
            # the window of test bases need to be updated (the first one should be deleted and a new one should be appened at the end)
            xTestupdate = xTest[1]
            popfirst!(xTest)

            mnow = j + ql - N
            if mnow < m
                # a new basis should be computed
                copyto!(xTestupdate, xTest[end-1])
                # Note that n = mnow - 1 and alpha = beta = N - 1/2 are the indexes for recurrence now
                ttr!(xTestupdate, xTest[end], GCx, mnow - 1, N - 1/2, rrtemp)

                push!(xTest, xTestupdate)
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

function GSBSPG_Chebyshev_NI_solve(::Type{T}, lincoeffs::Vector{Vector{T}}, K::BandedMatrix{T}, bc::AbstractVector{T}, fc::AbstractVector{T}) where {T}
    # construct the GSBSPG discretization (numerical integration) of operator lincoeffs[1]*u + lincoeffs[2]*u' + ... + lincoeffs[N+1]*u^{N} by Chebyshev T polynomials and solve by banded solver

    n = size(K, 2)
    if length(fc) >= n
        f = fc[1:n]
    else
        f = copy(fc)
        append!(f, Zeros(T, n - length(fc)))
    end

    # matrix
    L, u = GSBSPG_Chebyshev_NI(T, lincoeffs, K, bc, f)

    # standard solver (note that L is constructed with extra upper diagonals for factorization)
    ldiv!(lu!(L), u)
    
    # convert back to Chebyshev series x = K * x + bc
    Wlmul!(K, u)
    if length(u) > length(bc)
        axpy!(true, bc, view(u, 1:length(bc)))
    else
        axpy!(true, view(bc, 1:length(u)), u)
    end

    u
end

function ttr!(Pnm1::Vector{T}, Pn::Vector{T}, x::Vector{T}, n::Integer, alpha::Number, temp::Vector{T}) where {T}
    # perform three-term recurrence of orthogonal polynomials, i.e., P_{n+1}^{alpha} = a*x*P_{n}^{alpha} - b*x*P_{n-1}^{alpha}, where P_{n}^{alpha} is the diagonal Jacobi polynomial with index alpha and a and b are coefficients determined from n and alpha
    # the old polynomial Pnm1 is rewritten by the new polynomial, whose values on the same grid are stored (Gauss-Chebyshev points for example)

    # CAUTION: no assertation of the same size among Pnm1, Pn and x (as well as temp)
    # compute a and b (note that b is in opposite sign)
    a = ((n+alpha+1)/(n+1)) * ((2*(n+alpha)+1)/(n+2*alpha+1))
    b = -((n + alpha)/(n + 2*alpha + 1))*((n+alpha+1)/(n+1))

    # three-term recurrence in a neat expression
    broadcast!(*, temp, Pn, x)
    axpby!(a, temp, b, Pnm1)

    Pnm1
end

function Chebyshev_rhs_NI(::Type{T}, f::Function, m::Integer, n::Integer, k::Integer) where {T}
    # compute the rhs in GSBSPG, i.e., f = B_{k}^{k} fhat where fhat = (If, Qtilde)_N and If is the interpolation of f on quadrature points and Qtilde = h\Q is scaled test function

    # quadrature points and weights
    nq = n + k
    GCx, GCw = gausschebyshevt(nq)

    # values of f on quadrature points
    fx = f.(GCx)

    # test functions (unscaled Chebyshev polynomials)
    Tnm1 = ones(T, nq)
    Tn = copy(GCx)
    # and scaling factors (divided later)
    hT = Chebyshev_nfactor(T, n + k, 0)

    # iterate to compute the inner product
    fc = Vector{T}(undef, n + k)
    temp = Vector{T}(undef, nq)
    for i in eachindex(fc)
        # inner product
        broadcast!(*, temp, fx, Tnm1)
        fc[i] = dot(GCw, temp) / hT[i]

        # three-term recurrence for upcoming polynomials
        if i < n + k
            broadcast!(*, temp, GCx, Tn)
            Tnm1, Tn = Tn, Tnm1 
            axpby!(2, temp, -1, Tn)
        end
    end

    # multiplying B_{k}^{k}
    fc = Chebyshev_B(T, m, n + k, k, k) * fc

    fc
end

function Legendre_rhs_NI(::Type{T}, f::Function, m::Integer, n::Integer, k::Integer) where {T}
    # compute the rhs in GSBSPG, i.e., f = B_{k}^{k} fhat where fhat = (If, Qtilde)_N and If is the interpolation of f on quadrature points and Qtilde = h\Q is scaled test function

    # quadrature points and weights
    nq = n + k
    GLx, GLw = gausslegendre(nq)

    # values of f on quadrature points
    fx = f.(GLx)

    # test functions (unscaled Legendre polynomials)
    Tnm1 = ones(T, nq)
    Tn = copy(GLx)
    # and scaling factors (divided later)
    hT = Legendre_nfactor(T, n + k, 0)

    # iterate to compute the inner product
    fc = Vector{T}(undef, n + k)
    temp = Vector{T}(undef, nq)
    for i in eachindex(fc)
        # inner product
        broadcast!(*, temp, fx, Tnm1)
        fc[i] = dot(GLw, temp) / hT[i]

        # three-term recurrence for upcoming polynomials
        if i < n + k
            broadcast!(*, temp, GLx, Tn)
            Tnm1, Tn = Tn, Tnm1
            # Note that n = i is the index for recurrence now
            axpby!((2*i+1)/(i+1), temp, -i/(i+1), Tn)
        end
    end

    # multiplying B_{k}^{k}
    fc = Legendre_B(T, m, n + k, k, k) * fc

    fc
end