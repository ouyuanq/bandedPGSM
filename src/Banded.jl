# some useful functions for BandedMatrices
import LinearAlgebra: lmul!, ldiv!

function Sscale!(M::BandedMatrix, mu::Integer, bw::Integer)
    # overwirte M with DS * M * DS^{-1} where DS = diagonal([mu/mu mu/(mu+1) mu/(mu+2) ...], 0) is the scaling factors in conversion operator
    # only elements in (bw, bw) bandwidths are executed

    scalefactor = mu:mu+max(size(M)...)-1
    Mdata, Mu = M.data, M.u
    @inbounds for j in axes(M, 2)
        for i in max(j - bw, 1):min(j + bw, size(M, 1))
            Mij = inbands_getindex(Mdata, Mu, i, j) * (scalefactor[j] / scalefactor[i])
            inbands_setindex!(Mdata, Mu, Mij, i, j)
        end
    end

    M
end

function Slmul!(M::BandedMatrix, bw::Integer)
    # overwirte M with (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2)) * M
    # Note that the bandwidths of M is not enlarged
    # only elements in (bw, bw) bandwidths are executed

    Mdata, Mu = M.data, M.u
    @inbounds for j in axes(M, 2)
        for i = max(j - bw, 1):min(j + bw, size(M, 1))-2
            Mij = inbands_getindex(Mdata, Mu, i, j) - inbands_getindex(Mdata, Mu, i + 2, j)
            inbands_setindex!(Mdata, Mu, Mij, i, j)
        end
    end

    M
end

function Srdiv!(M::BandedMatrix, bw::Integer)
    # overwirte M with M / (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2))
    # Note that the bandwidths of M is not enlarged

    n = size(M, 2)
    Mdata, Mu = M.data, M.u
    @inbounds for j = 3:n
        j2 = j - 2
        colj = max(j - bw, 1):min(j + bw, size(M, 1))
        colj2 = max(j2 - bw, 1):min(j2 + bw, size(M, 1))
        for i in colj2 ∩ colj
            Mij = inbands_getindex(Mdata, Mu, i, j) + inbands_getindex(Mdata, Mu, i, j2)
            inbands_setindex!(Mdata, Mu, Mij, i, j)
        end
    end

    M
end

function Slmul!(A::BandedMatrix, lambda::Integer, mu::Integer, Alnow::Integer, Aunow::Integer)
    # left application of S_{mu-1}S_{mu-2}...S_{lambda} to M
    # Note that S_{mu} = (diagonal([1 1 1 ...], 0) + diagonal([-1 -1 -1 ...], 2)) * diagonal([mu/mu mu/(mu+1) mu/(mu+2) ...], 0), i.e., a bidiagonal matrix times a diagonal matrix
    # the bandwidths are enlarged and only Al lower bandwidth and Au upper bandwidth of A are considered

    @assert Alnow <= A.l && Aunow + 2*(mu-lambda) <= A.u "Not enough bandwidths of in-place operations for Slmul"

    Adata, Au = A.data, A.u
    szA2 = size(A, 2)
    for alpha = lambda:mu-1
        # scale
        if alpha > 0
            scalefactor = alpha:alpha+szA2-1
            @inbounds for j in axes(A, 2)
                for i in max(j-Aunow, 1):min(j+Alnow, szA2)
                    # why size(A, 2) instead of size(A, 1)? Since we have extra storage beyond szA1 rows and these rows can be used to generate elements within szA1 rows. szA2 = szA1 + N is assumed where N is the maximum order of differentiation
                    Aij = inbands_getindex(Adata, Au, i, j) * (alpha / scalefactor[i])
                    inbands_setindex!(Adata, Au, Aij, i, j)
                end  
            end
        else
            ldiv!(2, view(Adata, Au+1-Aunow:Au+1+Alnow, :))
            @inbounds for j in 1:min(1+Aunow, size(Adata, 2))
                Aij = inbands_getindex(Adata, Au, 1, j) * 2
                inbands_setindex!(Adata, Au, Aij, 1, j)
            end
        end

        # bidiagonal
        @inbounds for j in axes(A, 2)
            for i in max(j-Aunow-2, 1):min(j+Alnow, szA2)-2
                # size(A, 2) for the same reason above
                Aij = inbands_getindex(Adata, Au, i, j) - inbands_getindex(Adata, Au, i+2, j)
                inbands_setindex!(Adata, Au, Aij, i, j)
            end
        end

        Aunow += 2
    end

    A
end

function Wlmul!(W::BandedMatrix{T}, x::Vector{T}) where T
    # overwirte x with W * x, where W is banded lowertriangular matrix
    # x would be prolonged
    @assert W.u == 0 "Only lowertriangular matrix are allowed"
    Wl = W.l
    if size(W, 2) == length(x)
        # prolong x
        append!(x, Zeros(T, Wl))
    elseif size(W, 1) == length(x)
        # only first size(W, 2) elements of x are taken into consideration
        x[size(W, 2)+1:end] .= 0
    else
        error("Length of x does not equal to either side of W")
    end

    Wdata = W.data
    @inbounds for j in reverse(axes(W, 2))
        xj = x[j]
        x[j] = 0
        # axpy!(xj, view(Wdata, :, j), view(x, j:j+Wl))
        for i in reverse(colrange(W, j))
            Wij = inbands_getindex(Wdata, 0, i, j)
            x[i] = muladd(Wij, xj, x[i])
        end
    end

    x
end

function Wtlmul!(W::BandedMatrix{T}, x::Vector{T}) where T
    # overwirte x with W^T * x, where W is banded lowertriangular matrix
    # x would be truncated
    @assert W.u == 0 "Only lowertriangular matrix are allowed"
    @assert size(W, 1) == length(x) "Wt and x can not be multiplied"

    Wl = W.l
    Wdata = W.data
    @inbounds for j in axes(W, 2)
        x[j] *= inbands_getindex(Wdata, 0, j, j)
        for i in j+1:j+Wl
            Wij = inbands_getindex(Wdata, 0, i, j)
            x[j] = muladd(Wij, x[i], x[j])
        end
    end

    # truncate x to proper dimension
    deleteat!(x, size(W, 2)+1:length(x))

    x
end

function Wrmul!(A::BandedMatrix{T}, W::BandedMatrix{T}, Au::Integer) where T
    # overwirte A with A * W, where W is a lowertriangular matrix with lower bandwidth Wl and only Au superdiagonals of A are nonzero

    @assert W.u == 0 "Only lowertriangular matrix are allowed"
    @assert size(A, 2) == size(W, 1) "Incompatible dimensions of A and W"
    N = W.l
    Adata, Wdata, Wu = A.data, W.data, W.u
    szAdata1 = size(Adata, 1)
    Au1 = A.u+1-Au
    @inbounds for j in axes(Wdata, 2)
        lmul!(inbands_getindex(Wdata, Wu, j, j), view(Adata, :, j))
        for k = 1:N
            axpy!(inbands_getindex(Wdata, Wu, j+k, j), view(Adata, Au1:szAdata1-k, j+k), view(Adata, Au1+k:szAdata1, j))
        end
    end

    # discard the last N columns
    _BandedMatrix(view(Adata, :, 1:size(W, 2)), Base.OneTo(size(A, 1)), A.l, A.u)
end

function Wrmul!(A::Matrix{T}, W::BandedMatrix{T}) where T
    # overwirte A with A * W, where W is a lowertriangular matrix with lower bandwidth Wl and A is a dense matrix

    @assert W.u == 0 "Only lowertriangular matrix are allowed"
    @assert size(A, 2) == size(W, 1) "Incompatible dimensions of A and W"
    N = W.l
    Wdata, Wu = W.data, W.u
    @inbounds for j in axes(Wdata, 2)
        lmul!(inbands_getindex(Wdata, Wu, j, j), view(A, :, j))
        for k = 1:N
            axpy!(inbands_getindex(Wdata, Wu, j+k, j), view(A, :, j+k), view(A, :, j))
        end
    end

    # discard the last N columns
    view(A, :, 1:size(W, 2))
end

function Wtlmul!(W::BandedMatrix{T}, A::BandedMatrix{T}, Aunow::Integer) where T
    # overwirte A with W^T * A, where W is banded lowertriangular matrix and only Aunow superdiagonals of A are nonzero
    @assert W.u == 0 "Only lowertriangular matrix are allowed"
    @assert size(W, 1) == size(A, 1) "Wt and A can not be multiplied"
    @assert A.u >= Aunow + W.l "Not enough storage for in-place multiplication"

    Wl, szW2 = W.l, size(W, 2)
    Aunow1 = Aunow + Wl  # there would be Wl more superdiagonals
    Adata, Al, Au = A.data, A.l, A.u
    Wdata = W.data
    @inbounds for j in axes(A, 2)
        for i = max(1, j - Aunow1):min(szW2, j + Al)
            Aij = zero(T)
            ij = j - i
            for k = max(0, ij - Aunow1):min(Wl, Al + ij)
                Aij = muladd(Wdata[k+1, i], inbands_getindex(Adata, Au, k + i, j), Aij)
            end
            inbands_setindex!(Adata, Au, Aij, i, j)
        end
    end

    # discard the last N rows
    _BandedMatrix(Adata, Base.OneTo(szW2), Al, Au)
end

function Dlmul!(D::Diagonal{T}, A::BandedMatrix{T}) where T
    # in-place left multiplication of a diagonal matrix D to a banded matrix A
    # the results are stored in A

    @assert size(D, 2) == size(A, 1) "Incompatible dimension of D and A"
    d = D.diag
    Adata, Al, Au = A.data, A.l, A.u
    szAdata1, szA1 = size(Adata, 1), size(A, 1)
    @inbounds for j in axes(Adata, 2)
        Aj = view(Adata, max(1, Au+2-j):min(szAdata1, Au+1+szA1-j), j)
        dj = view(d, max(1, j-Au):min(szA1, j+Al))
        broadcast!(*, Aj, Aj, dj)
    end

    A
end

# left multiplication of a diagonal matrix for a banded matrix
function lmul!(D::Diagonal{T}, B::BandedMatrix{T}) where {T}
    # compute B = D \ B

    d = D.diag
    m = length(d)
    @assert m == size(B, 1) "The sizes of diagonal matrix and banded matrix are not compatible"

    Bdata = B.data
    Bu = B.u
    szBdata1 = size(Bdata, 1)
    szB1 = size(B, 1)
    for j in axes(Bdata, 2)
        colj = colrange(B, j)
        dj = view(d, colj)
        Bdataj = view(Bdata, max(1, Bu+2-j):min(szBdata1, Bu+1+szB1-j), j)
        broadcast!(*, Bdataj, Bdataj, dj)
    end

    B
end

# left division of a diagonal matrix for a banded matrix
function ldiv!(D::Diagonal{T}, B::BandedMatrix{T}) where {T}
    # compute B = D \ B

    d = D.diag
    m = length(d)
    @assert m == size(B, 1) "The sizes of diagonal matrix and banded matrix are not compatible"

    Bdata = B.data
    Bu = B.u
    szBdata1 = size(Bdata, 1)
    szB1 = size(B, 1)
    for j in axes(Bdata, 2)
        colj = colrange(B, j)
        dj = view(d, colj)
        Bdataj = view(Bdata, max(1, Bu+2-j):min(szBdata1, Bu+1+szB1-j), j)
        broadcast!(/, Bdataj, Bdataj, dj)
    end

    B
end

function mybandedmul!(A::BandedMatrix{T}, B::BandedMatrix{T}, C::BandedMatrix{T}, alpha::Number, beta::Number) where {T}
    # compute A = alpha * B * C + beta * A

    Au, Al, Bu, Bl, Cu, Cl = A.u, A.l, B.u, B.l, C.u, C.l
    @assert A.u >= B.u + C.u && A.l >= B.l + C.l "The destination A can not contain the result of B*C."
    @assert size(A, 1) == size(B, 1) && size(A, 2) == size(C, 2) && size(B, 2) == size(C, 1) "Incompatible dimensions of A + B*C"

    Adata, Bdata, Cdata = A.data, B.data, C.data
    szB1, szBdata1 = size(B, 1), size(Bdata, 1)
    szC1 = size(C, 1)
    # multiplication by each column
    for j in axes(Cdata, 2)
        for k in max(1, j-Cu):min(szC1, j+Cl)
            Ckj = Cdata[Cu+1+(k-j), j] * alpha
            Bk = view(Bdata, max(1, Bu+2-k):min(szBdata1, Bu+1+szB1-k), k)
            shift = Au+1-j
            Ajind = max(1, k-Bu)+shift:min(szB1, k+Bl)+shift
            Aj = view(Adata, Ajind, j)
            axpby!(Ckj, Bk, beta, Aj)
        end
    end

    A
end