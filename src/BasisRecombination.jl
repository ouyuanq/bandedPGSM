## matrices corresponding to basis recombination
## see our codes in Mathematica for analytically computing the coefficients

## test and trial basis (Chebyshev)
function secondDirichlet_ChebyshevW(::Type{T}, n::Integer) where T
    # the recombined basis for seconder order equation with homogeneous Dirichlet boundary conditions
    Wtrial = BandedMatrix{T}((0=>Fill(-1, n-2), -1=>Fill(0, n-2), -2=>Fill(1, n-2)), (n, n-2), (2, 0))
    Wtest = BandedMatrix{T}((0=>Fill(0, n-2), -1=>Fill(0, n-2), -2=>Fill(1, n-2)), (n, n-2), (2, 0))
    @inbounds for i = 1:n-2
        Wtest[i, i] = -((i+3)*(i+4)) / (i*(i+1))
    end
    broadcast!(/, Wtest.data, Wtest.data, (4:2:2*n-2)')
    d = Vector{T}(undef, n)
    con = pi / 2^3
    @inbounds for i in eachindex(d)
        d[i] = i*(i+2)*con
    end
    Omega = Diagonal(d)

    Wtrial, Wtest, Omega
end

function leftDirichlet_ChebyshevW(::Type{T}, n::Integer) where T
    # the recombined basis for first order equation with homogeneous left Dirichlet boundary conditions
    Wtrial = BandedMatrix{T}((0=>Fill(1, n-1), -1 => Fill(1, n-1)), (n, n-1), (1, 0))
    Wtest = BandedMatrix{T}((0=>Fill(0, n-1), -1 => Fill(1, n-1)), (n, n-1), (1, 0))
    @inbounds for i = 1:n-1
        Wtest[i, i] = (i+1) / i
    end
    d = Vector{T}(undef, n)
    d .= pi / 2
    Omega = Diagonal(d)

    Wtrial, Wtest, Omega
end

function tenthorder(::Type{T}, n::Integer) where T
    # transformation operator related to zeroth to fourth derivatives of left and right end point conditions

    R = BandedMatrix{T}(Zeros(T, n, n-10), (10, 0))
    Rdata = R.data
    # 10-th lower diagonal
    view(Rdata, 11, :) .= 1
    # 8-th lower diagonal
    diagnow = Vector{T}(8:n-3)
    broadcast!(/, diagnow, diagnow, 4:n-7)
    axpy!(-5, diagnow, view(Rdata, 9, :))
    # 6-th lower diagonal
    diagnow .= 6:n-5
    broadcast!(*, diagnow, diagnow, 9:n-2)
    broadcast!(/, diagnow, diagnow, 3:n-8)
    broadcast!(/, diagnow, diagnow, 4:n-7)
    axpy!(10, diagnow, view(Rdata, 7, :))
    # 4-th lower diagonal
    diagnow .= 8:n-3
    broadcast!(*, diagnow, diagnow, 9:n-2)
    broadcast!(/, diagnow, diagnow, 2:n-9)
    broadcast!(/, diagnow, diagnow, 3:n-8)
    axpy!(-10, diagnow, view(Rdata, 5, :))
    # 2-th lower diagonal
    diagnow .= 7:n-4
    broadcast!(*, diagnow, diagnow, 8:n-3)
    broadcast!(*, diagnow, diagnow, 9:n-2)
    broadcast!(/, diagnow, diagnow, 1:n-10)
    broadcast!(/, diagnow, diagnow, 3:n-8)
    broadcast!(/, diagnow, diagnow, 4:n-7)
    axpy!(5, diagnow, view(Rdata, 3, :))
    # 0-th lower diagonal (main diagonal)
    diagnow .= 6:n-5
    broadcast!(*, diagnow, diagnow, 7:n-4)
    broadcast!(*, diagnow, diagnow, 8:n-3)
    broadcast!(*, diagnow, diagnow, 9:n-2)
    broadcast!(/, diagnow, diagnow, 1:n-10)
    broadcast!(/, diagnow, diagnow, 2:n-9)
    broadcast!(/, diagnow, diagnow, 3:n-8)
    broadcast!(/, diagnow, diagnow, 4:n-7)
    axpy!(-1, diagnow, view(Rdata, 1, :))

    R
end