## matrices corresponding to basis recombination
## see our codes in Mathematica for analytically computing the coefficients

## test and trial basis (Chebyshev)
function secondDirichlet_Chebyshev(::Type{T}, n::Integer) where T
    # the recombined basis for second order equation with homogeneous Dirichlet boundary conditions
    R = BandedMatrix{T}((0=>Fill(-1, n-2), -1=>Fill(0, n-2), -2=>Fill(1, n-2)), (n, n-2), (2, 0))

    Q = BandedMatrix{T}((0=>Fill(0, n-2), -1=>Fill(0, n-2), -2=>Fill(1, n-2)), (n, n-2), (2, 0))
    for i in axes(Q, 2)
        Q[i, i] = -((i+3)*(i+4)) / (i*(i+1))
    end
    broadcast!(/, Q.data, Q.data, (4:2:2*n-2)')

    d = Vector{T}(undef, n)
    con = pi / 2^3
    for i in eachindex(d)
        d[i] = i*(i+2)*con
    end
    Omega = Diagonal(d)

    R, Q, Omega
end

function leftDirichlet_Chebyshev(::Type{T}, n::Integer) where T
    # the recombined basis for first order equation with homogeneous left Dirichlet boundary conditions
    R = BandedMatrix{T}((0=>Fill(1, n-1), -1 => Fill(1, n-1)), (n, n-1), (1, 0))

    Q = BandedMatrix{T}((0=>Fill(0, n-1), -1 => Fill(1, n-1)), (n, n-1), (1, 0))
    for i in axes(Q, 2)
        Q[i, i] = (i+1) / i
    end

    d = Vector{T}(undef, n)
    d .= pi / 2
    Omega = Diagonal(d)

    R, Q, Omega
end

function thirdrightNeumann_Chebyshev(::Type{T}, n::Integer) where T
    # the recombined basis for third order equation with homogeneous Dirichlet boundary conditions and right Neumann boundary condition
    R = BandedMatrix{T}((-1=>Fill(-1, n-3), -3=>Fill(1, n-3)), (n, n-3), (3, 0))
    for i in axes(R, 2)
        R[i, i] = (i+1) / i
        R[i+2, i] = -R[i, i]
    end

    Q = BandedMatrix{T}((-3=>Fill(1, n-3),), (n, n-3), (3, 0))
    for i in axes(Q, 2)
        Q[i, i] = -((i+4)*(i+5)*(i+6)*(i+7)) / (i*(i+1)*(i+2)*(i+3))
        Q[i+1, i] = -((i+6)*(i+7)) / ((i+1)*(i+2))
        Q[i+2, i] = ((i+4)*(i+7)) / ((i+2)*(i+3))
    end
    broadcast!(/, Q.data, Q.data, (8 .* (3:n-1))')

    d = Vector{T}(undef, n)
    con = pi / 2^7
    for i in eachindex(d)
        d[i] = i*(i+1)*(i+3)*(i+4)*con
    end
    Omega = Diagonal(d)

    R, Q, Omega
end

function thirdrightNeumann_Chebyshev_nondual(::Type{T}, n::Integer) where T
    # the recombined basis for third order equation with homogeneous Dirichlet boundary conditions and left Neumann boundary condition
    R = BandedMatrix{T}((-1=>Fill(-1, n-3), -3=>Fill(1, n-3)), (n, n-3), (3, 0))
    for i in axes(R, 2)
        R[i, i] = (i+1) / i
        R[i+2, i] = -R[i, i]
    end

    Q = BandedMatrix{T}((-3=>Fill(1, n-3),), (n, n-3), (3, 0))
    for i in axes(Q, 2)
        Q[i, i] = ((i+4)*(i+5)*(i+6)*(i+7)) / (i*(i+1)*(i+2)*(i+3))
        Q[i+1, i] = -((i+6)*(i+7)) / ((i+1)*(i+2))
        Q[i+2, i] = -((i+4)*(i+7)) / ((i+2)*(i+3))
    end
    broadcast!(/, Q.data, Q.data, (8 .* (3:n-1))')

    d = Vector{T}(undef, n)
    con = pi / 2^7
    for i in eachindex(d)
        d[i] = i*(i+1)*(i+3)*(i+4)*con
    end
    Omega = Diagonal(d)

    R, Q, Omega
end

function fourthDirichletNeumann_Chebyshev(::Type{T}, n::Integer) where T
    # the recombined basis for fourth order equation with homogeneous Dirichlet and Neumann boundary conditions on both ends
    R = BandedMatrix{T}((-4=>Fill(1, n-4), ), (n, n-4), (4, 0))
    for i in axes(R, 2)
        R[i, i] = (i+2) / i
        R[i+2, i] = -2*((i+1)/i)
    end

    Q = BandedMatrix{T}((-4=>Fill(1, n-4), ), (n, n-4), (4, 0))
    for i in axes(Q, 2)
        Q[i, i] = ((i+6)/(i+4))*((i+7)/(i+3))*((i+8)/(i+2))*((i+9)/(i+1))*((i+10)/i)
        Q[i+2, i] = -2*((i+5)/(i+4))*((i+9)/(i+3))*((i+10)/(i+2))
    end

    d = Vector{T}(undef, n)
    con = pi / (2^7 * 36)
    for i in eachindex(d)
        d[i] = i*(i+1)*(i+2)*(i+4)*(i+5)*(i+6)*con
    end
    Omega = Diagonal(d)

    R, Q, Omega
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