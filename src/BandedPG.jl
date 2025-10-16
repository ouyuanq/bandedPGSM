using Plots, Printf, DelimitedFiles, BenchmarkTools

using LinearAlgebra, BandedMatrices, FFTW, SemiseparableMatrices, FillArrays, ApproxFun, SpecialFunctions, FastGaussQuadrature
using BandedMatrices: BandedMatrices, inbands_getindex, inbands_setindex!, _BandedMatrix
import Polynomials.coeffs as pcoeffs
import Polynomials.fit as pfit

# basic functions
include("Operators.jl")
include("Chebyshev.jl")
include("Multiplications.jl")
include("BasisRecombination.jl")
include("MatrixPG.jl")
include("MatrixUS.jl")
include("Banded.jl")
include("Vandermonde.jl")

# GSBSPG
include("GSBSPG.jl")
include("GSBSPG_NI.jl")