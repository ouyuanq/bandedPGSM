using Plots, Printf, DelimitedFiles, BenchmarkTools

using LinearAlgebra, BandedMatrices, FFTW, SemiseparableMatrices, FillArrays, ApproxFun, SpecialFunctions, FastGaussQuadrature
using BandedMatrices: BandedMatrices, inbands_getindex, inbands_setindex!, _BandedMatrix

# basic functions
include("Operators.jl")
include("Chebyshev.jl")
include("Mutiplications.jl")
include("BasisRecombination.jl")
include("MatrixPG.jl")
include("MatrixUS.jl")
include("Banded.jl")

# GSBSPG
include("GSBSPG.jl")
include("GSBSPG_NI.jl")