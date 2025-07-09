module QuantumChemQC

#Module for creating quantum chemistry hamiltonians
#and use them for quantum computing purposes.
using Tullio
using LinearAlgebra
using SpecialFunctions
using GaussianBasis
using FileIO
using StaticArrays

using PauliOperators  # v2

include("io.jl")
include("integrals.jl")
include("param.jl")
include("scf.jl")
include("diis.jl")
include("utils.jl")

end # module QuantumChemQC
