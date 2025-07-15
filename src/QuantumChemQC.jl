module QuantumChemQC

#Module for creating quantum chemistry hamiltonians
#and use them for quantum computing purposes.
using Tullio
using LinearAlgebra
using SpecialFunctions
using GaussianBasis
using FileIO
using StaticArrays
using IterTools

using PauliOperators  # v2

include("io.jl")
include("integrals.jl")
include("param.jl")
include("scf.jl")
include("diis.jl")
include("utils.jl")
include("fermion_utils.jl")
include("qubit_utils.jl")
include("type_FermionOp.jl")

export FermionOp
export FermionOperator
end # module QuantumChemQC
