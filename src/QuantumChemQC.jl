module QuantumChemQC

#Module for creating quantum chemistry hamiltonians
#and use them for quantum computing purposes.

using Tullio, LinearAlgebra, SpecialFunctions, GaussianBasis,
      FileIO, StaticArrays, IterTools, Graphs, Printf

using PauliOperators   # v3

# Explicitly import functions that QuantumChemQC also extends.
import PauliOperators:
    coeff_clip!,
    weight_clip!,
    majorana_weight_clip!

using DBF
using Random
using StatsBase
using Statistics
using NPZ

include("io.jl")
include("integrals.jl")
include("param.jl")
include("scf.jl")
include("diis.jl")
include("utils.jl")
include("type_FermionOp.jl")
include("fermion_utils.jl")
include("qubit_utils.jl")
include("hamiltonians.jl")
include("qc_utils.jl") #submodule QCUtils
include("molecules.jl")
include("evolve.jl")
include("qsp.jl")

using .QCUtils     # make it accessible inside QuantumChemQC
using .Molecules

export FermionOp, FermionOperator
export Molecules
export coeff_clip!, weight_clip!, majorana_weight_clip!


end # module QuantumChemQC