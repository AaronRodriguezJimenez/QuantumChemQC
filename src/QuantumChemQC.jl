module QuantumChemQC

# Module for creating quantum chemistry Hamiltonians
# and using them for quantum-computing purposes.

using Tullio, LinearAlgebra, SpecialFunctions, GaussianBasis,
      FileIO, StaticArrays, IterTools, Graphs, Printf

using PauliOperators

# Explicitly import functions that QuantumChemQC extends.
import PauliOperators:
    coeff_clip!,
    weight_clip!,
    majorana_weight_clip!

using NPZ

include("io.jl")
include("integrals.jl")
include("param.jl")
include("scf.jl")
include("diis.jl")
include("utils.jl")
include("qubit_utils.jl")
include("hamiltonians.jl")
include("qsp.jl")

# Re-export the same PauliOperators function bindings.
export coeff_clip!,
       weight_clip!,
       majorana_weight_clip!

end # module QuantumChemQC