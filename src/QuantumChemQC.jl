module QuantumChemQC

#Module for creating quantum chemistry hamiltonians
#and use them for quantum computing purposes.

using Tullio, LinearAlgebra, SpecialFunctions, GaussianBasis,
      FileIO, StaticArrays, IterTools, Graphs, Printf
using PauliOperators   # v2
using DBF

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
include("lattices/Lattices.jl")  # defines submodule

using .QCUtils
using .Lattices                  # make it accessible inside QuantumChemQC
using .Molecules

export FermionOp, FermionOperator
export Molecules
export Lattices, LatticeBond, Lattice
export lattice2graph, dec2bin, bin2dec, bin2bonds!
export square_lattice


end # module QuantumChemQC