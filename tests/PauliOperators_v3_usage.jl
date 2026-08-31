#
# General usage of PauliOperators for evolving routines
#
using PauliOperators
using LinearAlgebra

# Create Pauli operators (3 qubits)
X = Pauli("XII")
Y = Pauli("IYI")
Z = Pauli("IIZ")

# Build a Hamiltonian
H = PauliSum(3, ComplexF64)
H[PauliBasis("XII")] = -1.0
H[PauliBasis("IYI")] = -0.5
H[PauliBasis("IIZ")] = -0.3

# Create a quantum state
ψ = Ket([1, 0, 1])  # |101⟩

# Compute expectation value
E = expectation_value(H, ψ)

# Tensor products
ZZ = Z ⊗ Z  # Create a 6-qubit operator

# Direct sums
H_total = H ⊕ H  # Combine two 3-qubit Hamiltonians

#Standard linear algebra.norm is extended for PauliSum and KetSum
norm(H)        # L2 norm (default)
norm(H, 1)     # L1 norm (sum of |c_k|)
norm(H, Inf)   # L∞ norm (max |c_k|)
#isapprox(H1, H2; atol=1e-10)  # Approximate equality

#
# Truncation Strategies
#
# Available strategies
strat = CoeffTruncation(1e-6)             # Drop terms with |c| < threshold
#strat = WeightTruncation(3)               # Keep terms with Pauli weight ≤ 3
#strat = XWeightTruncation(2)              # Keep terms with ≤ 2 off-diagonal (X/Y) factors
#strat = MajoranaWeightTruncation(4)       # Keep terms with Majorana weight ≤ 4
#strat = WeightDampedTruncation(α, ε)      # Drop terms with |c|·e^(-α·weight) ≤ ε
#strat = XWeightDampedTruncation(α, ε)     # Drop terms with |c|·e^(-α·x_weight) ≤ ε
#strat = AdaptiveTruncation(1000, 1e-8)    # Keep at most 1000 terms, min threshold 1e-8
#strat = CompositeTruncation(s1, s2, ...)  # Apply multiple strategies in sequence

#
# First-order Trotter
dt = 0.8
gens, angs = trotterize(H_total, dt; n_trotter=1)

# Second-order (symmetric) Trotter
#gens, angs = trotterize(H, dt; order=2)

# QDrift stochastic protocol
#gens, angs = qdrift(H, dt; n_samples=100)

println(" EXAMPLE WITH A GENERAL HAMILTONIAN ...")
println("# # # Generators and angles # # #")
println(gens)
println(angs)

#
# - - - EVOLUTION - - - 
#
# Create a quantum state (6qubits)
ψ = Ket([0, 0, 0, 0, 0, 0])  # |000000⟩

#Single-generator evolution
XX = X ⊗ Pauli("III")
O_sparse = SparsePauliVector(XX)
println("# # # O_sparse # # #")
display(O_sparse)

O_evolved = evolve(O_sparse, gens, angs)
         # gather back into a Dict-based PauliSum anytime
println("# # # O_evolved single-generator evolution # # #")
display(O_evolved)

#Sequence evolution with truncation
O_out = evolve(O_sparse, gens, angs;
               truncation=CoeffTruncation(1e-8),
               correction=EnergyCorrection(ψ))

println("# # # O_out sequence evolution # # #")               
display(O_out)

E = expectation_value(O_out, ψ)   # observables work directly on the flat form
println("ExpVal:", E)
#O = PauliSum(O_out)               # gather back into a Dict-based PauliSum anytime

using QuantumChemQC
using NPZ
println(" EXAMPLE WITH H2 HAMILTONIAN...")
# ============================================================
# 1. SETTINGS
# ============================================================
function build_H()
    # Get Molecular Hamiltonian
    data_path =  joinpath(@__DIR__, "h2-RHF_test_ham.npz")
    data = npzread(data_path)
    H0 = data["hc"][1]
    H1 = data["h1e"]
    H2 = data["h2e"]

    println("H0 shape: ", size(H0))
    println(typeof(H0))
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))

    Norbs = size(H1,1)  # number of spatial orbitals
    
    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=true)

    return H
end

# Initial operator, ket and Hamiltonian
O = Pauli(4, X=[1,3])
O = PauliSum(O)
O += Pauli(4, X=[2,3])
H2_ham = build_H()
ket, _ = QuantumChemQC.string_to_ket("0000")

println("H2 Hamiltonian")
display(H2_ham)
# First-order Trotter
# Time evolution parameters
n_intervals = 200
t = 50.0
dt = t/n_intervals

#Pruning parameters
evol_thresh = 0.01
max_weight = 4 

println()
println("Hamiltonian terms: ", length(H))
println("Initial operator terms: ", length(O))
println("Time step dt: ", dt)
println("Number of intervals: ", n_intervals)
println("Pruning threshold: ", evol_thresh)

gens, angs = trotterize(H2_ham, dt; n_trotter=n_intervals)
println("GENERATORS")
println(gens)
println("ANGLES")
println(angs)

# ------------------------------------------------------------
# Create the SparsePauliVector versions
# ------------------------------------------------------------
O_sparse = SparsePauliVector(O)

H_sparse = SparsePauliVector(H)

println("Rotations per interval: ", length(angles))
println("Total rotations: ", length(angles) * n_intervals)

#Sequence evolution with truncation
O_out = evolve(O_sparse, gens, angs;
               truncation=CoeffTruncation(1e-6),
               correction=EnergyCorrection(ket))

println("# # # O_out sequence evolution # # #")               
display(O_out)

E = expectation_value(O_out, ket)   # observables work directly on the flat form
println("ExpVal:", E)