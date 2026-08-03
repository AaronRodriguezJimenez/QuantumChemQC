#
# SparsePauliVector is an alternative storage engine for Sums of Paulis,
# designed for evolution-heavy workloads. It stores the terms as flat,
# stored, preallocated parallel arrays instead of a Dict:
# rotations become linear array sweeps, deduplications becomes a sort-merge,
# and the steady-state hot path allocates zero bytes.
# 
using PauliOperators

# Build with the convenient Dict-based API, then convert to flat storage
H = PauliSum(4, ComplexF64)
H[PauliBasis("ZZII")] = 1.0
H[PauliBasis("IZZI")] = 1.0
H[PauliBasis("XXII")] = 0.5
v = SparsePauliVector(H; T=Float64)   # real coefficients suffice for Hermitian operators

# Same evolution API as PauliSum, plus windowed merging
ψ = Ket([0, 0, 0, 0])
gens, angs = trotterize(H, 0.1; n_trotter=100)
evolve!(v, gens, angs;
        window=10,                        # dedup + truncate every 10 rotations
        truncation=CoeffTruncation(1e-8),
        correction=EnergyCorrection(ψ))

E = expectation_value(v, ψ)   # observables work directly on the flat form
O = PauliSum(v)               # gather back into a Dict-based PauliSum anytime

println("Hamiltonian:")
display(H)
println("SparsePauliVector:")
display(v)
println("Expectation value: ", E)

#
# Playing with SparisePauliVector:

coeff_clip!(v, 1e-4)   # clip small coefficients to zero
println("After coeff_clip!():")
display(v)