using QuantumChemQC
using PauliOperators
using NPZ
using LinearAlgebra
using SparseArrays


# ============================
# Pauli basis indexing
# ============================

@inline function index(p::PauliBasis{N}) where {N}
    return Int(p.z + (p.x << N) + 1)
end

@inline function index(p::Pauli{N}) where {N}
    return index(PauliBasis(p))
end

@inline function basis_from_index(::Val{N}, j::Integer) where {N}
    k = Int128(j - 1)
    mask = (Int128(1) << N) - 1
    z = k & mask
    x = k >> N
    return PauliBasis{N}(z, x)
end

# ============================
# Convert PauliSum <-> coefficient vector and back
# ============================
function coeff_vector(ps::PauliSum{N}; T=ComplexF64) where {N}
    D = Int(4)^N
    v = zeros(T, D)

    for (p, c) in ps
        v[index(p)] += convert(T, c)
    end

    return v
end

function paulisum_from_coeffs(v::AbstractVector, ::Val{N}; atol=1e-14) where {N}
    D = Int(4)^N
    @assert length(v) == D

    out = PauliSum(N, ComplexF64)

    for p in PauliBasis{N}
        c = v[index(p)]
        if abs(c) > atol
            out[p] = ComplexF64(c)
        end
    end

    return out
end

# ============================
# Build the Pauli-Liouville matrix
# Each column is the result of applyiing the Liouville superoperator
# to a Pauli basis element
# ============================
function singleton_paulisum(p::PauliBasis{N}; T=ComplexF64) where {N}
    out = PauliSum(N, T)
    out[p] = one(T)
    return out
end

function liouvillian_matrix(H_in::PauliSum{N}; atol=0.0) where {N}
    D = Int(4)^N
    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]

    for Pj in PauliBasis{N}
        col = index(Pj)

        Pj_sum = singleton_paulisum(Pj)
        LPj = 1im * commutator(H, Pj_sum)

        for (Pi, c) in LPj
            if abs(c) > atol
                push!(rows, index(Pi))
                push!(cols, col)
                push!(vals, c)
            end
        end
    end

    return sparse(rows, cols, vals, D, D)
end

# ============================
# Exact finie-time propagation
# ============================
function evolve_liouville(O::PauliSum{N}, H::PauliSum{N}, t;
                          atol=1e-12) where {N}

    L = liouvillian_matrix(H)
    v0 = coeff_vector(O)

    vt = exp(t * Matrix(L)) * v0

    Ot = paulisum_from_coeffs(vt, Val(N); atol=atol)

    return Ot, vt, L
end

# ============================
# Get input Operators
# ===========================
BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/tensors" # Canonical RHF
path = joinpath(BASE, "RHF_H2_R_0p7414_tensors.npz")
data = npzread(path) 
H1 = data["h1e"] 
Norbs = size(H1,1) 
Nqubits = 2*Norbs

H = QuantumChemQC.molecular_hamiltonian( Norbs, path, NOI = false, block = false, )

O = PauliSum(Pauli(Nqubits, Z = [4]))

t = 5.0

Ot, vt, L = evolve_liouville(O, H, t)

println("* * * * * * O(t) as PauliSum * * * * *")
display(Ot)

println("* * * * * * coefficient vector * * * * *")
display(vt)

# ================================
# Pauli Weight Distribution
# ================================
function pauli_weight_distribution(v::AbstractVector, ::Val{N};
                                   normalize=true,
                                   drop_identity=false) where {N}
    D = Int(4)^N
    @assert length(v) == D

    dist = zeros(Float64, N + 1)

    for p in PauliBasis{N}
        w = PauliOperators.weight(p)
        dist[w + 1] += abs2(v[index(p)])
    end

    if drop_identity
        dist[1] = 0.0
    end

    if normalize
        s = sum(dist)
        if s > 0
            dist ./= s
        end
    end

    return dist
end

function mean_pauli_weight(v::AbstractVector, ::Val{N}; drop_identity=false) where {N}
    dist = pauli_weight_distribution(v, Val(N); normalize=true,
                                     drop_identity=drop_identity)

    return sum(w * dist[w + 1] for w in 0:N)
end

dist = pauli_weight_distribution(vt, Val(Nqubits))
mw = mean_pauli_weight(vt, Val(Nqubits))

println("Weight distribution:")
for w in 0:Nqubits
    println("weight $w : ", dist[w + 1])
end

println("Mean Pauli weight = ", mw)