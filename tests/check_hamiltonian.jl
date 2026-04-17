using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots

function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/benzene_and_acenes/benz-RHF_integrals.npz"
    data = npzread(data_path)
    H0 = data["hc"][1]
    println("H0: ", H0)
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

# Helper: extract PauliBasis operators and coefficients from PauliSum
function extract_hamiltonian_coeffs_and_ops(H::PauliSum{N,T}) where {N,T}
    ops = PauliBasis{N}[]
    coeffs = Float64[]
    for (p, c) in H
        push!(ops, p)
        push!(coeffs, real(c))
    end
    return ops, coeffs
end

"""
- Operator atributeds o is the operator to be evolved under the hamiltonian H
"""
function check_ops(H::PauliSum{N,T}; threshold=1e-6) where {N,T}

    # Extract Pauli strings and coefficients
    H_ops, H_coeffs = extract_hamiltonian_coeffs_and_ops(H)
    println("Hamiltonian H has $(length(H_coeffs)) terms")

    # Count small coefficients
    small_H = count(abs.(H_coeffs) .< threshold)

    println("Terms in H with |coeff| < $threshold: $small_H")

    return (
        H_ops = H_ops,
        H_coeffs = H_coeffs,
        small_H = small_H,
    )
end

H = run_H()
check_ops(H)
Hmat = Matrix(H)
ket, _ = QuantumChemQC.string_to_ket("111000111000")
V0 = Vector(ket)

#1) Check reference state energy:
Ek = V0' * Hmat * V0
println("Reference state energy: ", Ek)
