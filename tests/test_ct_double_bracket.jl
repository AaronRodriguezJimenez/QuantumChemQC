using PauliOperators
using UnitaryPruning
using QuantumChemQC
using Tullio
using Printf
using IterTools
using LinearAlgebra
using Printf

"""
 Here we test the idea for the computation of the Correlation function
    C(t) = <psi|exp(-iHt)|psi> using double bracket notation.
    psi is a HF state in the qubit bitstring notation.
    H is the qubit Hamiltonian in PauliOperator format.
    t is a time parameter.

    It can be shown that:
    C(t) = <psi|exp(-iHt)|psi> = exp(-iϕt) + 2(<0|[A†, [exp(-iH't),A]]|0>)
    
    where ϕ is the eigenvalue from H|0> = ϕ|0>
"""


function eval_ct(N, t, H, A,phi)
    H_mat = Matrix(H)
    A_mat = Matrix(A)

    # Compute exp(-iHt)
    expmH = exp(-im*H_mat*t)
    comm1 = expmH * A_mat - A_mat * expmH
    comm2 = A_mat' * comm1 - comm1 * A_mat'

    zero_state = Ket(N,0)
    display(zero_state)
    V0 = Vector(zero_state)
    expval = V0' * comm2 * V0
    val = exp(-im*phi*t) + 2*expval[1]
    return val
end

# Define a simple test Hamiltonian and operator A
N = 2  # Number of qubits
H = 0.5 * Pauli(N, Z=[1]) + 0.3 * Pauli(N, Y=[1]) * Pauli(N, X=[2])
A = Pauli(N, X=[1])
phi = 0.5  # Example eigenvalue
t = 1.0  # Time parameter
ct_value = eval_ct(N, t, H, A, phi)


println("C(t) value:", ct_value)
println("Magnitude of C(t):", abs(ct_value))
println("Phase of C(t):", angle(ct_value))
println("Norm of C(t):", norm(ct_value))
println("Real part of C(t):", real(ct_value))
println("Imaginary part of C(t):", imag(ct_value))
