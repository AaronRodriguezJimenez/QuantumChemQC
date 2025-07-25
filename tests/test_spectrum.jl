using PauliOperators
using QuantumChemQC
using Tullio
using Printf
using IterTools
using LinearAlgebra
using Printf

"""
 Here we explore the construction of a qubit hamiltonian for H2
 under the Jordan-Wigner transformation.
"""

mol = """
0 1
 H  0.0000   0.000  0.000
 H  0.0000   0.000  0.740
"""

#Iitialize molecule
# bset - basis set object from GaussianBasis.
# p    - ENV parameters used for the Quantum Chemistry functions.
bset, p = QuantumChemQC.molecule(mol, "3-21g", spherical = false)
scf_obj = QuantumChemQC.SCF(mol, bset, p)

# Atomic integrals
ao_hcore = scf_obj.T + scf_obj.V
ao_eris = scf_obj.I

# Transform Atomic integrals into Molecular integrals
C = scf_obj.C # MO coeffs
mo_hcore , mo_eris = QuantumChemQC.ao2mo_coefficients(C, ao_hcore, ao_eris)

# Spinorbital tensors
core =0 
h0 = scf_obj.Enuc + core 
h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)
N = size(h1, 1)
qubit_paulis, qubit_coeffs = QuantumChemQC.qubit_hamiltonian(N, h0, h1, h2) 

#println("Total Terms are : ", length(qubit_paulis))
#println("Total Coeff are : ", length(qubit_coeffs))

println(" ")
println(" The Qubit operator is")
println(" ")

ham_dict = QuantumChemQC.combine_pauli_terms(qubit_paulis, qubit_coeffs)

for (pauli, coeff) in sort(collect(ham_dict); by=x->string(x[1]))
    println(@sprintf("%+.6f %+.6fim | %s", real(coeff), imag(coeff), string(pauli)))
end

println("Total number of operators: ", length(ham_dict))

"""
 Compute and print the eigenspectrum for the qubit hamiltonian
    
"""
function eigenspectrum_hop(Hdict; digits=10)
    # Build Hamiltonian matrix
    _, first_coeff = first(Hdict)
    first_pauli = first(keys(Hdict))
    dim = size(Matrix(first_pauli), 1)
    H_mat = zeros(ComplexF64, dim, dim)

    for (pauli, coeff) in Hdict
        H_mat .+= coeff .* Matrix(pauli)
    end

    # Compute and sort eigenvalues
    eigvals = sort(eigen(H_mat).values)

    println("  Eigenvalues of the Hamiltonian:")
    println("────────────────────────────────────")

    # Build format strings dynamically
    real_fmt_str = "λ%-2d = % .$(digits)f\n"
    complex_fmt_str = "λ%-2d = % .$(digits)f %+ .$(digits)fi\n"

    for (i, val) in enumerate(eigvals)
        buf = IOBuffer()
        if abs(imag(val)) < 1e-10
            Printf.format(buf, Printf.Format(real_fmt_str), i, real(val))
        else
            Printf.format(buf, Printf.Format(complex_fmt_str), i, real(val), imag(val))
        end
        print(String(take!(buf)))
    end

    println("────────────────────────────────────")
    return eigvals;
end


# Example usage
println(" ")
eigvals = eigenspectrum_hop(ham_dict)
