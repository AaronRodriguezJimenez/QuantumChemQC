using PauliOperators
using QuantumChemQC
using Tullio
using Printf
using IterTools
using LinearAlgebra

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
bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = false)
scf_obj = QuantumChemQC.SCF(mol, bset, p)

"""
 We have generated an scf object, which contains the information for the molecule
 obtained after a RHF SCF procedure, analogous to PySCF, we can access to its content.
 Our first task is to construct the fermionic operator.
"""
# Atomic integrals
ao_hcore = scf_obj.T + scf_obj.V
ao_eris = scf_obj.I

# Transform Atomic integrals into Molecular integrals
C = scf_obj.C # MO coeffs

"""
 Molecular orbitals are linear combinations of atomic orbitals forming an
 eigenbasis for the Fock operator. In order to transform the atomic integrals
 into molecular integrals, we need to transform them using the matrix which 
 describes those linear combinations which is composed by the MO coeffs.
 
 Working in the molecular orbital basis offer various advantages, from a practical
 point of view: the natural choice of reference state is expressed in this basis,
 and additionally the easiest Fermion to Qubit mapping works in this representation.
"""

mo_hcore = C' * ao_hcore * C 
mo_eris = Array{Float64, 4}(undef, size(ao_eris)...)
# Openfermion convention:
@tullio mo_eris[p,q,r,s] := C[μ, p] * C[λ, s] * C[ν, q] * C[σ, r] * ao_eris[μ, λ, ν, σ]

# If needed, we can convert between chemist and physicis notations by transposing axes:
#ao_eris_phys = permutedims(ao_eris_chem, (1,3,2,4))  # chemist → physicist
#println("One body integrals (H_core) in MO representation: ")
#display(mo_hcore)
#println("Two body integrals in MO representation: ")
#display(mo_eris)

"""
 In some cases depending on the storing of the two electron integrals, some codes 
 may take a prefered order different from others, that's the case for example when
 comparing the results from openfermion and those with pyscf. 
 Here we just wan to figure out how to proced and do not take the order in consideration.

 Now we can continue with the construction of the fermionic Hamiltonian.
 Here one can actually choose an active space, by selecting those occupied and
 virtual orbitals of interest. If there's interest in this, one needs to "cut" the 
 obtained integrals according to the selected space. 
 Here we do not take this in consideration.

 The matrices so far contain only the information of the spatial orbitals, we need
 to introduce the information of the spin in order to construct proper spinorbitals
 that will be associated with qubits.

 So, the next step is to take the molecular integrals and construct the coefficients
 for each one- and two-body spin orbital interaction appearing in the second quantized
 hamiltonian. This is nothing but an antisymmetrization procedure which will give 
 the spin-orbital TENSORS as a result:
"""
core =0 #In this case we do not take an Active Space so Ecore = 0
h0 = scf_obj.Enuc + core 
h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)
# Note: current implementation does not yield the interleaved version of the fermionic 
# hamiltonian. Matching with the default version from OpenFermion.

"""
 Now we can construc the fermionic hamiltonian as
"""
H_fermion1, H_fermion2 = QuantumChemQC.fermionic_hamiltonian(h0, h1, h2)

println(" ")
println(" The Fermionic operator is")
println(" ")
for term in H_fermion1
    display(term)
end

for term in H_fermion2
    display(term)
end

"""
 This is the Fermionic Hamiltonian, in order to transform to a qubit Hamiltonian
 we can recurr to the Jordan-Wigner transformation.
 In order to do this we will make use of the PauliOperators.jl package,
 which recurr to the symplectic representation of Pauli Strings.
 
 For uses into Pauli Propagation techniques it is more usefull to represent the 
 Hamiltonian as a function that delivers the operators and the associated coefficients
 more than delivering a single Qubit operator object.

"""

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
 Compute Hartree-Fock expectation value with the qubit Hamiltonian
 <ket_0 | H_qubit | ket_0>
"""
function compute_expectation(ham_dict, ket, bra)
    expval = 0.0 + 0.0im
    for (pauli, coeff) in ham_dict
        v_new = pauli * ket
        #display(v_new)
        braket = bra * v_new[2]
        expval += coeff * v_new[1] * braket
    end
    return expval
end

ket = Ket(N, 3) #HF STATE
bra = ket' 
expectation = compute_expectation(ham_dict, ket, bra)
println(" ")
println(" HF initial state:")
display(ket)
println(" ")
println("Expectation value $expectation")