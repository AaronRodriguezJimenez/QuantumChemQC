using PauliOperators
using QuantumChemQC
using Tullio

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
@tullio  mo_eris[p,q,r,s] := C[μ,p] * C[ν,q] * C[λ,r] * C[σ,s] * ao_eris[μ,ν,λ,σ]


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
core =0 #In this case we do not take an AS so Ecore = 0
h0 = scf_obj.Enuc + core 
h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)

"""
 Now we can construc the fermionic hamiltonian as
"""
# One-body terms: a†_p a_q
println("  ")
fzero = FermionOperator(ComplexF64(h0))
display(fzero)

n = size(h1, 1)
for p in 1:n
    for q in 1:n
        coeff = h1[p, q]
        if abs(coeff) > 1e-7
            term = [(p, true), (q, false)]  # 1-based indexing
            #  println("a$p a$q $term $coeff")
            fop = FermionOperator(term, coeff)
            display(fop)
        end
    end
end


# Two-body terms: a†_p a†_q a_r a_s
for p in 1:n, q in 1:n, r in 1:n, s in 1:n
    coeff = h2[p, q, r, s]
    if abs(coeff) > 1e-7
        term = [(p, true), (q, true), (r, false), (s, false)]
        #  println("a$p a$q a$r a$s $term $coeff")
        fop = FermionOperator(term, coeff)
        display(fop)
    end
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

println(" Qubit Operator Paulis")
for term in qubit_paulis
    display(term)
end

for coeff in qubit_coeffs
    display(coeff)
end
println("Total Terms are : ", length(qubit_coeffs))