using PauliOperators
using QuantumChemQC
using Tullio
using Printf
using IterTools

"""
 Here we test the functions designed to check the properties of the 2-body
 operators.
"""

mol = """
0 1
 H  0.0000   0.000  0.000
 H  0.0000   0.000  0.740
"""

#Iitialize molecule
bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = false)
scf_obj = QuantumChemQC.SCF(mol, bset, p)

# Atomic integrals
ao_hcore = scf_obj.T + scf_obj.V
ao_eris = scf_obj.I

# Transform Atomic integrals into Molecular integrals
C = scf_obj.C # MO coeffs

mo_hcore = C' * ao_hcore * C 
mo_eris = Array{Float64, 4}(undef, size(ao_eris)...)
# Note ao_eris are in chemist notation, so the following line:
@tullio  mo_eris[p,q,r,s] := C[μ,p] * C[ν,q] * C[λ,r] * C[σ,s] * ao_eris[μ,ν,λ,σ]
# will produce discrepancies with the result from openfermion
# To fix it, we need to perform:
#@tullio mo_eris[p,q,r,s] := C[μ,p] * C[λ,r] * C[ν,q] * C[σ,s] * ao_eris[μ,λ,ν,σ]
mo_eris = QuantumChemQC.chem_to_phys(mo_eris)



#The following matches with openfermion ordering
println("One body integrals (H_core) in AO representation: ")
display(ao_hcore)
println("Two body integrals in AO representation: ")
display(ao_eris)
println("One body integrals (H_core) in MO representation: ")
display(mo_hcore)
println("Two body integrals in MO representation: ")
display(mo_eris)

# Spinorbitals tensors:
core =0 #Not using active space
h0 = scf_obj.Enuc + core 
h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)


# Thesting the tensors to TEST

println("Test permutation symmetries for the two-body tensor")
sym_hold = QuantumChemQC.all_permutation_symmetries_hold(ao_eris)
type_eris = QuantumChemQC.find_index_order(ao_eris)
println("- Atomic: $sym_hold, type: $type_eris")

sym_hold = QuantumChemQC.all_permutation_symmetries_hold(mo_eris)
type_eris = QuantumChemQC.find_index_order(mo_eris)
println("- Molecular: $sym_hold, type: $type_eris")

sym_hold = QuantumChemQC.all_permutation_symmetries_hold(h2)
type_eris = QuantumChemQC.find_index_order(h2)
println("- Spinorbital: $sym_hold, type: $type_eris")


