using QuantumChemQC
using Printf

geom_H2 = Molecules.H2_geometry(0.74)
println("Molecule Geometry: ", geom_H2)

geom_H2 = [
    "H" => (0.0, 0.0, 0.0),
    "H" => (0.0, 0.0, 0.74),
]

mol_H2 = Molecules.Molecule(geom_H2; skip_FCI=false)
println("Molecule Info: ")
println(" Number of spatial orbitals: ", mol_H2.N)
println(" Number of spin orbitals:    ", mol_H2.n)
println(" Number of electrons:        ", mol_H2.Ne)
println(" Hartree-Fock Energy:        ", mol_H2.E_HF)
println(" FCI Energy:                 ", mol_H2.E_FCI)  # 

println("One Body Integrals h1: \n")
println(mol_H2.h1)
println("Two Body Integrals h2: \n")
println(mol_H2.h2)

H_mol = Molecules.molecular_hamiltonian(mol_H2)
println("Molecular Hamiltonian: \n")
println(H_mol)

# Construct the fermionic Hamiltonian components
N = mol_H2.n
h0 = mol_H2.h0
h1 = mol_H2.h1
h2 = mol_H2.h2

# Here we construct the qubit Hamiltonian using Jordan-Wigner mapping
# this can be done more efficiently if instrad we construct the PauliSum Hamiltonian directly
# and clip it.
#qubit_paulis, qubit_coeffs = QuantumChemQC.qubit_hamiltonian(N, h0, h1, h2) 
H  = QuantumChemQC.PauliSum_hamiltonian(N, h0, h1, h2)
println(" ")
println(" The Qubit operator is")
println(" ")
display(H)

println("Total number of operators: ", length(H))