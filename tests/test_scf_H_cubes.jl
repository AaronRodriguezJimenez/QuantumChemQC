using QuantumChemQC


println("= = = = = = = = = = = = = = = = = = = = =")
println("Testing SCF module")
println(" ")
println("Testing H (2x2x2) cube at 0.74 Angstroms")
println("= = = = = = = = = = = = = = = = = = = = =")

# Build a cube of H atoms using the Molecules module
mol3 = Molecules.build_h_cube_mol(0.74, 2)  # 2x2x2
println(mol3)
#Iitialize molecule

#Initialize molecule
bset3, params = QuantumChemQC.molecule(mol3, "sto-3g", spherical = false)

#Run an scf and get an SCF object:
scf_obj3 = QuantumChemQC.SCF(mol3, bset3, params)
println("Total energy is:", QuantumChemQC.total_energy(scf_obj3))

# = = = Qubit Hamiltonian = = =
H, N = QuantumChemQC.get_qubit_hamiltonian(scf_obj3)
println("Qubit Hamiltonian (N=$N qubits):")
display(H)


println("= = = = = = = = = = = = = = = = = = = = =")
println("Testing SCF module")
println(" ")
println("Testing H (3x3x3) cube at 0.74 Angstroms")
println("= = = = = = = = = = = = = = = = = = = = =")

# Build a cube of H atoms using the Molecules module
mol4 = Molecules.build_h_cube_mol(0.74, 4)  # 4x4x4
println(mol4)
#Iitialize molecule

#Initialize molecule
bset4, p4 = QuantumChemQC.molecule(mol4, "sto-3g", spherical = false)

#Run an scf and get an SCF object:
scf_obj4 = QuantumChemQC.SCF(mol4, bset4, p4)
println("Total energy is:", QuantumChemQC.total_energy(scf_obj4))

# = = = Qubit Hamiltonian = = =
H, N = QuantumChemQC.get_qubit_hamiltonian(scf_obj4)
println("Qubit Hamiltonian (N=$N qubits):")
display(H)
