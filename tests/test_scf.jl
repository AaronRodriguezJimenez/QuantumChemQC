using QuantumChemQC


println("= = = = = = = = = = = = = = = = = = = = =")
println("Testing SCF module")
println(" ")
println("Testing H2 molecule at 0.74 Angstroms")
println("= = = = = = = = = = = = = = = = = = = = =")

mol = """
0 1
 H  0.0000   0.000  0.000
 H  0.0000   0.000  0.740
"""

#Iitialize molecule
bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = false)
println(typeof(bset))
println(typeof(p))
#Run an scf and get an SCF object:

scf_obj = QuantumChemQC.SCF(mol, bset, p)
println("Total energy is:", QuantumChemQC.total_energy(scf_obj))

println(" ")
println("Testing H8 cube at 0.74 Angstroms")
println("= = = = = = = = = = = = = = = = = = = = =")

# Build a cube of H atoms using the Molecules module
mol2 = Molecules.build_h_cube_mol(0.74, 2)  # 2x2x2 cube of H atoms

println(mol2)
#Iitialize molecule

#Initialize molecule
bset2, p2 = QuantumChemQC.molecule(mol2, "sto-3g", spherical = false)

#Run an scf and get an SCF object:
scf_obj2 = QuantumChemQC.SCF(mol2, bset2, p2)
println("Total energy is:", QuantumChemQC.total_energy(scf_obj2))

println(" ")
println("Testing H64 cube at 0.74 Angstroms")
println("= = = = = = = = = = = = = = = = = = = = =")

# Build a cube of H atoms using the Molecules module
mol3 = Molecules.build_h_cube_mol(0.74, 4)  # 4x4x4 cube of H atoms
println(mol3)
#Iitialize molecule

#Initialize molecule
bset3, p3 = QuantumChemQC.molecule(mol3, "sto-3g", spherical = false)

#Run an scf and get an SCF object:
scf_obj3 = QuantumChemQC.SCF(mol3, bset3, p3)
println("Total energy is:", QuantumChemQC.total_energy(scf_obj3))
