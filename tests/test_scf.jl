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

