using QuantumChemQC
using LinearAlgebra

mol = """
0 1
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

# Initial example: computing molecular integrals for a given molecule.
# bset contains the information of the basis BasisSet
# p contains the environment parameters
bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = true)

S, T, V, H, I = QuantumChemQC.compute_integrals(bset, p)

# Overlap Matrix
println("- - - - S matrix - - - ")
display(S)

# Kinetic Energy Matrix
println("- - - - T matrix - - - ")
display(T)

# Nuclear-Electron Matrix
println("- - - - V matrix - - - ")
display(V)

# H = T+V Matrix
println("- - - - H matrix - - - ")
display(H)

# ERIs Matrix
println("- - - - I matrix - - - ")
display(I)