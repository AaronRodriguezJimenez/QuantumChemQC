using QuantumChemQC
using LinearAlgebra

mol = """
0 1
 H  0.0000   0.751  -0.485
 H  0.0000  -0.751  -0.485
"""

bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = true)

S, T, V, H, I = QuantumChemQC.compute_integrals(bset, p)

println("- - - - S matrix - - - ")
display(S)

println("- - - - T matrix - - - ")
display(T)

println("- - - - V matrix - - - ")
display(V)

println("- - - - H matrix - - - ")
display(H)

println("- - - - I matrix - - - ")
display(I)