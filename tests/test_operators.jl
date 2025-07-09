using QuantumChemQC
using PauliOperators

"""
 Here we show the syntaxis from PauliOperators
"""

N =  4 # State a given number of Qubits

# Create a ket
ket = Ket(N, 1)
println(ket)

# Create a bra
bra = Bra(N, 1)
println(bra)

# Perform their product
braket =  bra*ket
println(braket)

# Or their dyad
dyad = ket*bra
println(dyad)

#Define a Pauli
p = Pauli(N, Z=[1])
println(p)

#Create sum of PauliOperators
q = Pauli(N, X=[1])
psumq = p + q
println(psumq)