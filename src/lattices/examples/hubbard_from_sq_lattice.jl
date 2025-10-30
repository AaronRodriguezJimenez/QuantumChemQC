using QuantumChemQC

"""
 Construct a Hubbard model from a square lattice.
"""
Nx = 2
Ny = 2

#sites = siteinds("Electron", N; conserve_qns = true)

lattice = square_lattice(Nx, Ny; yperiodic = false)

println("Square lattice with ", length(lattice), " bonds:")
for bond in lattice
    println(bond)
end
