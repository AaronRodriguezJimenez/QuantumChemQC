using Revise
using QuantumChemQC
using Compose
import Colors
using Graphs, GraphPlot


lx = 2
ly = 1
N = lx*ly  # total number of sites
Jx = 1.0
Jy = 0.0
Jz = 1.0
H= QuantumChemQC.heisenberg_1D(N, Jx, Jy, Jz; x=0, y=0, z=0)

##### Lets convert the PauliSum H to a lattice #####
# Create X_lattice and Z_lattice for X and Z terms in H and plot them
X_lattice = Lattice()
Z_lattice = Lattice()

display(H)

for (p, c) in H
#for p in gen
    display(p)
    println(p.z, " ", p.x)
    b = bin2bonds!(Vector{LatticeBond}(), dec2bin(p.x, N))
    println("X -> bond: ", b)
    for lb in b
        push!(X_lattice, lb)
    end
    b = bin2bonds!(Vector{LatticeBond}(), dec2bin(p.z, N))
    println("Z -> bond: ", b)
    for lb in b
        push!(Z_lattice, lb)
    end
 #   println("Coefficient: ", c)
    println("-----")
end
println("Lattice with ", length(X_lattice), " X bonds and ", length(Z_lattice), " Z bonds:")
gX = lattice2graph(X_lattice)
for bond in X_lattice
    println("X bond: ", bond)
end

gZ = lattice2graph(Z_lattice)
for bond in Z_lattice
    println("Z bond: ", bond)
end

println("X Graph has ", nv(gX), " vertices and ", ne(gX), " edges.")
println("Z Graph has ", nv(gZ), " vertices and ", ne(gZ), " edges.")

nodefillcX = Colors.colorant"skyblue1"
nodefillcZ = Colors.colorant"darkseagreen"

draw(PDF("Heis_gX.pdf", 6inch, 6inch), gplot(gX, nodefillc=nodefillcX, layout=shell_layout, nodelabel=1:nv(gX), title="X Lattice"))
draw(PDF("Heis_gZ.pdf", 6inch, 6inch), gplot(gZ, nodefillc=nodefillcZ, layout=shell_layout, nodelabel=1:nv(gZ), title="Z Lattice"))
#draw(PDF("gX.pdf", 6inch, 6inch), gplot(gX, layout=spectral_layout, nodelabel=1:nv(gX), title="X Lattice"))
#draw(PDF("gZ.pdf", 6inch, 6inch), gplot(gZ, layout=spectral_layout, nodelabel=1:nv(gZ), title="Z Lattice"))

println("Fermi-Hubbard Hamiltonian on a $lx x $ly lattice:")
display(H)