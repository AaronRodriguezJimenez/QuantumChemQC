using .Lattices
using Graphs, GraphPlot


"""
 Test for lattice types
"""
# Create some sites
bond_1 = LatticeBond(1, 2, 0.0, 0.0, 0.0, 0.0, "")
bond_2 = LatticeBond(2, 3, 1.0, 0.0, 1.0, 1.0, "type1")
bond_3 = LatticeBond(3, 1, 0.0, 1.0, 0.0, 0.0, "type2") # triangular bond  
println("bond_1: ", bond_1)
println("bond_2: ", bond_2)
println("bond_3: ", bond_3)

# Create a lattice
lattice = Lattice([bond_1, bond_2, bond_3])
println("Lattice: ", lattice)

# Creat a square_lattice
sq_latt = square_lattice(3, 3, yperiodic=false)
println("Square Lattice (3x3): ")
for bond in sq_latt
    println(bond)
end

# Convert to graph
g = lattice2graph(sq_latt)
draw(PDF("gSquare.pdf", 6inch, 6inch), gplot(g, nodefillc=nodefillcX, nodelabel=1:nv(g), title="Square Lattice"))
