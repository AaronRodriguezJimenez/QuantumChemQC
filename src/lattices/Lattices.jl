"""
 This module generates lattice Hamiltonians and allowt to plot them 
 in a graph representation.

 Main IDEAS:
 A lattice in physics is a collection of discrete sites with defined neighbors
 interactions.
 Meanwhile, a graph is a collection of vertices and edges connecting them.

 The grap can be considered to be a data structure that encoded lattice connectivity
 independent of the geometry.

 - Each site/orbital in the lattice <-> vertex in the graph.
 - Each hopping path/bond/interaction <-> edge in the graph.
 - Each edge weight <-> interaction strength (hopping integral, exchange coupling, etc).
"""
module Lattices

using Graphs
using GraphPlot
using Karnak
using NetworkLayout
using Colors

include("./lattice_types.jl")
include("./lattice_to_graph.jl")
include("./lattice_collections.jl")

export LatticeBond, Lattice# generate_lattice, hamiltonian, plot_lattice
export lattice2graph, dec2bin, bin2dec,bin2bonds!
export square_lattice

end # module Lattices

