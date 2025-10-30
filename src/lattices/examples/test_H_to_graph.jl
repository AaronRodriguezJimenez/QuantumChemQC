using Revise
using QuantumChemQC
using Compose
using Graphs, GraphPlot

H, N = QuantumChemQC.H2_hamiltonian(0.74)
println("H2 Hamiltonian at d=0.74 Angstrom:")
display(H)

"""
  Function for transforming form decimal to binary representation
"""
function dec2bin(n::Union{Int, Int128}, N::Union{Int, Int64})
    b = zeros(Int, N)
    for i in 1:N
        b[N-i+1] = n % 2
        n = n ÷ 2
    end
    return b
end

"""
 Function for trasforming the binary represetation to
 to LatticeBond
 """
function bin2bonds!(acc::Vector{LatticeBond}, b::AbstractVector{<:Integer}; mode::Symbol=:path)
    idx = findall(!iszero, b)
    n = length(idx)
    n < 2 && return acc

    if mode === :path
        @inbounds for k in 1:n-1
            push!(acc, LatticeBond(idx[k], idx[k+1]))
        end
    elseif mode === :clique
        @inbounds for i in 1:n-1, j in i+1:n
            push!(acc, LatticeBond(idx[i], idx[j]))
        end
    elseif mode === :star
        c = idx[1]
        @inbounds for j in 2:n
            push!(acc, LatticeBond(c, idx[j]))
        end
    else
        throw(ArgumentError("unknown mode=$mode"))
    end
    return acc
end

"""
 Function for transforming from LatticeBond to a graph from Graphs
"""
function lattice2graph(lattice::Lattice)
    g = Graph()
    # First, find the maximum site number to determine number of vertices
    max_site = 0
    for bond in lattice
        max_site = max(max_site, bond.s1, bond.s2)
    end
    # Add vertices
    for i in 1:max_site
        add_vertex!(g)
    end
    # Add edges based on bonds
    for bond in lattice
        add_edge!(g, bond.s1, bond.s2)
    end
    return g
end


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

draw(PDF("gX.pdf", 6inch, 6inch), gplot(gX, layout=shell_layout, nodelabel=1:nv(gX), title="X Lattice"))
draw(PDF("gZ.pdf", 6inch, 6inch), gplot(gZ, layout=shell_layout, nodelabel=1:nv(gZ), title="Z Lattice"))