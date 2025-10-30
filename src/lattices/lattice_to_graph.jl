using Graphs, GraphPlot


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
  Function for transforming from binary to decimal representation
"""
function bin2dec(b::AbstractVector{<:Integer})
    n = 0
    N = length(b)
    for i in 1:N
        n += b[N-i+1] * (2^(i-1))
    end
    return n
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