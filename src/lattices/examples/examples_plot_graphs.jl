using Graphs
using GraphPlot
using Karnak
using NetworkLayout
using Colors
g = barabasi_albert(60, 1)

#g = Graph(3) # graph with 3 vertices
# make a triangle
#add_edge!(g, 1, 2)
#add_edge!(g, 1, 3)
#add_edge!(g, 2, 3)

# Using GraphPlot: simple plotting of graphs
#gplot(g, nodelabel=1:60)


# Plots with Karnak: The Karnak.jl package integrates the Luxor.jl 2D graphics package, and uses NetworkLayout.jl for calculating layouts. 
@drawsvg begin
    background("black")
    sethue("grey40")
    fontsize(8)
    drawgraph(g, 
        layout=stress, 
        vertexlabels = 1:nv(g),
         vertexfillcolors = 
            [RGB(rand(3)/2...) 
               for i in 1:nv(g)]
    )
end 600 400