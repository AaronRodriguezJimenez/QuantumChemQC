using Revise
using QuantumChemQC
using Compose
import Colors
using Graphs, GraphPlot
using JLD2
using Printf

lx = 7
ly = 7
N = 2*lx*ly  # total number of spin-orbitals
t = 0.1
U = 0.09
H= QuantumChemQC.fermi_hubbard_2D_snake(lx, ly, t, U, snake_ordering=false)

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

layout=(args...)->spring_layout(args...; C=2.0, MAXITER=10)
#draw(PDF("gX.pdf", 6inch, 6inch), gplot(gX, nodefillc=nodefillcX, layout=shell_layout, nodelabel=1:nv(gX), title="X Lattice"))
#draw(PDF("gZ.pdf", 6inch, 6inch), gplot(gZ, nodefillc=nodefillcZ, layout=shell_layout, nodelabel=1:nv(gZ), title="Z Lattice"))
draw(PDF("gX.pdf", 6inch, 6inch), gplot(gX,  layout=layout, nodelabel=1:nv(gX), title="X Lattice"))
draw(PDF("gZ.pdf", 6inch, 6inch), gplot(gZ, layout=layout, nodelabel=1:nv(gZ), title="Z Lattice"))

"""
  Now we can use this graph representation to explore the connectivity of the evolved Hamiltonian 
"""

# remove common invisible/unwanted whitespace characters
function clean_path(s::AbstractString)
    s2 = strip(s)                                  # remove leading/trailing regular whitespace
    # replace NBSP with normal space, remove zero-width and BOM characters
    s2 = replace(s2, '\u00A0' => ' ')             # NO-BREAK SPACE -> space
    s2 = replace(s2, '\u200B' => "")              # ZERO WIDTH SPACE
    s2 = replace(s2, '\u200C' => "")              # ZERO WIDTH NON-JOINER
    s2 = replace(s2, '\uFEFF' => "")              # BOM / ZERO WIDTH NO-BREAK SPACE
    # normalize Unicode to NFC (helpful if mixing composed/decomposed chars)
    try
        s2 = Unicode.normalize(s2, :NFC)
    catch
        # some Julia installs may not have Unicode.normalize; ignore if unavailable
    end
    return s2
end

# diagnostic printer for debugging paths
function show_path_diagnostics(wf::AbstractString)
    println("repr(path) = ", repr(wf))
    println("abspath(path) = ", repr(abspath(wf)))
    println("isdir(parent)?, isfile(path)?")
    wf_abs = abspath(wf)
    parent_dir = dirname(wf_abs)
    println(" parent_dir repr: ", repr(parent_dir))
    println(" isdir(parent_dir): ", isdir(parent_dir))
    println(" isfile(path): ", isfile(wf_abs))
    # show each char with codepoint for parent_dir (only if not too long)
    println("\nParent dir characters (index, char, U+hex, isspace):")
    j = 1
    for c in collect(parent_dir)
        @printf("%3d: '%s' U+%04X  isspace=%s\n", j, c, Int(c), isspace(c))
        j += 1
    end
end

function load_H(workfile::AbstractString, force_iostream::Bool=true)
    # 1) clean the raw input path
    wf_clean = clean_path(workfile)
    wf_abs = abspath(wf_clean)
    parent_dir = dirname(wf_abs)

    println("Attempting to load Hamiltonian from: ", repr(wf_abs))
    println("Parent directory (repr): ", repr(parent_dir))

    # 2) existence checks and safer diagnostics
    if !isdir(parent_dir)
        println("Parent directory does not exist according to Julia (isdir=false).")
        println("Showing diagnostics to help find hidden characters:")
        show_path_diagnostics(workfile)

        println("\nSuggested terminal debug (copy/paste) to run in your shell:")
        println("  ls -la \"$(parent_dir)\"")
        println("Or list the parent of the parent (in case of hidden characters):")
        println("  ls -la \"$(dirname(parent_dir))\"")
        error("Parent directory not found: " * repr(parent_dir))
    end

    if !isfile(wf_abs)
        println("File not found at path above; listing files in parent directory:")
        for f in sort(readdir(parent_dir))
            println("  $f")
        end
        error("File does not exist: " * repr(wf_abs))
    end

    # 3) open and read H
    try
        if force_iostream
            jldopen(wf_abs, "r"; iotype=IOStream) do f
                println("Opened JLD2 with IOStream.")
                if haskey(f, "H")
                    return read(f, "H")
                else
                    error("Key 'H' not found in JLD2 file.")
                end
            end
        else
            jldopen(wf_abs, "r") do f
                println("Opened JLD2 (default iotype).")
                if haskey(f, "H")
                    return read(f, "H")
                else
                    error("Key 'H' not found in JLD2 file.")
                end
            end
        end
    catch err
        println("Failed to open/read JLD2 file: ", err)
        rethrow(err)
    end
end

# Read the Hamiltonian from file
workdir = "/Users/admin/VSCProjects/QuantumChemQC/src/lattices/examples/"
hamiltonian_filename = "H_evolved_U_0.09_threshold_0.01_wmax_2_wtype_0.jld2"  # no trailing space
#workfile = joinpath(workdir, hamiltonian_filename)
workfile = "/Users/admin/VSCProjects/QuantumChemQC/src/lattices/examples/H_evolved_U_0.09_threshold_0.01_wmax_2_wtype_0.jld2"

H_loaded = load_H(workfile, true)

# ===========================================================

for (p, c) in H_loaded
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

layout=(args...)->spring_layout(args...; C=2.0, MAXITER=100)
#draw(PDF("gX_evolved.pdf", 6inch, 6inch), gplot(gX, nodefillc=nodefillcX, layout=shell_layout, nodelabel=1:nv(gX), title="X Lattice"))
#draw(PDF("gZ_evolved.pdf", 6inch, 6inch), gplot(gZ, nodefillc=nodefillcZ, layout=shell_layout, nodelabel=1:nv(gZ), title="Z Lattice"))
draw(PDF("gX_evolved.pdf", 6inch, 6inch), gplot(gX, layout=layout, nodelabel=1:nv(gX), title="X Lattice"))
draw(PDF("gZ_evolved.pdf", 6inch, 6inch), gplot(gZ, layout=layout, nodelabel=1:nv(gZ), title="Z Lattice"))

println("INITIAL Fermi-Hubbard Hamiltonian on a $lx x $ly lattice:")
display(H)
println("Evolved Hamiltonian under DBF (Loaded Hamiltonian from file):")
display(H_loaded)