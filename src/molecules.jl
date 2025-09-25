
"""
 Selected Molecular models to be used as toy models for testing and demonstration.
"""

module Molecules

using Printf
export build_h_cube_mol

"""
 build_h_cube_mol(a::Float64, n::Int) -> String

    Function to build a cube of hydrogen atoms with side length a (in Angstroms)
    returns a string that can be used as input to QuantumChemQC.molecule

    a -> distance between atoms in Angstroms.
    n -> number of atoms along each edge of the cube.
"""
function build_h_cube_mol(a::Float64, n::Int)
    lines = ["0 1"]  # charge and multiplicity
    for i in 0:n-1
        for j in 0:n-1
            for k in 0:n-1
                x, y, z = i*a, j*a, k*a
                push!(lines, @sprintf(" H  %8.4f  %8.4f  %8.4f", x, y, z))
            end
        end
    end
    return join(lines, "\n")
end

end # module Molecules