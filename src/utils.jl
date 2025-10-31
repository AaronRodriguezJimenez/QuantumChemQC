"""
 # # #  General Utilities
"""

function nuclear_repulsion_energy(bset, p)

    E_nuc = 0.0
    for i = 1:p.natoms
        for j = i+1:p.natoms
            E_nuc +=
                bset.atoms[i].Z * bset.atoms[j].Z /
                (norm(bset.atoms[i].xyz - bset.atoms[j].xyz) * 1.88973)
        end
    end

    return E_nuc

end

function Z_to_symbol(Z)
    dict = Dict(
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        10 => "Ne",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        21 => "Sc",
        22 => "Ti",
        23 => "V",
        24 => "Cr",
        25 => "Mn",
        26 => "Fe",
        27 => "Co",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        31 => "Ga",
        32 => "Ge",
        33 => "As",
        34 => "Se",
        35 => "Br",
        36 => "Kr",
        37 => "Rb",
        38 => "Sr",
        39 => "Y",
        40 => "Zr",
        41 => "Nb",
        42 => "Mo",
        43 => "Tc",
        44 => "Ru",
        45 => "Rh",
        46 => "Pd",
        47 => "Ag",
        48 => "Cd",
        49 => "In",
        50 => "Sn",
        51 => "Sb",
        52 => "Te",
        53 => "I",
        54 => "Xe",
        83 => "Bi",
    )

    return dict[Z]
end

function coeff_clip!(ps::PauliSum{N}; thresh=1e-16) where {N}
    filter!(p->abs(p.second) > thresh, ps)
end

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