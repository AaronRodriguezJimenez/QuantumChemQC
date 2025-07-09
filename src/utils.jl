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

"""
 This function converts the one and two-body integrals (in MO basis) into tensor form
 is based on openfermion.ops.representations.get_tensors_from_integrals
"""
function get_spin_orbital_tensors(hcore::Matrix{Float64}, eri::Array{Float64,4})
    n_orb = size(hcore, 1)
    n_spin = 2 * n_orb

    h1 = zeros(n_spin, n_spin)
    h2 = zeros(n_spin, n_spin, n_spin, n_spin)

    for p in 1:n_orb
        for q in 1:n_orb
            # Populate 1-body terms
            h1[2p-1, 2q-1] = hcore[p, q]  # alpha-alpha
            h1[2p,   2q]   = hcore[p, q]  # beta-beta

            for r in 1:n_orb
                for s in 1:n_orb
                    val = eri[p, q, r, s] / 2

                    # Mixed spin
                    h2[2p-1, 2q,   2r,   2s-1] += val
                    h2[2p,   2q-1, 2r-1, 2s]   += val

                    # Same spin
                    h2[2p-1, 2q-1, 2r-1, 2s-1] += val
                    h2[2p,   2q,   2r,   2s]   += val
                end
            end
        end
    end

    return h1, h2
end
