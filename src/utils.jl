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
 is based on openfermion.ops.representations.get_tensors_from_integrals, we follow  
 the indexing convention
 p = 2i      → β spin of orbital i
 p = 2i - 1  → α spin of orbital i
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


"""
 # # # Qubit utilities
"""

"""
 Under the Symplectic representation, operators are mapped as follows:
        a_j^dagger -> Z_0 .. Z_{j-1} (X_j - iY_j) / 2
        a_j        -> Z_0 .. Z_{j-1} (X_j + iY_j) / 2
"""
function JW_creator_mapping(N, i::Int64)
    # Compute a^dagger_i term
    ax_term = Pauli(2^(i-1)-1, 2^(i-1), N)
    ay_term = Pauli(2^(i)-1, 2^(i-1), N)
    a_dagg_i = 0.5 * (ax_term - ay_term)
    return a_dagg_i
end

function JW_annihilator_mapping(N, j::Int64)
    # Compute a_j term
    bx_term = Pauli(2^(j-1)-1, 2^(j-1), N)
    by_term = Pauli(2^(j)-1, 2^(j-1), N)
    a_j = 0.5 * (bx_term + by_term)
    return a_j
end

"""
 The following is the generator function that can produce the qubit Hamiltoinan
    in the JW representation, as input it uses the resulting molecular orbitals
    from a SCF calculation expressed in MOs.
"""
function qubit_hamiltonian(N::Int64, h0::Float64, h1::Matrix{Float64}, h2::Array{Float64, 4})
 
    generators = Vector{Pauli{N}}()
    parameters = Vector{ComplexF64}()

    one_e_term = PauliSum(N)
    two_e_term = PauliSum(N)

    # One-body terms: a†_p a_q
    n = size(h1, 1)
    for p in 1:n
        for q in 1:n
            coeff = h1[p, q]
            a_dag_i = JW_creator_mapping(N, p)
            a_j = JW_annihilator_mapping(N, q)
          #  display(a_dag_i*a_j)
          #  println("---")
            one_e_term += coeff*(a_dag_i * a_j)
        end
    end
    #println(length(one_e_term))

    for (pauli, coeff) in one_e_term
        abs(coeff) > 1e-7 || continue
        push!(generators, Pauli(pauli))
        push!(parameters, coeff)
        println("Coeff associated to paulis in one_e_term $coeff")
        println(Pauli(pauli))
    end

    # Two-body terms: a†_p a†_q a_r a_s
    for p in 1:n, q in 1:n, r in 1:n, s in 1:n
        # Apply physicist's antisymmetrization
        coeff = h2[p,q,r,s] 
        a_dag_p = JW_creator_mapping(N, p)
        a_dag_q = JW_creator_mapping(N, q)
        a_r = JW_annihilator_mapping(N, r)
        a_s = JW_annihilator_mapping(N, s)
        two_e_term += coeff* (a_dag_p * a_dag_q * a_r * a_s)
    end

    #println(length(two_e_term))
    for (pauli, coeff) in two_e_term
        abs(coeff) > 1e-7 || continue
        push!(generators, Pauli(pauli))
        push!(parameters, coeff)
        println("Coeff associated to paulis in two_e_term $coeff")
        println(Pauli(pauli))
    end

    return generators, parameters
end