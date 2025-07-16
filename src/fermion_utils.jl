"""
 # # #  Fermion Utilities
"""

"""
 This function converts the one and two-body integrals (in MO basis) into tensor form
(one- and two-body coefficients) is based on 
openfermion.ops.representations.get_tensors_from_integrals, 
we follow the indexing convention
 p = 2i      → α spin of orbital i
 p = 2i - 1  → β spin of orbital i
"""
function get_spin_orbital_tensors(hcore::Matrix{Float64}, eri::Array{Float64,4})
    n_orb = size(hcore, 1)
    n_spin = 2 * n_orb

    h1 = zeros(n_spin, n_spin)
    h2 = zeros(n_spin, n_spin, n_spin, n_spin)

    #logical mask, indexing by a boolean array selects elements at the indices where its values are true.
    mask(val) = val * (abs(val) > 1e-10)

    for p in 1:n_orb
        for q in 1:n_orb
            # Populate 1-body terms
            val = mask(hcore[p, q])
            h1[2p, 2q] = val #hcore[p, q]  # alpha-alpha
            h1[2p-1, 2q-1]   = val #hcore[p, q]  # beta-beta

            for r in 1:n_orb
                for s in 1:n_orb
                    val = mask(eri[p, q, r, s] / 2)

                    # Mixed spin
                    h2[2p,   2q-1, 2r-1, 2s]   += val
                    h2[2p-1, 2q,   2r,   2s-1] += val

                    # Same spin
                    h2[2p,   2q,   2r,   2s]   += val      # α-α-α-α
                    h2[2p-1, 2q-1, 2r-1, 2s-1] += val      # β-β-β-β
                end
            end
        end
    end

    return h1, h2
end

"""
 Permutation symmetries satisfied by a rank-4 tensor of real two-body
 integrals in chemist's index order as depictedin order of appearance
 in Molecular Electronic Structure Theory by Helgaker et. al (MEST)
"""
const ChemIndexPermutations = (
    PERM_1    = ((1, 2), (3, 4)),     # (0,1) → (2,3)       MEST(1.4.17)
    PERM_2_AB = ((1,), (2,)),         # (0,) → (1,)         MEST(1.4.38)
    PERM_2_AC = ((3,), (4,)),         # (2,) → (3,)         MEST(1.4.38)
    PERM_2_AD = ((1, 3), (2, 4)),     # (0,2) → (1,3)       MEST(1.4.38)
    PERM_3    = ((1, 2), (4, 3)),     # (0,1) → (3,2)       PERM_2_AB and PERM_1
    PERM_4    = ((1, 2, 3), (3, 4, 2)), # (0,1,2) → (2,3,1) PERM_2_AC and PERM_1 
    PERM_5    = ((1, 2, 3), (4, 3, 2))  # (0,1,2) → (3,2,1) PERM_2_AD and PERM_1
)

"""
Mimics numpy.moveaxis: Moves axes from `source` to `dest` in `tensor`.
"""
function moveaxes(tensor::AbstractArray, source::Tuple{Vararg{Int}}, dest::Tuple{Vararg{Int}})
    orig_axes = collect(1:ndims(tensor))

    # Map to 0-based internally like NumPy for clarity (if needed)
    # But we stick to 1-based here
    perm = copy(orig_axes)

    # Build new axis order by removing `source` and inserting at new positions
    # Start from compressed axes with `source` removed
    perm = filter(x -> x ∉ source, perm)

    # Normalize dest indices in case they go beyond current length
    insert_positions = dest

    # Insert sources in destination positions, adjusting for shift during insertions
    for (i, (s, d)) in enumerate(zip(source, insert_positions))
        idx = min(d, length(perm) + 1)  # cap insertion point
        insert!(perm, idx, s)
    end

    return permutedims(tensor, perm)
end


"""
 This function returns wether the provided tensor remains identical
 under the provided permutation.
 Args:
 - two_body_tensor: is the tensor to test 
 - permutation: The source and destination indices of the axis permutation 
"""
function _check_two_body_symmetry(
    tensor::Array{<:Real,4},
    src::Tuple{Vararg{Int}},
    dst::Tuple{Vararg{Int}}; 
    rtol::Float64 = 1e-5, 
    atol::Float64 = 1e-8
) :: Bool
    permuted_tensor = moveaxes(tensor, src, dst)
    return isapprox(tensor, permuted_tensor; rtol=rtol, atol=atol)
end

"""
 This function return whether a tensor has the required symmetries to represent two-electron terms.
    returns true if the whole rank-4 tensor has the required symmetries for coefficients of the two-electron terms.
"""
function all_permutation_symmetries_hold(two_body_tensor::Array{Float64,4}; rtol=1e-5, atol=1e-8) :: Bool
    #rtol = 1e-5  #Numerical tolerance for relative comparison
    #atol = 1e-8  #Numerical tolerance for absolute comparison
    for (name, (src, dst)) in pairs(ChemIndexPermutations)
        if !_check_two_body_symmetry(two_body_tensor, src, dst; rtol=rtol, atol=atol)
            println("Failed symmetry: ", name)
            return false
        end
    end
    return true
end

# 
# - - - Pending, add:
# transformations between physics and chemist 
# symmetrization
@enum IndexType begin
    CHEMIST
    PHYSICIST
    INTERMEDIATE
    UNKNOWN
end

"""
 Convert the rank-four tensor representing the two-body integrals from physicists' index
 to chemists' index order: i,j,k,l -> i,l,j,k
"""
function phys_to_chem(tensor::Array{Float64,4})
    return permutedims(tensor, (1, 3, 2, 4))
end

"""
 Convert the rank-four tensor representing the two-body integrals from chemists' index
 to physicists' index order: i,j,k,l -> i,l,j,k
 Note: when switching from chemist's notation to phycisist's notation, the symmetry Structure
 changes. In particular PERM_2_AB, which swaps the first two indices, p-q, or axis 1 and 2 in
 Julia. This symmetry holds in chemists' notation but does not generally holds in physicists'
 notation.
"""
function chem_to_phys(tensor::Array{Float64,4})
    return permutedims(tensor, (1, 3, 2, 4))  
end

function find_index_order(two_body_tensor::Array{Float64, 4}; rtol=1e-5, atol=1e-8)::IndexType
    if all_permutation_symmetries_hold(two_body_tensor; rtol=rtol, atol=atol)
        return CHEMIST
    end

    t_phys = phys_to_chem(two_body_tensor)
    if all_permutation_symmetries_hold(t_phys; rtol=rtol, atol=atol)
        return PHYSICIST
    end

    t_intermediate = phys_to_chem(t_phys)
    if all_permutation_symmetries_hold(t_intermediate; rtol=rtol, atol=atol)
        return INTERMEDIATE
    end

    return UNKNOWN
end