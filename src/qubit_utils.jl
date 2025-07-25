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
 Returns a dictionary with the combined terms.
"""
function combine_pauli_terms(generators::Vector{Pauli{N}}, parameters::Vector{ComplexF64}) where N
    combined = Dict{Pauli{N}, ComplexF64}()

    for (p, coeff) in zip(generators, parameters)
        if haskey(combined, p)
            combined[p] += coeff
        else
            combined[p] = coeff
        end
    end

    #filter out negligible terms
    filtered = Dict{Pauli{N}, ComplexF64}()
    for (p, c) in combined
        if abs(real(c)) > 1e-8 || abs(imag(c)) > 1e-8
            filtered[p] = c
        end
    end
    return combined #filtered
end

"""
 Returns the combined terms as resulting ordered vectors of generators and coefficients.
"""
function combine_terms(generators::Vector{Pauli{N}}, parameters::Vector{ComplexF64}) where N
    combined = Dict{Pauli{N}, ComplexF64}()
    geners = Vector{Pauli{N}}()
    params = Vector{Float64}()

    for (p, coeff) in zip(generators, parameters)
        if haskey(combined, p)
            combined[p] += coeff
        else
            combined[p] = coeff
        end
    end

    #filter out negligible terms
    for (p, c) in combined
        if abs(real(c)) > 1e-8 || abs(imag(c)) > 1e-8
            push!(geners, p)
            push!(params, c) 
        end
    end
    return geners, params
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

    # One-body terms
    for p in 1:N
        for q in 1:N
            hval = h1[p, q]
            if abs(hval) > 1e-7
                a_dag = JW_creator_mapping(N, p)
                a = JW_annihilator_mapping(N, q)
                one_e_term += hval * (a_dag * a)
            end
        end
    end

    # Two-body terms
    for p in 1:N, q in 1:N, r in 1:N, s in 1:N
        coeff = h2[p,q,r,s] - h2[p,q,s,r] #antisymmetry
        abs(coeff) > 1e-10 || continue

        A = JW_creator_mapping(N, p) *
            JW_creator_mapping(N, q) *
            JW_annihilator_mapping(N, s) *
            JW_annihilator_mapping(N, r)

        A_hc = JW_creator_mapping(N, r) *
               JW_creator_mapping(N, s) *
               JW_annihilator_mapping(N, q) *
               JW_annihilator_mapping(N, p)

        two_e_term += -0.25* coeff * (A + A_hc)
    end
    

    # Add constant term h0 (identity)
    if abs(h0) > 1e-10
        push!(generators, Pauli(N))  # Identity
        push!(parameters, h0)
    end

    # Collect Pauli terms
    for (pauli, coeff) in one_e_term
        abs(coeff) > 1e-10 || continue
        push!(generators, Pauli(pauli))
        push!(parameters, coeff)
    end

    for (pauli, coeff) in two_e_term
        abs(coeff) > 1e-10 || continue
        push!(generators, Pauli(pauli))
        push!(parameters, coeff)
    end

    return generators, parameters
end