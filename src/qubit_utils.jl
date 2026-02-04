"""
 # # # Qubit utilities
"""

"""
  Improved JW mapper
"""
# --- Helper: Jordan-Wigner Mapping ---
function fermion_op(N::Int, i::Int; dagger::Bool=true)
    # Indices in the math are 1-based (1 to N), but bit operations are 0-based.
    bit_idx = i - 1
    one_val = Int128(1)
    two = Int128(2)
    # Z-string mask: 1s from bit 0 to i-2
    z_mask = two^(bit_idx) - one_val
    # Active bit mask: The specific qubit i
    x_mask = two^(bit_idx)
    # Pauli{N}(phase, z_mask, x_mask)
    ax_term = PauliBasis(Pauli{N}(1, z_mask, x_mask))
    # Y acts as both X and Z on the same qubit.
    ay_term = PauliBasis(Pauli{N}(1, z_mask | x_mask, x_mask))
    if dagger
        # a^dagger = 0.5 * (X - iY)
        return 0.5 * (ax_term - im * ay_term)
    else
        # a = 0.5 * (X + iY)
        return 0.5 * (ax_term + im * ay_term)
    end
end

function molecular_hamiltonian_uhf(N_spin_orbitals::Int, 
                                   path_h1_a::String, path_h1_b::String,
                                   path_h2_aa::String, path_h2_bb::String, path_h2_ab::String; 
                                   tol=1e-10)
    
    # 1. Determine Spatial Dimensions
    N_spatial = div(N_spin_orbitals, 2)
    println("Building Hamiltonian from Separate Spin Blocks...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # 2. Load Integrals
    println("  -> Loading Integrals...")
    h1_a = npzread(path_h1_a)
    h1_b = npzread(path_h1_b)
    h2_aa = npzread(path_h2_aa)
    h2_bb = npzread(path_h2_bb)
    h2_ab = npzread(path_h2_ab)

    # 3. Pre-calculate Operators
    println("  -> Generating Operator Cache...")
    # ops_dag[i] is a^dagger_i
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

    H = PauliSum(N_spin_orbitals, ComplexF64)
    
    # --- PART 1: ONE-BODY TERMS ---
    println("  -> Processing 1-Body terms (Alpha & Beta)...")
    # for p in 1:N_spatial, q in 1:N_spatial
    for q in 1:N_spatial, p in 1:N_spatial

        val = h1_a[p, q]
        if abs(val) > tol
            # Alpha indices: 1 to N_spatial
            sum!(H, val * ops_dag[p] * ops_col[q])
        end

        val = h1_b[p, q]
        if abs(val) > tol
            # Beta indices: N_spatial+1 to 2*N_spatial
            P, Q = p + N_spatial, q + N_spatial
            sum!(H, val * ops_dag[P] * ops_col[Q])
        end
    end

    # Free memory for 1-body
    h1_a = nothing
    h1_b = nothing

    # --- PART 2: TWO-BODY TERMS ---
    # Interaction Form: 1/2 * sum V_pqrs * a_p^dag a_r^dag a_s a_q
    println("  -> Processing 2-Body terms (AA, BB, AB)...")

    # for s in 1:N_spatial, r in 1:N_spatial, q in 1:N_spatial, p in 1:N_spatial
    for p in 1:N_spatial, q in 1:N_spatial, r in 1:N_spatial, s in 1:N_spatial

        # --- 1. Alpha-Alpha (aa|aa) ---
        val_aa = h2_aa[p, q, r, s]
        if abs(val_aa) > tol
            sum!(H, 0.5 * val_aa * ops_dag[p] * ops_dag[r] * ops_col[s] * ops_col[q])
        end

        # --- 2. Beta-Beta (bb|bb) ---
        val_bb = h2_bb[p, q, r, s]
        if abs(val_bb) > tol
            # Shift all indices: P, Q, R, S
            P, Q = p + N_spatial, q + N_spatial
            R, S = r + N_spatial, s + N_spatial
            sum!(H, 0.5 * val_bb * ops_dag[P] * ops_dag[R] * ops_col[S] * ops_col[Q])
        end

        # --- 3. Alpha-Beta (ab|ab) ---
        val_ab = h2_ab[p, q, r, s]
        if abs(val_ab) > tol
            # Shift only Beta indices (r -> R, s -> S)
            R, S = r + N_spatial, s + N_spatial
            # NO 0.5 factor (distinguishable sets)
            sum!(H, val_ab * ops_dag[p] * ops_dag[R] * ops_col[S] * ops_col[q])
        end
    end

    # Clear memory after the massive loop finishes
    h2_aa = nothing
    h2_bb = nothing
    h2_ab = nothing


    println("  -> Hamiltonian Construction Complete.")
    coeff_clip!(H) 
    return H
end

function molecular_hamiltonian(N_spin_orbitals::Int, path::String; tol=1e-10)
    # 1. Determine Spatial Dimensions
    N_spatial = div(N_spin_orbitals, 2)
    println("Building RHF Hamiltonian Directly (Memory Efficient)...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # 2. Load Spatial Integrals (Small arrays)
    # h1 is N x N
    # h2 is N x N x N x N
    H_op = npzread(path)
    H0 = H_op["hc"][1]
    h1_spatial = H_op["h1e"]
    h2_spatial = H_op["h2e"]

    # 3. Pre-calculate Operators
    # We store these to avoid re-allocating Pauli objects constantly.
    println("  -> Generating Operator Cache...")
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

    H = PauliSum(N_spin_orbitals, ComplexF64)
    
    if isnan(H0)
        println("WARNING: H0 term is NaN, skiping term...")
        elseif iszero(H0)
            println("H0 term is 0.0, skipping term...")
        else
            H += H0 * Pauli(N_spin_orbitals) 
    end 

    # --- PART 1: ONE-BODY TERMS ---
    # Hamiltonian: sum h_pq ( a^dag_p_alpha a_q_alpha + a^dag_p_beta a_q_beta )
    println("  -> Processing 1-Body terms...")
    
    for p in 1:N_spatial, q in 1:N_spatial
        val = h1_spatial[p, q]
        
        if abs(val) > tol
            # Alpha Term (indices 1...N)
            sum!(H, val * ops_dag[p] * ops_col[q])
            
            # Beta Term  (indices N+1...2N)
            # Shift indices by N_spatial
            P, Q = p + N_spatial, q + N_spatial
            sum!(H, val * ops_dag[P] * ops_col[Q])
        end
    end
    # Clear h1 to free memory (optional, but good practice)
    h1_spatial = nothing 

    # --- PART 2: TWO-BODY TERMS ---
    # Hamiltonian: 1/2 * sum (pq|rs) * [Terms]
    # We iterate spatial indices once and add contributions to all 4 spin sectors.
    println("  -> Processing 2-Body terms...")
    
    count = 0
    for p in 1:N_spatial, q in 1:N_spatial, r in 1:N_spatial, s in 1:N_spatial
        val = h2_spatial[p, q, r, s]
        temp = PauliSum(N_spin_orbitals, ComplexF64)
        # OPTIMIZATION: Check tolerance immediately
        if abs(val) > tol
            # Pre-calculate shifted indices for Beta
            P, Q = p + N_spatial, q + N_spatial
            R, S = r + N_spatial, s + N_spatial

            # 1. Alpha-Alpha Sector (p, q, r, s all Alpha)
            # Term: a_p^dag a_r^dag a_s a_q
            # term_aa = ops_dag[p] * ops_dag[r] * ops_col[s] * ops_col[q]
            # sum!(H, 0.5 * val * term_aa)
            sum!(temp, 0.5 * val * ops_dag[p] * ops_dag[r] * ops_col[s] * ops_col[q])
            


            # 2. Beta-Beta Sector (P, Q, R, S all Beta)
            # term_bb = ops_dag[P] * ops_dag[R] * ops_col[S] * ops_col[Q]
            # sum!(H, 0.5 * val * term_bb)
            sum!(temp, 0.5 * val * ops_dag[P] * ops_dag[R] * ops_col[S] * ops_col[Q])



            # 3. Mixed Spin Sectors
            # Alpha-Beta: p,q are Alpha; R,S are Beta
            # Term: a_p^dag a_R^dag a_S a_q
            # term_ab = ops_dag[p] * ops_dag[R] * ops_col[S] * ops_col[q]
            # sum!(H, 0.5 * val * term_ab)
            sum!(temp, 0.5 * val * ops_dag[p] * ops_dag[R] * ops_col[S] * ops_col[q])



            # Beta-Alpha: P,Q are Beta; r,s are Alpha
            # Term: a_P^dag a_r^dag a_s a_Q
            # term_ba = ops_dag[P] * ops_dag[r] * ops_col[s] * ops_col[Q]
            # sum!(H, 0.5 * val * term_ba)
            sum!(temp, 0.5 * val * ops_dag[P] * ops_dag[r] * ops_col[s] * ops_col[Q])

            count += 1
            sum!(H, temp)
        end
    end
    
    # # Force Garbage Collection
    # h2_spatial = nothing 
    # GC.gc()

    println("     (Added $count significant spatial interaction terms)")
    coeff_clip!(H)
    return H
end


#####################################################
"""
 OLD JW mapper
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

function get_qubit_hamiltonian(scf_obj)
    # Atomic integrals
    ao_hcore = scf_obj.T + scf_obj.V
    ao_eris = scf_obj.I

    # Transform Atomic integrals into Molecular integrals
    C = scf_obj.C # MO coeffs
    mo_hcore , mo_eris = QuantumChemQC.ao2mo_coefficients(C, ao_hcore, ao_eris)

    # Spinorbital tensors
    core =0 
    h0 = scf_obj.Enuc + core 
    h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)
    N = size(h1, 1)
    qubit_paulis, qubit_coeffs = QuantumChemQC.qubit_hamiltonian(N, h0, h1, h2) 

    gen, coeff = QuantumChemQC.combine_terms(qubit_paulis, qubit_coeffs)

    H = PauliSum(N)
    H = foldl(+, (c*p for (p,c) in zip(gen, coeff)); init=PauliSum(N))
    coeff_clip!(H)

return H, N;
end

function PauliSum_hamiltonian(N::Int64, h0::Float64, h1::Matrix{Float64}, h2::Array{Float64, 4})
    H = PauliSum(N, Float64)
    #one_e_term = PauliSum(N, Float64)
    #two_e_term = PauliSum(N, Float64)
    
    # H0 term (identity)
    H += h0 * Pauli(N)
    
    # One-body terms
    for p in 1:N
        for q in 1:N
            hval = h1[p, q]
            if abs(hval) > 1e-7
                a_dag = JW_creator_mapping(N, p)
                a = JW_annihilator_mapping(N, q)
                #one_e_term += hval * (a_dag * a)
                H += hval * (a_dag * a)
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

        #two_e_term += -0.25* coeff * (A + A_hc)
        H += -0.25* coeff * (A + A_hc)
    end

    coeff_clip!(H)
    return H
end