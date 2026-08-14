"""
 # # # Qubit utilities
"""
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
    in the JW representation, as input it uses the resulting spinorbital tensors, 
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

"""
  Creates Molecular Hamiltoinan (PauliSum) using spinorbital tensors or Molecular integrals in the physics notation.
  N :: Total number of spinorbitals.
  h0 :: Repulsion term
  h1 :: One particle tensor
  h2 :: Two particle tensor
  MoInts :: bool, controls if the transformation is made with spinorbitals (False) or with Molecular integrals (true)
  NOI :: bool, controls if include the all-I term (nuclear generally)
  Returns a PauliSum Hamiltonian with interleaved ordering.
"""
function PauliSum_hamiltonian(N::Int64, h0::Float64, h1::Matrix{Float64}, h2::Array{Float64, 4}; NOI=true, MoInts=false)
    H = PauliSum(N, Float64)
    #one_e_term = PauliSum(N, Float64)
    #two_e_term = PauliSum(N, Float64)
    
    # H0 term (identity)
    if !NOI
        H += h0 * Pauli(N)
    end

    if MoInts
        N = 2*N #N_spin_orbitals
        h1, h2 = QuantumChemQC.get_spin_orbital_tensors(h1, h2)
    end

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

function molecular_hamiltonian_uhf(N_spatial::Int, path::String; tol=1e-10, NOI=true, block=true)
    coeff_thresh_clip = 1e-10
    N_spin_orbitals = 2 * N_spatial
    println("Building UHF Hamiltonian (memory-efficient)...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # Load integrals container (assume npzread/path returns a Dict-like object)
    H_op = npzread(path)

    println("Building Hamiltonian from Separate Spin Blocks...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # Load integrals (spatial)
    println("  -> Loading Integrals...")
    h1_a = H_op["h1_a"]
    h1_b = H_op["h1_b"]
    h2_aa = H_op["h2aa"]
    h2_bb = H_op["h2bb"]
    h2_ab = H_op["h2ab"]


    H = PauliSum(N_spin_orbitals, ComplexF64)

    # optional constant term
    if NOI
        H0 = 0
    else
        H0 = H_op["hc"][1]
    end

    if isnan(H0)
        println("WARNING: H0 term is NaN, skiping term...")
        elseif iszero(H0)
            println("H0 term is 0.0, skipping term...")
        else
            H += H0 * Pauli(N_spin_orbitals) 
    end 

    # Index mapping: block vs interleaved
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end

    # Pre-calculate Fermionic operators (per spin-orbital)
    println("  -> Generating Operator Cache...")
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

    # PART 1: ONE-BODY TERMS
    println("  -> Processing 1-Body terms (Alpha & Beta)...")
    for p in 1:N_spatial, q in 1:N_spatial
    #for q in 1:N_spatial, p in 1:N_spatial
        val_a = h1_a[p, q]
        if abs(val_a) > tol
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(H, val_a * ops_dag[ia] * ops_col[ja])
        end

        val_b = h1_b[p, q]
        if abs(val_b) > tol
            ib = idx_beta(p)
            jb = idx_beta(q)
            sum!(H, val_b * ops_dag[ib] * ops_col[jb])
        end
    end
    coeff_clip!(H)

    # free 1-body integrals
    h1_a = nothing
    h1_b = nothing

    # PART 2: TWO-BODY TERMS
    # Interaction form (spatial integrals):
    #   aa: 1/2 * sum_{pqrs} h2_aa[p,q,r,s] a_p^† a_r^† a_s a_q  (and similarly for bb)
    #   ab: sum_{pqrs} h2_ab[p,q,r,s] a_p^† a_R^† a_S a_q  (no 1/2 for distinguishable spin sectors)
    println("  -> Processing 2-Body terms (AA, BB, AB)...")

    for p in 1:N_spatial, q in 1:N_spatial, r in 1:N_spatial, s in 1:N_spatial
        #temp = PauliSum(N_spin_orbitals, ComplexF64)

        # --- Alpha-Alpha (aa|aa) ---
        val_aa = h2_aa[p, q, r, s]
        if abs(val_aa) > tol
            ip = idx_alpha(p); iq = idx_alpha(q); ir = idx_alpha(r); is_ = idx_alpha(s)
            # 1/2 factor for same-spin two-body terms
            sum!(H, 0.5 * val_aa * ops_dag[ip] * ops_dag[ir] * ops_col[is_] * ops_col[iq])
        end

        # --- Beta-Beta (bb|bb) ---
        val_bb = h2_bb[p, q, r, s]
        if abs(val_bb) > tol
            ip = idx_beta(p); iq = idx_beta(q); ir = idx_beta(r); is_ = idx_beta(s)
            # 1/2 factor for same-spin two-body terms
            sum!(H, 0.5 * val_bb * ops_dag[ip] * ops_dag[ir] * ops_col[is_] * ops_col[iq])
        end

        # --- Alpha-Beta (ab|ab) ---
        # Alpha-Beta: p,q are Alpha; R,S are Beta
        # Term: a_p^dag a_R^dag a_S a_q
        # here p,q are alpha indices, r,s are beta indices (distinguishable => no 1/2)
        val_ab = h2_ab[p, q, r, s]
        if abs(val_ab) > tol
            ia_p = idx_alpha(p); ia_q = idx_alpha(q)
            ib_r = idx_beta(r);  ib_s = idx_beta(s)
            # create alpha p, annihilate alpha q; create beta r, annihilate beta s
            sum!(H, val_ab * ops_dag[ia_p] * ops_dag[ib_r] * ops_col[ib_s] * ops_col[ia_q])
        end



        # accumulate
    end

    # free two-body integrals
    h2_aa = nothing
    h2_bb = nothing
    h2_ab = nothing

    println("  -> Hamiltonian Construction Complete.")
    coeff_clip!(H, thresh=coeff_thresh_clip)
    return H
end


"""
 Fast molecular Hamiltonian (PauliSum) computed from Molecular Integrals in chemist notation.
 N_spin_orbitals :: Number of spinorbitals is always the double of the number of MOs (N_spatial).
 Path :: Path to the npz file which stores the integrals.
 NOI :: Allow to disable IIII.. constant term (treated as zero when NOI is true)
 block :: Structure of the qubit hamiltonian by default, if block=false, it creates an interleaved hamiltonian.
"""
function molecular_hamiltonian(N_spatial::Int, path::String; tol=1e-10, NOI=true, block=true)
    coeff_thresh_clip = 1e-6
    N_spin_orbitals = 2*N_spatial
    # 1. Determine Spatial Dimensions
    println("Building RHF Hamiltonian Directly (Memory Efficient)...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # 2. Load Spatial Integrals (Small arrays)
    # h1 is N x N
    # h2 is N x N x N x N
    H_op = npzread(path)

    if NOI #Accept III... term (not necessary for AS calculations)
        H0 = 0 #So this will be skipped
    else
        H0 = H_op["hc"][1] 
    end

    h1_spatial = H_op["h1e"]
    h2_spatial = H_op["h2e"]

    # 3. Pre-calculate Operators
    # We store these to avoid re-allocating Pauli objects constantly.
    println("  -> Generating Operator Cache...")
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

     # --- INDEX MAPPING: block vs interleaved ---
    # block=true:  alpha indices  = 1:N_spatial, beta indices = N_spatial+1:2N_spatial
    # block=false: interleaved  = [alpha,beta,alpha,beta,...]  => alpha = 2*p-1, beta = 2*p
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end

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
            # Alpha Term
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(H, val * ops_dag[ia] * ops_col[ja])
            
            # Beta Term
            ib = idx_beta(p)
            jb = idx_beta(q)
            sum!(H, val * ops_dag[ib] * ops_col[jb])
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
            # map indices according to chosen ordering
            ap = idx_alpha(p); aq = idx_alpha(q); ar = idx_alpha(r); as_ = idx_alpha(s)
            bp = idx_beta(p);  bq = idx_beta(q);  br = idx_beta(r);  bs = idx_beta(s)

            # 1. Alpha-Alpha Sector (p, q, r, s all Alpha)
            # Term: a_p^dag a_r^dag a_s a_q
            sum!(temp, 0.5 * val * ops_dag[ap] * ops_dag[ar] * ops_col[as_] * ops_col[aq])         

            # 2. Beta-Beta Sector (P, Q, R, S all Beta)
            sum!(temp, 0.5 * val * ops_dag[bp] * ops_dag[br] * ops_col[bs] * ops_col[bq])

            # 3. Mixed Spin Sectors
            # Alpha-Beta: p,q are Alpha; R,S are Beta
            sum!(temp, 0.5 * val * ops_dag[ap] * ops_dag[br] * ops_col[bs] * ops_col[aq])

            # Beta-Alpha: P,Q are Beta; r,s are Alpha
            sum!(temp, 0.5 * val * ops_dag[bp] * ops_dag[ar] * ops_col[as_] * ops_col[bq])

            count += 1
            sum!(H, temp)
        end
    end
    
    # # Force Garbage Collection
    # h2_spatial = nothing 
    # GC.gc()

    println("     (Added $count significant spatial interaction terms)")
    PauliOperators.coeff_clip!(H, coeff_thresh_clip)
    return H
end

"""
   Dipole moment in JW representation
   This function uses integrals from a Restricted-SCF type calculation
"""
function R_dipole_moment_op(N_spatial::Int, path::String; tol=1e-10, NOI=true, block=true)
    coeff_thresh_clip = 1e-8
    N_spin_orbitals = 2*N_spatial
    # 1. Determine Spatial Dimensions
    println("Building Dipole Operator Directly (Memory Efficient)...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # 2. Load Spatial Integrals (Small arrays)
    # h1 is N x N
    # h2 is N x N x N x N
    Op = npzread(path)
    dip_op = Op["dip_op"]

    # 3. Pre-calculate Operators
    # We store these to avoid re-allocating Pauli objects constantly.
    println("  -> Generating Operator Cache...")
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

     # --- INDEX MAPPING: block vs interleaved ---
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end

    D = PauliSum(N_spin_orbitals, ComplexF64)

    # --- PART 1: ONE-BODY TERMS ---
    # Hamiltonian: sum h_pq ( a^dag_p_alpha a_q_alpha + a^dag_p_beta a_q_beta )
    println("  -> Processing 1-Body terms...")
    dipX = dip_op[1, :, :]
    dipY = dip_op[2, :, :]
    dipZ = dip_op[3, :, :]

    #- - - r_x terms - - -
    for p in 1:N_spatial, q in 1:N_spatial
        val = dipX[p, q]
        if abs(val) > tol
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(D, val * ops_dag[ia] * ops_col[ja])
        end
    end
    coeff_clip!(D, thresh=coeff_thresh_clip)

    #- - - r_y terms - - -
    for p in 1:N_spatial, q in 1:N_spatial
        val = dipY[p, q]
        if abs(val) > tol
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(D, val * ops_dag[ia] * ops_col[ja])
        end
    end
    coeff_clip!(D)

    #- - - r_z terms - - -
    for p in 1:N_spatial, q in 1:N_spatial
        val = dipZ[p, q]
        if abs(val) > tol
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(D, val * ops_dag[ia] * ops_col[ja])
        end
    end

    # Clear h1 to free memory (optional, but good practice)
    dip_op = nothing
    coeff_clip!(D, thresh=coeff_thresh_clip)

    return D
end

"""
   Dipole moment compoentn (x, y or z) JW representation
   This function uses integrals from a Restricted-SCF type calculation
    the path should point to a npz file which contains the dipole moment integrals in the format 
    (N, N) where the first dimension corresponds to the x, y or z components of the dipole moment.
"""
function xyz_dipole_moment_op(N_spatial::Int, path::String; tol=1e-10, NOI=true, block=true)
    coeff_thresh_clip = 1e-8
    N_spin_orbitals = 2*N_spatial
    # 1. Determine Spatial Dimensions
    println("Building Dipole Operator Directly (Memory Efficient)...")
    println("  -> System: $N_spatial Spatial Orbitals ($N_spin_orbitals Spin Orbitals)")

    # 2. Load Spatial Integrals (Small arrays)
    # h1 is N x N
    # h2 is N x N x N x N
    Op = npzread(path)
    dip_op = Op["dip_op"]

    # 3. Pre-calculate Operators
    # We store these to avoid re-allocating Pauli objects constantly.
    println("  -> Generating Operator Cache...")
    ops_dag = [fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

     # --- INDEX MAPPING: block vs interleaved ---
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end

    D = PauliSum(N_spin_orbitals, ComplexF64)

    # --- PART 1: ONE-BODY TERMS ---
    # Hamiltonian: sum h_pq ( a^dag_p_alpha a_q_alpha + a^dag_p_beta a_q_beta )
    println("  -> Processing 1-Body terms...")

    #- - - r_x terms - - -
    for p in 1:N_spatial, q in 1:N_spatial
        val = dip_op[p, q]
        if abs(val) > tol
            ia = idx_alpha(p)
            ja = idx_alpha(q)
            sum!(D, val * ops_dag[ia] * ops_col[ja])
        end
    end
    coeff_clip!(D, thresh=coeff_thresh_clip)

    # Clear h1 to free memory (optional, but good practice)
    dip_op = nothing
    coeff_clip!(D, thresh=coeff_thresh_clip)

    return D
end

"""
 Function that builds the HOMO-LUMO excitation operator based on a given input ket.
 This function can retrieve a spin_conserving operator, based on the block or interlieved 
 encoding used in the Hamiltonian construction. CURRENTLY WORKS FOR SINGLET REF. STATES
 * * *
 NOTE: Function in current development, based only in the alpha-channel excitation (spin conserving),
       but it can contain the beta-channel excitation via uncommenting the lines below.
 * * *
"""
function homo_lumo_excitation_op(ket::Ket{N}; spin_conserving=true, block=true) where N
    coeff_thresh_clip = 1.0e-6
    N_spatial = Int(N ÷ 2)

    # --- INDEX MAPPING: block vs interleaved ---
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end
    
    # --- Operator cache ---
    ops_dag = [QuantumChemQC.fermion_op(N, i, dagger=true) for i in 1:N]
    ops_col = [QuantumChemQC.fermion_op(N, i, dagger=false) for i in 1:N]

    #1) Read Ket occupations
    v = ket.v
    occ = digits(v, base=2, pad=N)  

    #2) Determine HOMO and LUMO indices
    function find_homo_lumo(idx_map)
        occ_indices = Int[]
        virt_indices = Int[]

        for p in 1:N_spatial
            i = idx_map(p)
            if occ[i] == 1
                push!(occ_indices, p)
            else
                push!(virt_indices, p)
            end
        end
        
        isempty(occ_indices) && return nothing
        isempty(virt_indices) && return nothing

        homo = maximum(occ_indices)
        lumo = minimum(virt_indices)

        return homo, lumo
    end

    #3) Build HOMO-LUMO excitation operator: O + h.c.
    O = PauliSum(N, ComplexF64)
    #Establish possible excitation based on spin_cons
    if spin_conserving
        # α channel
        hl_alpha = find_homo_lumo(idx_alpha)
        if hl_alpha !== nothing
            h, l = hl_alpha
            sum!(O, ops_dag[idx_alpha(l)] * ops_col[idx_alpha(h)])
            sum!(O, ops_dag[idx_alpha(h)] * ops_col[idx_alpha(l)]) #h.c
        end

        # β channel
        #hl_beta = find_homo_lumo(idx_beta)
        #if hl_beta !== nothing
        #    h, l = hl_beta
        #    sum!(O, ops_dag[idx_beta(l)] * ops_col[idx_beta(h)])
        #    sum!(O, ops_dag[idx_beta(h)] * ops_col[idx_beta(l)]) #h.c
        #end

    else
        # Global HOMO/LUMO (ignore spin)
        occ_indices = findall(x -> x ==1, occ)
        virt_indices = findall(x -> x ==0, occ)
        println("OCC indx ", occ_indices)
        println("VIRT indx ", virt_indices)

        if !isempty(occ_indices) && !isempty(virt_indices)
            h = maximum(occ_indices)
            l = minimum(virt_indices)

            println("HOMO idx ", h)
            println("LUMO idx ", l)

            sum!(O, ops_dag[l] * ops_col[h])
            sum!(O, ops_dag[h] * ops_col[l]) #h.c
        end
    end

    coeff_clip!(O, thresh=coeff_thresh_clip)
    return O
end