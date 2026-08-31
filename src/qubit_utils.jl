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
    PauliOperators.coeff_clip!(H)

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
    PauliOperators.coeff_clip!(H, coeff_thresh_clip)
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