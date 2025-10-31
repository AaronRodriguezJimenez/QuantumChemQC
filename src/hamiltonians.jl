# - - - - - - - - - - - - -
# Spin Hamiltonians
# - - - - - - - - - - - - -

function heisenberg_1D(N, Jx, Jy, Jz; x=0, y=0, z=0)
    H = PauliSum(N, Float64)
    for i in 0:N-1
        H += -2*Jx * Pauli(N, X=[i+1,(i+1)%(N)+1])
        H += -2*Jy * Pauli(N, Y=[i+1,(i+1)%(N)+1])
        H += -2*Jz * Pauli(N, Z=[i+1,(i+1)%(N)+1])
    end 
    for i in 1:N
        H += x * Pauli(N, X=[i])
        H += y * Pauli(N, Y=[i])
        H += z * Pauli(N, Z=[i])
    end
    coeff_clip!(H) 
    return H
end

# - - - - - - - - - - - - -
# Fermionic Hamiltonians
# - - - - - - - - - - - - -
"""
   HELPERS
 The following functions are designed to perform the Jordan-Wigner mapping.
"""
# --- helpers for JWmapping ---
@inline ubit(i::Int) = UInt128(1) << (i-1)                         # bit at site i (1-based)
@inline umask_lt(i::Int) = i==1 ? UInt128(0) : (ubit(i) - UInt128(1))  # bits < i
@inline umask_le(i::Int) = umask_lt(i) | ubit(i)                      # bits ≤ i

# --- JW mapping (original real ± form; your coeff() supplies +i on ZX sites) ---
function JWmapping(N; i::Int, j::Int)
    1 <= i <= N || throw(DimensionMismatch("site i=$i out of 1:$N"))
    1 <= j <= N || throw(DimensionMismatch("site j=$j out of 1:$N"))

    # X pieces with Z-strings
    ax = Pauli{N}(1, reinterpret(Int128, umask_lt(i)), reinterpret(Int128, ubit(i)))  # Z^{<i} X_i
    bx = Pauli{N}(1, reinterpret(Int128, umask_lt(j)), reinterpret(Int128, ubit(j)))  # Z^{<j} X_j

    # "Y" pieces = Z^{≤i} X_i, Z^{≤j} X_j  (no explicit im; your coeff() turns ZX into iY)
    ay = Pauli{N}(1, reinterpret(Int128, umask_le(i)), reinterpret(Int128, ubit(i)))
    by = Pauli{N}(1, reinterpret(Int128, umask_le(j)), reinterpret(Int128, ubit(j)))

    # c†_i = (X_i - Y_i)/2,  c_j = (X_j + Y_j)/2   in your convention
    c_dagg_i = 0.5 * (ax - ay)
    c_j      = 0.5 * (bx + by)

    return c_dagg_i * c_j
end

"""
    fermi_hubbard_2D(Lx, Ly, t, U; reverse_ordering=false)

Construct generators and parameters for the 2D spinful Hubbard model on Lx×Ly
(physical sites). Each physical site has two spin-orbitals (up, down), so
total qubits N must equal 2 * Lx * Ly.

Returns (generators::Vector{Pauli{N}}, parameters::Vector{Float64}).
"""
function fermi_hubbard_2D(Lx::Int, Ly::Int, t::Float64, U::Float64)
    Nsites = Lx * Ly
    N_total = 2 * Nsites   # Total number of fermionic modes (spin up and down)
    H = PauliSum(N_total, Float64)

    if 2 * Nsites != N_total
        throw(ArgumentError("Total qubits N must equal 2 * Lx * Ly. Got N=$N_total, Lx*Ly=$Nsites"))
    end

    up(j) = 2*j - 1
    dn(j) = 2*j
    linear_index(x,y) = (x - 1) * Ly + y   # x in 1:Lx, y in 1:Ly

    # HOPPING: loop nearest-neighbour pairs once, add c_i^† c_j + c_j^† c_i (both spins)
    for x in 1:Lx, y in 1:Ly
        jsite = linear_index(x, y)
         # neighbor +x (right in x)
        if x < Lx
            isite = linear_index(x + 1, y)
            for spin in (up, dn)
                m = spin(jsite)   # mode index for j
                n = spin(isite)   # mode index for i
                term = JWmapping(N_total, i=m, j=n) + JWmapping(N_total, i=n, j=m)
                H += -t * term
            end
        end
         # neighbor +y (right in y)
         if y < Ly
            isite = linear_index(x, y + 1)
            for spin in (up, dn)
                m = spin(jsite)
                n = spin(isite)
                term = JWmapping(N_total, i=m, j=n) + JWmapping(N_total, i=n, j=m)
                H += -t * term
            end
        end
    end

    for i in 1:Nsites
        a_up = 2*i - 1   # spin-up orbital index
        a_dn = 2*i       # spin-down orbital index
        interaction_term = U *JWmapping(N_total, i=a_up, j=a_up) * JWmapping(N_total, i=a_dn, j=a_dn)

        H += interaction_term
    end

    # Filter zero coefficients
    coeff_clip!(H)

    return H
end


"""
 - - - Fermi-Hubbard model 2D (snake/zizag ordering) - - -
 The following function constructs the Hamiltonian for the 2D Fermi-Hubbard model
    using a snake-like (zigzag) ordering of the lattice sites.
  Constucts the model on a Lx x Ly lattive of physical sites,
  each site has two spin-orbitals (up, down), so the total number of fermionic modes is 
  N_total = 2 * Lx * Ly.
  If 'snake_ordering' is true, the physical sites follow a row-wise "snake/zigzag" pattern.
  Returns a PauliSum representing the Hamiltonian.
"""
function fermi_hubbard_2D_snake(Lx::Int, Ly::Int, t::Float64, U::Float64; snake_ordering::Bool=false)
    Nsites  = Lx * Ly
    N_total = 2 * Nsites
    H = PauliSum(N_total, Float64)

    # Spin-orbital indices for site j (1-based site indexing)
    up(j) = 2*j - 1
    dn(j) = 2*j

    # Map (row x, col y) -> site index j in [1, Nsites]
    @inline function site_index(x::Int, y::Int)
        if !snake_ordering
            # Row-major
            return (x - 1) * Ly + y
        else
            # Snake (zigzag) row-major:
            # odd rows (x=1,3,...) go left->right
            # even rows (x=2,4,...) go right->left
            if isodd(x)
                return (x - 1) * Ly + y
            else
                return x * Ly - (y - 1)
            end
        end
    end

    # HOPPING: nearest-neighbor pairs (right and down) — add h.c.; both spins
    for x in 1:Lx, y in 1:Ly
        jsite = site_index(x, y)

        # neighbor in +x (next row)
        if x < Lx
            isite = site_index(x + 1, y)
            for spin in (up, dn)
                m = spin(jsite)
                n = spin(isite)
                term = JWmapping(N_total, i=m, j=n) + JWmapping(N_total, i=n, j=m)
                H += -t * term
            end
        end

        # neighbor in +y (next column)
        if y < Ly
            isite = site_index(x, y + 1)
            for spin in (up, dn)
                m = spin(jsite)
                n = spin(isite)
                term = JWmapping(N_total, i=m, j=n) + JWmapping(N_total, i=n, j=m)
                H += -t * term
            end
        end
    end

    # On-site interaction U * n_up * n_dn
    for i in 1:Nsites
        a_up = 2*i - 1
        a_dn = 2*i
        interaction_term = U * JWmapping(N_total, i=a_up, j=a_up) *
                               JWmapping(N_total, i=a_dn, j=a_dn)
        H += interaction_term
    end

    coeff_clip!(H)
    return H
end

"""
  Function returnign the Hamiltonian for H2 at a given interatomic distance
    d -> The interatomic distance
"""
function H2_hamiltonian(d)
    mol = """
    0 1
    H  0.0000   0.000  0.000
    H  0.0000   0.000  $d
    """
  
    bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = false)
    scf_obj = QuantumChemQC.SCF(mol, bset, p)

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