using PyCall
using PauliOperators
using DBF
using Printf
using QuantumChemQC
using NPZ

"""
 Here we read the precomputed Active Space spinorbital tensors for N2 from pyqctools,
    and use PauliOperators and DBF to compute its ground state energy.

    Relevant information of the calculation
    CASSCF energy = -107.581589622102
    CASCI E = -107.581589622102  E(CI) = -27.9562154436841  S^2 = 0.0000000
    DBF must converge to E(CI) = -27.9562154436841
    Active space: 10 electrons in 8 orbitals (16 spin orbitals)
"""
function run()
    println("PyCall.python = ", PyCall.python)   # should be /Users/admin/VSCProjects/py4julia/bin/python
    importlib = pyimport("importlib")
    importlib.reload(pyimport("pyqctools"))

    # import the module (lowercase)
    pq = pyimport("pyqctools")
    println("pyqctools imported:", pq)
    Hfcns = pyimport("pyqctools.ham_fcns")


    # Get geometry for H2 molecule and define parameters
    geom = pq.geometries.N2(3.0)    

    println("Molecule geometry from pyqctools:", geom)
    
    # GET MOLECULAR ORBITALS
    mol = pyimport("pyscf.gto").Mole()
    mol.atom = geom
    mol.basis = "sto-3g"
    mol.spin = 0
    mol.build()    

    # RUN SCF
    mf = pyimport("pyscf.scf").UHF(mol).run()
    E_HF = mf.e_tot
    println("Hartree-Fock Energy:", E_HF)
    
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/DBF//N2_tensors_ccsd/3.0_tensors.npz" #uccsd Natural orbs

    data = npzread(data_path)
    H0 = data["hc"][1]
    H1 = data["h1_a"]
    H2 = data["h2aa"]

    println("H0 shape: ", size(H0))
    println(H0)
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))

    Norbs = size(H1,1)  # number of spatial orbitals
    @time H = QuantumChemQC.molecular_hamiltonian_uhf(Norbs, data_path, NOI=false, block=false)

    #QuantumChemQC.coeff_clip!(H, thresh=1e-6)
    return H
end

#- - - Helpers for getting the occupation from bitstring - - - -
# return a PauliOperators Ket equivalent to a given bitstring
function string_to_ket(bits::String)
    b = collect(bits)
    v = parse.(Int128, b)
    N = length(v)
    out = 0
    count = 0

    for bit in v
        if bit%2 == 1
            out += 2^count
        end
        count +=1
    end
    ket = Ket{N}(out)
    return ket, v, out
end

# return 0 or 1 for bit i (1-based indexing, i in 1:N)
get_bit(k::Union{Ket{N}, Bra{N}}, i::Integer) where N = Int((k.v >> (i-1)) & Int128(1))

# return a Vector{Int} of indices where the bit is 1 (1-based)
function get_on_bits(v::Int128)
    out = Int[]
    w = v
    pos = 1
    while w != 0
        if (w & Int128(1)) != 0
            push!(out, pos)
        end
        w >>= 1
        pos += 1
    end
    return out
end

# occupation vector of length N (Int elements 0/1)
function occvec(k::Union{Ket{N}, Bra{N}}) where N
    [ get_bit(k, i) for i in 1:N ]
end

# Now use DBF to compute ground state energy
function dbf_gstate(H::PauliSum{N, T}) where {N,T}

    # = = = DBF parameters = = = 
    max_iter=20
    conv_thresh=1e-10
    evolve_coeff_thresh=1e-10
    grad_coeff_thresh=1e-10
    energy_lowering_thresh=1e-10
    max_rots_per_grad=1
    ket, occ, kidx  = string_to_ket("11111111111111000000") #Leading CAS/sto3g configuration
    #ket, occ, kidx  = string_to_ket("1111100011111000") #Leading CAS/sto3g configuration
    #ket, occ, kidx  = string_to_ket("11111110001111111000") #Leading CAS/sto3g configuration
    #ket, occ, kidx  = string_to_ket("11111111111111000000") #Leading CAS/sto3g configuration SINGLET
    #ket, occ, kidx  = string_to_ket("11111111101010101010") #Leading CAS/sto3g configuration SINGLET
    #ket, occ, kidx  = string_to_ket("1111111111000000") #Leading CAS/sto3g configuration SINGLET
    #ket, occ, kidx  = string_to_ket("1111111110100000") #Leading CAS/sto3g configuration TRIPLET
    #ket, occ, kidx  = string_to_ket("1111111010101000") #Leading CAS/sto3g configuration QUINTET
    #ket, occ, kidx  = string_to_ket("11111111101010101010") #Leading CAS/sto3g configuration SEPTET
    
    println(occ)
    println(N)

    println("Initial state:")
    display(ket)
    #e1 = expectation_value(H,ket)
    #@printf(" Reference = %12.8f\n", e1)

    # Transform H to make |000> the most stable bitstring
    for i in 1:N
        if occ[i] == 1
            H = Pauli(N, X=[i]) * H * Pauli(N, X=[i])
        end
    end

    H0 = deepcopy(H)
    #Hmat = Matrix(H)
    #evals = eigvals(Hmat)
    #@show minimum(evals)

    ψ = Ket([0 for i in 1:N])
    display(ψ)

    e1 = expectation_value(H,ψ)
    @printf(" Reference = %12.8f\n", e1)

    g = Vector{PauliBasis{N}}([]) 
    θ = Vector{Float64}([]) 

    println("\n ########################")
    res = DBF.dbf_groundstate(H0, ψ,
                    max_iter= max_iter,
                    conv_thresh= conv_thresh, 
                    evolve_coeff_thresh= evolve_coeff_thresh,
                    grad_coeff_thresh= grad_coeff_thresh,
                    energy_lowering_thresh= energy_lowering_thresh,
                    max_rots_per_grad = max_rots_per_grad,
                    n_body=1
                    )                        
    
    H = res["hamiltonian"]
    #gi = res["generators"]
    #θi = res["angles"]
        

    e2 = real(expectation_value(H,ψ))

    #- - - Get Calculation Details
    println("\n===============================")
    println("= = = DBF Parameters = = =")
    println("max_iter : ", max_iter) 
    println("conv_thresh : ",conv_thresh)
    println("evolve_coeff_thresh :", evolve_coeff_thresh)
    println("grad_coeff_thresh : ", grad_coeff_thresh)
    println("energy_lowering_thresh : ", energy_lowering_thresh)

    println("Number of qubits: ", N)
    println("Input Hamiltonian #terms: ", length(H))

    # Initial energy
    println("\n===============================")
    println("\n Initial energy:")
    println("Index of min energy bitstring: ", kidx)
    println("Initial state:")
    display(ket)
    println("Occuation vector: ", occ)
    @printf("<H> = %12.8f <U'HU> = %12.8f \n", e1, e2)

    #println("* * * * * * * * * * * * * * * * * ")
    #println("Evolved Hamiltonian")
    #display(H)
                    
    return
end


H = run()
#println("* * * * * * * * * * * * * * * * * ")
#println("Initial Hamiltonian 😆")
#display(H)
#ket, occ, kidx  = string_to_ket("11111110001111111000") #Block
ket, occ, kidx  = string_to_ket("11111111111111000000") #Leading CAS/sto3g configuration SINGLET
println("Initial state:")
display(ket)
e1 = expectation_value(H,ket)
@printf(" Reference = %12.8f\n", e1)

dbf_gstate(H)

#=
# - - - Testing HF energy - - - 

# assume h0, h1_alpha, h1_beta, h2_aa, h2_bb, h2_ab exist in scope
# and have shapes (nmo,nmo) and (nmo,nmo,nmo,nmo).
# COmpute UHF energy from integrals:
function compute_hf_energy(h0, h1_alpha, h1_beta, h2_aa, h2_bb, h2_ab,
                           Nalpha, Nbeta)

    occ_a = 1:Nalpha
    occ_b = 1:Nbeta

    # 1) One-electron part
    E1 = 0.0
    for i in occ_a
        E1 += h1_alpha[i,i]
    end
    for i in occ_b
        E1 += h1_beta[i,i]
    end

    # 2) Same-spin (direct - exchange)
    E2_same = 0.0
    for i in occ_a, j in occ_a
        E2_same += 0.5 * (h2_aa[i,i,j,j] - h2_aa[i,j,j,i])
    end
    for i in occ_b, j in occ_b
        E2_same += 0.5 * (h2_bb[i,i,j,j] - h2_bb[i,j,j,i])
    end

    # 3) Alpha-Beta
    E2_cross = 0.0
    for i in occ_a, j in occ_b
        E2_cross += h2_ab[i,i,j,j]
    end

    E_elec = E1 + E2_same + E2_cross
    return h0 + E_elec
end

"""
    test_hamiltonian_expectation(H, h0, h1_alpha, h1_beta,
                                 h2_aa, h2_bb, h2_ab,
                                 N_spatial, Nalpha, Nbeta,
                                 ket, occ, kidx; block=true)

Compare the HF energy computed from MO integrals with the PauliSum expectation value on `ket`.

- H: PauliSum Hamiltonian (what you evaluate with `expectation_value(H, ket)`)
- h0: constant term (nuclear repulsion or saved constant)
- h1_alpha, h1_beta: one-electron MO matrices (nmo×nmo)
- h2_aa, h2_bb, h2_ab: two-electron MO integrals in chemist notation (nmo×nmo×nmo×nmo)
- N_spatial: number of spatial MOs
- Nalpha, Nbeta: number of occupied alpha/beta MOs
- ket, occ, kidx: output from `string_to_ket(...)` for the same determinant
- block: true => block ordering (α then β), false => interleaved (αβ αβ ...)
"""
function test_hamiltonian_expectation(H, h0, h1_alpha, h1_beta,
                                      h2_aa, h2_bb, h2_ab,
                                      N_spatial::Int, Nalpha::Int, Nbeta::Int,
                                      ket, occ, kidx; block::Bool=true)

    # sizes
    N = N_spatial
    Ns = 2 * N

    # index mapping (must match how H was constructed)
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N
        ordering = "block (all α then all β)"
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
        ordering = "interleaved (α,β,α,β...)"
    end

    println("------------------------------------------------------------")
    @printf("Testing Hamiltonian expectation (N_spatial=%d, Nspin=%d) using %s ordering\n", N, Ns, ordering)
    println("Constructing reference HF energy from integrals...")

    # 1-electron contribution (diagonal in MO basis for a determinant)
    E1 = 0.0
    for i in 1:Nalpha
        E1 += h1_alpha[i,i]
    end
    for j in 1:Nbeta
        E1 += h1_beta[j,j]
    end

    # 2-electron same-spin (direct - exchange) with 1/2 prefactor
    E2_same = 0.0
    # alpha-alpha
    for i in 1:Nalpha, j in 1:Nalpha
        E2_same += 0.5 * ( h2_aa[i,i,j,j] - h2_aa[i,j,j,i] )
    end
    # beta-beta
    for i in 1:Nbeta, j in 1:Nbeta
        E2_same += 0.5 * ( h2_bb[i,i,j,j] - h2_bb[i,j,j,i] )
    end

    # alpha-beta cross term (no exchange, no 1/2)
    E2_cross = 0.0
    for i in 1:Nalpha, j in 1:Nbeta
        E2_cross += h2_ab[i,i,j,j]
    end

    E_elec = E1 + E2_same + E2_cross
    E_direct = h0 + E_elec

    println("Reference integrals breakdown:")
    @printf("  H0 (constant) = %20.12f\n", h0)
    @printf("  E1 (one-body diag) = %20.12f\n", E1)
    @printf("  E2_same (aa+bb direct-exchange) = %20.12f\n", E2_same)
    @printf("  E2_cross (alpha-beta) = %20.12f\n", E2_cross)
    @printf("  E_total_from_integrals = %20.12f\n", E_direct)

    # Build psi occupancy from expected Nalpha/Nbeta and the idx mapping,
    # so we can check consistency with `occ` / `kidx` returned from string_to_ket
    psi = zeros(Int, Ns)
    for i in 1:Nalpha
        psi[idx_alpha(i)] = 1
    end
    for j in 1:Nbeta
        psi[idx_beta(j)] = 1
    end

    # Convert given occ/kidx into a vector for comparison if available
    # occ is typically a list of occupied spin-orbital indices or pair indices from string_to_ket
    # kidx might encode mapping too — but just compare bit counts & distribution
    n_occ_from_psi = sum(psi)
    n_occ_from_occ = isempty(occ) ? -1 : length(occ)   # if occ empty or not provided
    @printf("Occupations check: expected (Nalpha,Nbeta) = (%d,%d)  -> total occ = %d\n", Nalpha, Nbeta, n_occ_from_psi)
    if n_occ_from_occ != -1
        @printf(" Occ list length returned by string_to_ket = %d\n", n_occ_from_occ)
    end

    # Compute Pauli expectation using existing API
    println("\nComputing expectation value from PauliSum H on provided ket ...")
    e_pauli = expectation_value(H, ket)   # uses user's API
    @printf("  PauliSum expectation (from H,ket) = %20.12f\n", e_pauli)

    # Report mismatch
    delta = e_pauli - E_direct
    @printf("\nDifference (Pauli - integrals) = %20.12e\n", delta)

    if abs(delta) < 1e-8
        println("OK: PauliSum expectation matches integrals-based energy (within tolerance).")
    else
        println("WARNING: mismatch detected. Possible culprits (in order):")
        println("  - bitstring/spin-orbital ordering mismatch (block vs interleaved)")
        println("  - when constructing H: wrong operator ordering for two-body terms or missing H0")
        println("  - fermion_op / Jordan-Wigner mapping uses a different canonical ordering")
        println("Run the per-orbital occupation checks and small operator tests if needed.")
    end

    println("------------------------------------------------------------")
    return (E_direct=E_direct, E_pauli=e_pauli, delta=delta,
            breakdown=(H0=h0, E1=E1, E2_same=E2_same, E2_cross=E2_cross))
end

data_path = "/Users/admin/PycharmProjects/pyQCTools/tests/N2_test_uhf/1.1_tensors.npz"
data = npzread(data_path)
h0 = data["hc"][1]
h1_alpha = data["h1_a"]
h1_beta = data["h1_b"]
h2_aa = data["h2aa"]
h2_bb = data["h2bb"]
h2_ab = data["h2ab"]

Nalpha = 7   # number of alpha electrons (occupied alpha MOs)
Nbeta  = 7   # number of beta electrons  (occupied beta MOs)
occ_a = 1:Nalpha
occ_b = 1:Nbeta

E_test = compute_hf_energy(h0, h1_alpha, h1_beta,
                           h2_aa, h2_bb, h2_ab,
                           Nalpha, Nbeta)

println("Energy from integrals = ", E_test)

#test_hamiltonian_expectation(H, h0, h1_alpha, h1_beta,
#                             h2_aa, h2_bb, h2_ab,
#                             10, Nalpha, Nbeta,
#                             ket, occ, kidx; block=false)

=#