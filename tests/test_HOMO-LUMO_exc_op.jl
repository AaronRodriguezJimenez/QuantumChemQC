using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots

"""
 Here we test the implementation/usage of the HOMO-LUMO excitation operator as evolved under the molecular Hamiltonian.
"""
function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/dipole_moment/h2-RHF_tensors.npz"
    data = npzread(data_path)
    H0 = data["hc"][1]
    H1 = data["h1e"]
    H2 = data["h2e"]

    println("H0 shape: ", size(H0))
    println(typeof(H0))
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))

    Norbs = size(H1,1)  # number of spatial orbitals
    
    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=true)

    return H
end

function get_N(d::Dict{PauliBasis{N}, T}) where {N, T}
    return N
end

#- - - Get Molecular Hamiltonian
H = run_H()
N = get_N(H)
println("N:", N)

#- - - Define Reference ket
ket, _ = QuantumChemQC.string_to_ket("1010")

function occ_operator(N_spin_orbitals; block=true)
    coeff_thresh_clip = 1e-6
    N_spatial = Int(N_spin_orbitals/2)

    #N = number of spin orbitals
    # 1. Pre-calculate Operators
    # We store these to avoid re-allocating Pauli objects constantly.
    println("  -> Generating Operator Cache...")
    ops_dag = [QuantumChemQC.fermion_op(N_spin_orbitals, i, dagger=true) for i in 1:N_spin_orbitals]
    ops_col = [QuantumChemQC.fermion_op(N_spin_orbitals, i, dagger=false) for i in 1:N_spin_orbitals]

     # --- INDEX MAPPING: block vs interleaved ---
    if block
        idx_alpha = p -> p
        idx_beta  = p -> p + N_spatial
    else
        idx_alpha = p -> 2*p - 1
        idx_beta  = p -> 2*p
    end

    # 2. Construct Operator
    println("  -> Constructing Operator...")
    NOp =  PauliSum(N_spin_orbitals, ComplexF64)
    
    for p in 1:N_spatial
        NOp += ops_dag[idx_alpha(p)] * ops_col[idx_alpha(p)]
        NOp += ops_dag[idx_beta(p)] * ops_col[idx_beta(p)]
    end

    coeff_clip!(NOp, thresh=coeff_thresh_clip)

    return NOp
end

O = occ_operator(N, block=true)
display(O)

println("Initial eigenstate")
display(ket)
Oket = O * ket
println("- - - N|1010> = $Oket")

#=
println("Initial eigenstate |0> ")
display(ket)
println("- - - <0|O|0> - - - -")
function compute_ref_expval(H, ket)
    ref_expval = 0
    for (p, c) in O
        ref_expval += c * expectation_value(p, ket)
    end
    return ref_expval
end

expval_O = compute_ref_expval(O, ket)
=#

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
       # hl_beta = find_homo_lumo(idx_beta)
       # if hl_beta !== nothing
       #     h, l = hl_beta
       #     sum!(O, ops_dag[idx_beta(l)] * ops_col[idx_beta(h)])
       #     sum!(O, ops_dag[idx_beta(h)] * ops_col[idx_beta(l)]) #h.c
       # end

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

println("* * * * * * P2 * * * * * * ")
println("# HOMO-LUMO Excitation Operator (Spin-conserving)")
ket, _ = QuantumChemQC.string_to_ket("11110000") #Block encoding
O = homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
ket_2 = O*ket
display(ket_2)
println("# HOMO-LUMO Excitation Operator (Spin-non-conserving)")    
O = homo_lumo_excitation_op(ket, spin_conserving=false, block=false)
display(O)
println(O*ket)

println("* * * * * * P3 * * * * * * ")
println("# HOMO-LUMO Excitation Operator (Spin-conserving)")
ket, _ = QuantumChemQC.string_to_ket("111111000000") #Block encoding
O = homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
ket_2 = O*ket
display(ket_2)
println("# HOMO-LUMO Excitation Operator (Spin-non-conserving)")    
O = homo_lumo_excitation_op(ket, spin_conserving=false, block=false)
display(O)
println(O*ket)

println("* * * * * * P4 * * * * * * ")
println("# HOMO-LUMO Excitation Operator (Spin-conserving)")
ket, _ = QuantumChemQC.string_to_ket("1111111100000000") #Block encoding
O = homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
ket_2 = O*ket
display(ket_2)
println("# HOMO-LUMO Excitation Operator (Spin-non-conserving)")    
O = homo_lumo_excitation_op(ket, spin_conserving=false, block=false)
display(O)
println(O*ket)

println("* * * * * * P5 * * * * * * ")
println("# HOMO-LUMO Excitation Operator (Spin-conserving)")
ket, _ = QuantumChemQC.string_to_ket("11111111110000000000") #Block encoding
O = homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
ket_2 = O*ket
display(ket_2)
println("# HOMO-LUMO Excitation Operator (Spin-non-conserving)")    
O = homo_lumo_excitation_op(ket, spin_conserving=false, block=false)
display(O)
println(O*ket)

println("* * * * * * P6 * * * * * * ")
println("# HOMO-LUMO Excitation Operator (Spin-conserving)")
ket, _ = QuantumChemQC.string_to_ket("111111111111000000000000") #Block encoding
O = homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
ket_2 = O*ket
display(ket_2)
println("# HOMO-LUMO Excitation Operator (Spin-non-conserving)")    
O = homo_lumo_excitation_op(ket, spin_conserving=false, block=false)
display(O)
println(O*ket)