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

    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/VSCProjects/QuantumChemQC/tests/1.5_tensors.npz"
    data = npzread("/Users/admin/VSCProjects/QuantumChemQC/tests/1.5_tensors.npz")
    H0 = 0.00# data["hc"][1]
    H1 = data["h1e"]
    H2 = data["h2e"]

    println("H0 shape: ", size(H0))
    println(typeof(H0))
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))

    N_spin_orbitals = size(H1,1)  # number of spin orbitals

    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(N_spin_orbitals, data_path)
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
    max_iter=10
    conv_thresh=1e-5
    evolve_coeff_thresh=1e-4
    grad_coeff_thresh=1e-10
    energy_lowering_thresh=1e-10

    # = = = Set Reference State = = = 
    Nparticles = 10 # total electrons in active space
    #ket, occ = DBF.particle_ket(N, Nparticles, 0.0; mode=:first)
    #kidx = argmin([real(expectation_value(H, Ket{N}(ψi))) for ψi in 1:2^N])
    #ket = Ket{N}(kidx)
    #occ = occvec(ket)
    ket, occ, kidx  = string_to_ket("1111111111000000") #Leading CAS configuration

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

    #g = Vector{PauliBasis{N}}([]) 
    #θ = Vector{Float64}([]) 

    println("\n ########################")
    res = DBF.dbf_groundstate(H0, ψ,
                    max_iter= max_iter,
                    conv_thresh= conv_thresh, 
                    evolve_coeff_thresh= evolve_coeff_thresh,
                    grad_coeff_thresh= grad_coeff_thresh,
                    energy_lowering_thresh= energy_lowering_thresh,
                    checkfile = "Hchkf")   
    
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
                    
    return
end

H = run()
println(typeof(H))
dbf_gstate(H)