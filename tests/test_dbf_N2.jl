using PyCall
using PauliOperators
using DBF
using Printf
using QuantumChemQC
"""
 Here we call to pyqctools, to compute the molecular integrals for 
 constructing the Hamiltonian of ethylene. Then we use PauliOperators 
 and DBF to compute its ground state energy.
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
    geom = pq.geometries.N2()    

    println("Molecule geometry from pyqctools:", geom)
    
    # GET MOLECULAR ORBITALS
    mol = pyimport("pyscf.gto").Mole()
    mol.atom = geom
    mol.basis = "sto-3g"
    mol.build()    
    #println("Nao:")
    #println(mol.nao_nr())

    # RUN SCF
    mf = pyimport("pyscf.scf").RHF(mol).run()
    E_HF = mf.e_tot
    println("Hartree-Fock Energy:", E_HF)
    
    # Ger spinorbitals tensors from pyqctools
    H0, H1, H2 = Hfcns.get_tensors(mol, mf, localized=false)
    n = size(H1,1)  # number of spin orbitals
    N_total = n

    H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    return H, N_total
end

#- - - Helpers for getting the occupation from bitstring - - - -
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
function dbf_gstate(H, N)

    Nparticles = 14 # total electrons in molecule
    ket, occ = DBF.particle_ket(N, Nparticles, 0.0; mode=:first)
    
    #kidx = argmin([real(expectation_value(H,Ket{N}(ψi))) for ψi in 1:2^N])
    #ket = Ket{N}(kidx)
    #occ = occvec(ket)

    println("Initial state:")
    display(ket)

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
                     max_iter=400, conv_thresh=1e-5, 
                    evolve_coeff_thresh=1e-4,
                    grad_coeff_thresh=1e-10,
                    energy_lowering_thresh=1e-10)                         
    
    H = res["hamiltonian"]
    gi = res["generators"]
    θi = res["angles"]    

    e2 = real(expectation_value(H,ψ))
    @printf("<H> = %12.8f <U'HU> = %12.8f \n", e1, e2)

    # Diplay evolved Hamiltonian
    #println("\n Evolved Hamiltonian:")
    #display(H)
                    
    return
end

H, N_total = run()

#display(H)
println("Number of qubits: ", N_total)
println("Input Hamiltonian #terms: ", length(H))

dbf_gstate(H, N_total)

# Initial energy
#kidx = argmin([real(expectation_value(H,Ket{N_total}(ψi))) for ψi in 1:2^N_total])
#ψ = Ket{N_total}(kidx)
#println("Index of min energy bitstring: ", kidx)
#display(ψ)  
#e0 = real(expectation_value(H, ψ))
#@printf(" Initial <H> = %12.8f \n", e0)
#println("Occuation vector: ", occvec(ψ))