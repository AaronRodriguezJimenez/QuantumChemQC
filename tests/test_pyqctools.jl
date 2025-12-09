using PyCall
using PauliOperators
using DBF
using Printf
"""
 Here we call to pyqctools, which computes a molecular Hamiltonian and 
 returns it to Julia. Then we use PauliOperators and DBF to compute its 
 ground state energy.
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
    geom = pq.geometries.H2()
    #geom = pq.geometries.LiH()
    #geom = pq.geometries.BeH2()
    #geom = pq.geometries.H2O()
#    geom = pq.ethylene()

    println("Molecule geometry from pyqctools:", geom)
    
    bas = "cc-pVTZ" #"sto-3g"
    mult = 1
    ch = 0
    H_op, N_total = Hfcns.get_hamiltonian(geom, bas, mult, ch, print_table=false)
    #println("Hamiltonian operator from pyqctools:", H_op)
    println("Number of qubits:", N_total)

    # Transform to qiskit format
    Hqiskit = Hfcns.jw_to_sparse_pauli_op(H_op)
    #println("Hamiltonian in Qiskit format:", Hqiskit)
    println("Number of Pauli terms:", length(Hqiskit))

    # Construct the QUBIT Hamiltonian with PauliOperators
    gens = Vector{PyObject}()
    coeffs = Vector{ComplexF64}()
    for (i,term) in enumerate(Hqiskit.paulis)
        #println("Term: ", term, " Coefficient: ", Hqiskit.coeffs[i]) 
        push!(gens, term)
        push!(coeffs, Hqiskit.coeffs[i])
    end

    H = PauliSum(N_total, Float64)
    for (gens, coeff) in zip(gens, coeffs)
        pstring = String(gens.to_label())
        H += real(coeff) * Pauli(pstring)
    end

    #println("Hamiltonian in PauliOperators format:", H)
    #display(H)
    #println("N_total = ", N_total)
    return H, N_total
end

# Now use DBF to compute ground state energy
function dbf_gstate(H, N)

    Nparticles = 2 # total electrons in molecule
    ket, occ = DBF.particle_ket(N, Nparticles, 0.0; mode=:first)
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
                     max_iter=400, conv_thresh=1e-4, 
                    evolve_coeff_thresh=1e-3,
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

dbf_gstate(H, N_total)