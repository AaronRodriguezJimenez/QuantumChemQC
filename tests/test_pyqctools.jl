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


    # Get geometry for H2O molecule and define parameters
    geom = pq.geometries.H2O()

    println("Water geometry from pyqctools:", geom)
    
    bas = "sto-3g"
    mult = 1
    ch = 0
    H_op, N_total = Hfcns.get_hamiltonian(geom, bas, mult, ch, print_table=false)
    println("Hamiltonian operator from pyqctools:", H_op)
    println("Number of qubits:", N_total)

    # Transform to qiskit format
    Hqiskit = Hfcns.jw_to_sparse_pauli_op(H_op)
    println("Hamiltonian in Qiskit format:", Hqiskit)

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

    println("Hamiltonian in PauliOperators format:", H)
    display(H)

    return H, N_total
end

# Now use DBF to compute ground state energy
function dbf_gstate(H, N)

    Nparticles = 10 # total electrons in H2O
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

    e0 = expectation_value(H,ψ)
    @printf(" Reference = %12.8f\n", e0)

    g = Vector{PauliBasis{N}}([]) 
    θ = Vector{Float64}([]) 

    println("\n ########################")
    @time H, g, θ = DBF.groundstate_diffeq_test(H0, ψ, n_body=3, 
                                verbose=1, 
                                max_iter=1000, conv_thresh=1e-2, 
                                evolve_coeff_thresh=1e-2,
                                grad_coeff_thresh=1e-4,
                                stepsize=.01)
    
    # @save "out_$(i).jld2" N ψ H0 H g θ
    return

    println("\n Now reroptimize with higher accuracy:")
    @show length(θ)
    Ht = deepcopy(H0)
    err = 0
    ecurr = expectation_value(Ht,ψ)
    @printf(" Initial energy: %12.8f %8i\n", ecurr, length(Ht))
    for (i,gi) in enumerate(g)
            
        θj, costi = DBF.optimize_theta_expval_test(Ht, gi, ψ, verbose=0)
        Ht = DBF.evolve(Ht, gi, θj)
        θ[i] = θj
        
        e1 = expectation_value(Ht,ψ)
        DBF.coeff_clip!(Ht, thresh=1e-5)
        e2 = expectation_value(Ht,ψ)

        err += e2 - e1
        if i%100 == 0
            @printf(" Error: %12.8f\n", err)
            e0, e2 = DBF.pt2(Ht, ψ)
            @printf(" E0 = %12.8f E2 = %12.8f EPT2 = %12.8f \n", e0, e2, e0+e2)
            e0, e, v, basis = DBF.cepa(Ht, ψ, thresh=1e-6, tol=1e-2, verbose=0)
            e0, e, v, basis = DBF.fois_ci(Ht, ψ, thresh=1e-6, tol=1e-2, verbose=0)
        end
    end    
    ecurr = expectation_value(Ht,ψ)
    @printf(" ecurr %12.8f err %12.8f %8i\n", ecurr, err, length(Ht))
    return 
end

H, N_total = run()

dbf_gstate(H, N_total)