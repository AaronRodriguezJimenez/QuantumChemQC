using PauliOperators
using UnitaryPruning
using QuantumChemQC
using Tullio
using Printf
using IterTools
using LinearAlgebra
using Printf


"""
In this test we will run the UnitaryPruning on the H2 Hamiltonian, the task here is to check
 is to see if the current implementation allows to perform the evolution of the Hamiltonian
 and the possible effect of a given pool of operators.
"""
function simple_pool(N) 
    # N is the number of qubits
    # Create a pool of single Pauli operators using all Z strings combinations for 
    # N qubits  
    pool = PauliSum(N)
    for i in 1:N
        for zstr in IterTools.product(i:N:N) #This will generate only Z strings of weight 1
            println(zstr)
            pool += Pauli(N, Z=collect(zstr))
        end
    end
    return pool
end

# Test the pool generation
function test_simple_pool()
    N = 4
    pool = simple_pool(N)
    @show length(pool)
    for p in pool
        display(p)
    end
end
#test_simple_pool()

"""
  Function returnign the Hamiltonian for H2 at a given interatomic distance
    d -> The interatomic distance
"""
function H2_hamiltonian_diss(d)
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
return gen, coeff, N;
end

"""
 The following function computes the gradients of the energy computed as the expectation value of a commutator
  < [H, A] > with respect to a state |psi>
  H -> Hamiltonian as a PauliSum
  A -> Operator as a Pauli
  psi -> State as a Ket
"""
function energy_gradient(H::PauliSum{N}, A::PauliSum{N}, psi::Ket) where N
    # collect our results here...
    expval = zero(ComplexF64)
    commutator = PauliSum(N)

    for (h, c) in H
        for (a, d) in A
            if !PauliOperators.commute(h, a)
                comm = 2 * h * a # This includes the correct ±i phase
                println("Commutator of $h and $a is: ")
                display(comm)
                sum!(commutator, c * d * comm)
            end
        end
    end

    println("The full commutator is: ")
    display(commutator)
    println("- - - - - - - - - - - - - ")
    println("Now computing the expectation value of the commutator...")
    for (oi, coeff) in commutator
        expval += coeff * PauliOperators.expectation_value(oi, psi)
        println("Expectation value of the commutator is: $expval")
    end
    return imag(expval)
end

#- - - Test UnitaryPrunning ADAPT-VQE BASED
function run(; d = 0.74, w_type = "Majorana", max_weight=1)

    H_term, H_coeff, N_qubits = H2_hamiltonian_diss(d)
    ket = Ket(N_qubits, 3) #HF STATE
    display(ket)

    #Initial operator and angle
    generators = Vector{Pauli{N_qubits}}()
    parameters = Vector{Float64}()

    o_init = Pauli(N_qubits, X=[1,4]) #Initial operator for evolution
    angle_init = 1.0
    push!(generators, o_init)
    push!(parameters, angle_init)

    #pool = simple_pool(N_qubits)
    pool = PauliSum(N_qubits)
    pool += angle_init * o_init

    println("Selected pool is: ")
    display(pool)


    #println("Hamiltonian for H2 at d = $d :")
    #for (pauli, coeff) in zip(H_term, H_coeff)
    #    println(@sprintf("%+.6f %+.6fim | %s", real(coeff), imag(coeff), string(pauli)))
    #end

    println("- - - - - - - - - - - - - - - - -")
    hamiltonian = PauliSum(N_qubits)
    for (p,c) in zip(H_term, H_coeff)
        hamiltonian += c * p
        #push!(generators, p)
        #push!(parameters, c)
    end
    println("Hamiltonian as PauliSum:")
    display(hamiltonian)

    #Call to bfs bfs_evolution and perform initial propagation of the Hamiltonian
    evol_op = UnitaryPruning.bfs_evolution_vqe(generators, parameters, hamiltonian, thresh=1e-3)
    println("Evolved Hamiltonian after Unitary Pruning:")
    display(evol_op)

    # Let's add a new operator to the pool
    pool += 1.0 * Pauli(N_qubits, X=[1], Z=[2]) 

    #Compute the energy gradient with respect to a pool of operators
    gradient = energy_gradient(evol_op, pool, ket)
    println("Energy gradient with respect to the initial operator is: $gradient")

    # Compute energy of the transformed Hamiltonian
    energy = PauliOperators.expectation_value(evol_op, ket)
    println("CYCLE #1, Energy of the transformed Hamiltonian is: $energy")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    #Test GSD pool
    gsd_pool = ADAPT_PP.base.qubitexcitationpool(N_qubits)
    println("GSD pool generated with $(length(gsd_pool)) operators.")
    
    gradients = Vector{ComplexF64}()
    new_generators = Vector{Pauli{N_qubits}}()
    new_coeffs = Vector{Float64}()
    for (i, op) in enumerate(gsd_pool)
        println("Operator $i:")
        display(op)
        A = PauliSum(N_qubits)
        A += 1.0 * op


        #println("Target and source orbitals: ", target_and_source)
        println("- - - - - - - - - - - - - - - - -")
        gradient = energy_gradient(evol_op, A, ket)
        println("Energy gradient with respect to the initial operator is: $gradient")
        if gradient != 0.0
            println("Non-zero gradient found for operator $i:")
            display(op)
            push!(gradients, gradient)
            push!(new_generators, op)
            push!(new_coeffs, 1.0)
        end
        
    end
    println("All gradients computed:")
    for (i, grad) in enumerate(gradients)
        println("Gradient $i: $grad")
    end

    println("New pool with non-zero gradients has $(length(new_generators)) operators.")
    println("Operators in the new pool:")
    for (i, op) in enumerate(new_generators)
        println("Operator $i:")
        display(op)
        push!(gsd_pool, op) #Add selected operator to the pool
    end
    #Redo evolution with one od the new operators
    evol_op = UnitaryPruning.bfs_evolution_vqe(new_generators, new_coeffs, evol_op, thresh=1e-3)
    println("Evolved Hamiltonian after Unitary Pruning:")
    display(evol_op)
    
    # Compute energy of the transformed Hamiltonian
    energy = PauliOperators.expectation_value(evol_op, ket)
    println("CYCLE #2, Energy of the transformed Hamiltonian is: $energy")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    # Another cycle of gradient evaluation`
    gradients = Vector{ComplexF64}()
    new_generators = Vector{Pauli{N_qubits}}()
    new_coeffs = Vector{Float64}()
    for (i, op) in enumerate(gsd_pool)
        println("Operator $i:")
        display(op)
        A = PauliSum(N_qubits)
        A += 1.0 * op


        #println("Target and source orbitals: ", target_and_source)
        println("- - - - - - - - - - - - - - - - -")
        gradient = energy_gradient(evol_op, A, ket)
        println("Energy gradient with respect to the initial operator is: $gradient")
        if gradient != 0.0
            println("Non-zero gradient found for operator $i:")
            display(op)
            push!(gradients, gradient)
            push!(new_generators, op)
            push!(new_coeffs, 1.0)
        end
        
    end
    println("All gradients computed:")
    for (i, grad) in enumerate(gradients)
        println("Gradient $i: $grad")
    end

    println("New pool with non-zero gradients has $(length(new_generators)) operators.")
    println("Operators in the new pool:")
    for (i, op) in enumerate(new_generators)
        println("Operator $i:")
        display(op)
    end
    #Redo evolution with one od the new operators
    evol_op = UnitaryPruning.bfs_evolution_vqe(new_generators, new_coeffs, evol_op, thresh=1e-3)
    println("Evolved Hamiltonian after Unitary Pruning:")
    display(evol_op)
    
    # Compute energy of the transformed Hamiltonian
    energy = PauliOperators.expectation_value(evol_op, ket)
    println("CYCLE #3, Energy of the transformed Hamiltonian is: $energy")

end 

run(d = 0.740,  w_type = "Pauli", max_weight=4)

