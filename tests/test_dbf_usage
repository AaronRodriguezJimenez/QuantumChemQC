using QuantumChemQC
using DBF
using PauliOperators
using Tullio
using Printf
using IterTools
using DBF
using Random
using LinearAlgebra
using JLD2

function run()
    """
     H2 molecule ground staste using DBF .
    """
    #Initialize molecular Hamiltonian
    d = 0.740  # interatomic distance in Angstrom
    H, N = QuantumChemQC.H2_hamiltonian(d)
    
    display(H)
    # - - - Using DBF - - - 
    #
    """
     Here we perform a test calculation using an approximation to the projector: single particle, two particle, etc.
    """

    ket = Ket(N, 3) #HF STATE
    display(ket)
    occ = [1,1,0,0]  #Occupation vector for HF state
    # Transform H to make |000> the most stable bitstring
    for i in 1:N
        if occ[i] == 1
            H = Pauli(N, X=[i]) * H * Pauli(N, X=[i])
        end
    end
    
    H0 = deepcopy(H)
    println(" Transformed H:")
    display(H)
    Hmat = Matrix(H)
    evals = eigvals(Hmat)
    @show minimum(evals)

    ψ = Ket([0 for i in 1:N])
    display(ψ)

    e0 = expectation_value(H,ψ)
    @printf(" Reference = %12.8f\n", e0)

    g = Vector{PauliBasis{N}}([]) 
    θ = Vector{Float64}([]) 

    println("\n ########################")
    @time H, g, θ = DBF.groundstate_diffeq_test(H0, ψ, n_body=1, 
                                verbose=1, 
                                max_iter=1000, conv_thresh=1e-2, 
                                evolve_coeff_thresh=1e-8,
                                grad_coeff_thresh=1e-8,
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
   

end


run()