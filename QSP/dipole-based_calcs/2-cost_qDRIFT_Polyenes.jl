"""
 Here we compare estimates for qDRIFT without performing calculations
    given an epsilon, compute the number of Meas/run, lambda and tau
"""
#
using LinearAlgebra
using PauliOperators
using QuantumChemQC
using Statistics
using Plots
using Printf
using Random
using StatsBase
using NPZ


"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P1-RHF_dip_mo.npz"
    data = npzread(data_path)
    d_mo = data["dip_op"]
    println("d_mo shape: ", size(d_mo))
    println(typeof(d_mo))

    N = size(d_mo)[2]  # number of spatial orbitals
    println("N:", N)
    @time D = QuantumChemQC.R_dipole_moment_op(N, data_path, block=false)
#    display(D)
    return D
end

function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P1-RHF_integrals.npz"
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
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=false)

    return H
end

# Helper: extract PauliBasis operators and coefficients from PauliSum
function extract_hamiltonian_coeffs_and_ops(H::PauliSum{N,T}) where {N,T}
    ops = PauliBasis{N}[]
    coeffs = Float64[]
    for (p, c) in H
        push!(ops, p)
        push!(coeffs, float(c))
    end
    return ops, coeffs
end

"""
    qdrift_propagator(H::PauliSum, t::Real, eps::Real;
                      seed::Union{Int,Nothing}=nothing,
                      plot::Bool=true, print_selection::Bool=true)

Perform a qDRIFT approximation for e^{i H t} where H is a PauliSum.

Returns a NamedTuple with:
  V         = simulated qDRIFT unitary (matrix)
  U_exact   = exact unitary exp(i * t * H) (matrix)
  sample    = vector of sampled indices (1-based)
  probs     = sampling distribution p_j
  λ         = sum(abs(h_j))
  τ         = per-step angle t*λ/N
  N         = number of samples used
  ops, coeffs = original lists of PauliBasis and coefficients
"""
# qdrift_propagator: use an explicit rng (or seed) for all random choices
function qdrift_estimations(ket, o::PauliSum{N,T}, H::PauliSum{N,T},
                           thresh::Real, tot_measurements::Int,
                           t::Real, eps::Real;
                           rng::AbstractRNG = Random.default_rng(),
                           seed::Union{Integer,Nothing}=nothing,
                           plot::Bool=true, print_selection::Bool=true) where {N,T}

    # If caller passed an explicit seed, override rng with a seeded MersenneTwister.
    # This makes reproducing a run trivial by passing the same seed later.
    if seed !== nothing
        rng = MersenneTwister(Int(seed))
    end

    # Extract Pauli strings and coefficients
    ops, coeffs = extract_hamiltonian_coeffs_and_ops(H)
    #@printf("Hamiltonian has %d terms \n", length(coeffs))
    L = length(coeffs)
    if L == 0
        error("Hamiltonian has no terms.")
    end

    # λ = sum_j |h_j| and probabilities p_j = |h_j| / λ
    λ = sum(abs.(coeffs))
    if λ == 0.0
        error("Sum of absolute coefficients is zero (λ=0).")
    end
    probs = abs.(coeffs) ./ λ

    # Number of samples (N); 
    Nsamples = max(1, ceil(Int, 2 * (λ * t)^2 / eps))
    #println("Number of Samples: ", Nsamples)
    # i.i.d. sampling of indices 1..L with replacement according to probs using provided rng
    sample_list = sample(rng, 1:L, Weights(probs), Nsamples; replace=true)

    # Per-step angle tau = t * λ / N
    τ = t * λ / Nsamples

    # Measurement scheduling
    N_tau = length(sample_list)
    N_meas = min(tot_measurements, N_tau)
    meas_lst = unique(round.(Int, LinRange(1, N_tau, N_meas)))
    sort!(meas_lst)
    #@printf("Will take %d measurements at indices (first 20): %s\n",
    #        length(meas_lst), string(meas_lst[1:min(end,20)]))

    if plot
        plt_bar = bar!(probs, legend=false, xlabel="Term index j", ylabel="p_j",
            title = @sprintf("qDRIFT sampling distribution (λ = %.6g)", λ))
        savefig(plt_bar, "qdrift_sampling_distribution.pdf")
    end

    # If the provider RNG was a MersenneTwister seeded from a seed variable it will be reproducible.
    # Return the seed if it is a MersenneTwister, otherwise return nothing.
    used_seed = isa(rng, MersenneTwister) ? copy(rng.seed) : nothing

    return (sample=sample_list, nsamples=Nsamples, λ=λ, τ=τ, seed=used_seed)
end

# Hamiltonian and operators
D = run_dip()
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("1100")     #P1 (4qubits)
#ket, _ = QuantumChemQC.string_to_ket("11110000")  #P2 (8qubits)
#ket, _ = QuantumChemQC.string_to_ket("111111000000")  #P3 (12qubits)
#ket, _ = QuantumChemQC.string_to_ket("1111111100000000")  #P4 (16qubits)
#ket, _ = QuantumChemQC.string_to_ket("11111111110000000000")  #P5 (20qubits)
#ket, _ = QuantumChemQC.string_to_ket("111111111111000000000000")  #P6 (24qubits)

t = 50
reps = 1 
thresh = 1e-4            # Evolution threshold for evolve
n_meas = 100             #Number of measurements

s = rand(RandomDevice(), UInt32)      # RandomDevice() uses OS entropy
seed = Int(s)
rng = MersenneTwister(seed)


eps_lst = [0.07876, 0.01, 0.1, 0.15, 0.2, 0.25, 0.30, 
           0.35, 0.4, 0.45, 0.5, 0.55, 0.60, 0.65, 0.7, 0.75, 0.8, 
           0.85, 0.9, 0.95, 1.0, 1.5, 2.0, 2.5, 3.0]
println("Eps   Meas/run   lambda   tau")
for eps in eps_lst
    # Perform qDRIF estimation
    res = qdrift_estimations(ket, D, H, thresh, n_meas, t, eps;
                                rng=rng, plot=false, print_selection=false)
    @printf("%.4f   %d    %.6f    %.6f\n", eps,  res.nsamples, res.λ, res.τ)
    
end
