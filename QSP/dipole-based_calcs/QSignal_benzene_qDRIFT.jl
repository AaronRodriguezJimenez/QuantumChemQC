using QuantumChemQC
using PauliOperators
using Printf
using NPZ
using LinearAlgebra
using Statistics
using Random
using StatsBase
using Plots

"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/benzene_and_acenes/benz-RHF_dip_mo.npz"
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
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/benzene_and_acenes/benz-RHF_integrals.npz"
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
    return ket, out
end

"""
 Evolve function from DBF code
"""
function evolve!(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
    _cos = cos(θ)
    _sin = 1im*sin(θ)
    sin_branch = PauliSum(N)
    for (p,c) in O
        if PauliOperators.commute(p,G) == false
            # replace sum! with more efficient version
            # sum!(sin_branch, c*_sin*G*p)
            tmp = c*_sin*G*p
            curr = get(sin_branch, PauliBasis(tmp), 0.0) + PauliOperators.coeff(tmp)
            sin_branch[PauliBasis(tmp)] = curr 
            O[p] *= _cos
        end
    end
    sum!(O, sin_branch)
    return O 
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
function qdrift_propagator(ket, o::PauliSum{N,T}, H::PauliSum{N,T},
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
    @printf("Hamiltonian has %d terms \n", length(coeffs))
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
    println("Number of Samples: ", Nsamples)
    # i.i.d. sampling of indices 1..L with replacement according to probs using provided rng
    sample_list = sample(rng, 1:L, Weights(probs), Nsamples; replace=true)

    # Per-step angle tau = t * λ / N
    τ = t * λ / Nsamples
    println(τ)

    # Prepare storage
    Wt = deepcopy(o)
    W  = deepcopy(o)
    rCtvals = Float64[]
    iCtvals = Float64[]

    # Measurement scheduling
    N_tau = length(sample_list)
    N_meas = min(tot_measurements, N_tau)
    meas_lst = unique(round.(Int, LinRange(1, N_tau, N_meas)))
    sort!(meas_lst)
    #@printf("Will take %d measurements at indices (first 20): %s\n",
    #        length(meas_lst), string(meas_lst[1:min(end,20)]))

    time_grid = collect(range(0.0, stop=t, length=N_meas))

    # Iterate qDRIFT steps
    meas_idx = 1
    next_meas = meas_lst[meas_idx]
    sample_idx = 1
    for idx in sample_list
        Pi = ops[idx]
        s = sign(coeffs[idx])
        theta = 2 * τ * s
        pb = PauliBasis(Pi)
        evolve!(Wt, pb, theta)
        coeff_clip!(Wt, thresh=thresh)
        WWt = W * Wt
        
        if sample_idx == next_meas
            expval = expectation_value(WWt, ket)
            push!(rCtvals, real(expval))
            push!(iCtvals, imag(expval))

            meas_idx += 1
            if meas_idx <= length(meas_lst)
                next_meas = meas_lst[meas_idx]
            else
                next_meas = typemax(UInt128)
            end
        end
        sample_idx += 1
    end

    if print_selection
        println("=== qDRIFT Selection (first 200 shown) ===")
        println(sample_list[1:min(end,200)])
        println("... (total samples = $Nsamples)")
    end

    if plot
        plt_bar = bar!(probs, legend=false, xlabel="Term index j", ylabel="p_j",
            title = @sprintf("qDRIFT sampling distribution (λ = %.6g)", λ))
        savefig(plt_bar, "qdrift_sampling_distribution.pdf")
    end

    # If the provider RNG was a MersenneTwister seeded from a seed variable it will be reproducible.
    # Return the seed if it is a MersenneTwister, otherwise return nothing.
    used_seed = isa(rng, MersenneTwister) ? copy(rng.seed) : nothing

    return (sample=sample_list, probs=probs, time_grid=time_grid, nsamples=N_tau,
            λ=λ, τ=τ, N=Nsamples, RCt=rCtvals, ICt=iCtvals, seed=used_seed)
end

"""
 averaged_qdrift: create a fresh RNG per run (random, independent selections each call)
"""
function averaged_qdrift(n_runs, ket, o::PauliSum{N,T}, H::PauliSum{N,T},
                         thresh::Real, tot_measurements::Int,
                         t::Real, eps::Real) where {N,T}
    sum_RCt = zeros(tot_measurements)
    sum_ICt = zeros(tot_measurements)
    t_grid = zeros(tot_measurements)

    seeds = Vector{Union{Int,Nothing}}(undef, n_runs)  # to store seeds used per run (for reproducibility)
    elapsed_times = Vector()

    for r in 1:n_runs
        # Get a fresh, unpredictable 32-bit seed from OS entropy
        s = rand(RandomDevice(), UInt32)      # RandomDevice() uses OS entropy
        seed = Int(s)
        rng = MersenneTwister(seed)          # rng for this run
        t1 = time()
        res = qdrift_propagator(ket, o, H, thresh, tot_measurements, t, eps;
                                rng=rng, plot=false, print_selection=false)
        
        elapsed_time = time() - t1
        #println("Elapsed time: ", elapsed_time, " seconds")
        push!(elapsed_times, elapsed_time)
        
        sum_RCt .+= res.RCt
        sum_ICt .+= res.ICt
        t_grid = res.time_grid
        seeds[r] = seed   
        @printf("Seed: %d  Run: %d  Samples:  %d   Time:  %12.8f ", seed, r,  res.nsamples, elapsed_time)

        rRES = sum_RCt ./ r
        iRES = sum_ICt ./ r

        open("/Users/admin/VSCProjects/QuantumChemQC/Benzene_qdrift__Q_signals_rep_$r.txt", "w") do io
            println(io, "- - - Benzene (6e,6o) output qDRIFT signals, epsilon= $eps  - - - ")
            @printf(io, "Seed: %d  Run: %d  Samples:  %d   Time:  %12.8f\n", seed, r,  res.nsamples, elapsed_time)
            @printf(io, "dt   Re(C(t))    Im(C(t))\n")
            for (i, interval) in enumerate(t_grid)
                @printf(io, "%.4f     %.6f    %.6f\n", interval, rRES[i], iRES[i])
            end        
        end
    end
    
    println("Elapsed times") 
    for (i,t) in enumerate(elapsed_times)
        println(i," ", t)
    end


    return 
end

#########################################
D = run_dip()
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("111111000000")

t = 50
eps = 3.0
reps = 25 
thresh = 1e-5 # Evolution threshold for evolve
n_meas = 100 #Number of measurements

averaged_qdrift(reps, ket, D, H, thresh, n_meas, t, eps)