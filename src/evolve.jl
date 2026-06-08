using JLD2

"""
 # # # Operator evolution related functions
"""
function evolve!(O::PauliSum{N,T}, G::PauliBasis{N}, θ::Real) where {N,T}
    cθ = cos(θ)
    sθ = 1im * sin(θ)
    added = PauliSum(N)

    for (p, c) in O
        if !PauliOperators.commute(p, G)
            tmp = c * sθ * G * p
            key = PauliBasis(tmp)
            added[key] = get(added, key, 0.0) + PauliOperators.coeff(tmp)
            O[p] *= cθ
        end
    end

    sum!(O, added)
    return O
end

function old_evolve!(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
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

"""
 Converts the PauliSum Operator into generators and angles vectors
 H :: A PauliSum Hamiltonian
"""
function gens_from_H(H::PauliSum{N,T}) where {N,T}
    generators = Vector{PauliBasis{N}}([])
    angles = Vector{Float64}([])
    for (p,c) in H
        #println(c)
        push!(generators, p)
        push!(angles, c)
    end
    return generators, angles;
end


"""
Pauli-propagation time series with pruning controls.
Using First order Trotterization
THIS FUNCTION INCLUDES THE PRUNING ERROR CORRECTION.
Arguments:
- generators, angles  : from UnitaryPruning.heisenberg_1D(...)
- o                   : PauliSum (used as both V and W by default)
- ket                 : state for estimator (your current approach)
- thresh              : magnitude threshold (|coeff|)
- n_intervals : number of intervals,  n_intervals =tot_time/dt
"""
function QSP_evolution_op(ket, o::PauliSum{N,T}, 
                          H::PauliSum{N,T}, n_intervals, dt;
                          thresh::Float64=1e-3,
                          max_weight=N,
                          checkfile=nothing) where {N,T}

    Wt = deepcopy(o)
    W  = deepcopy(o)

    #Preallocate
    rCtvals = Vector{Float64}(undef, n_intervals + 1)
    iCtvals = Vector{Float64}(undef, n_intervals + 1)

    out = Dict(
        "ReCt" => Vector{Float64}(undef, n_intervals + 1),
        "ImCt" => Vector{Float64}(undef, n_intervals + 1),
        "intervals" => Vector{Float64}(undef, n_intervals + 1),
        "norms" => Vector{Float64}(undef, n_intervals + 1),
    )

    # t = 0
    WW = W * W
    expval0 = expectation_value(WW, ket)

    rCtvals[1] = real(expval0)
    iCtvals[1] = imag(expval0)

    out["ReCt"][1] = rCtvals[1]
    out["ImCt"][1] = iCtvals[1]
    out["intervals"][1] = 0.0
    out["norms"][1] = norm(Wt)

    # Precompute generators
    generators, angles = gens_from_H(H)
    nt = length(angles)

    # Precompute PauliBasis
    pbs = [PauliBasis(g) for g in generators]

    println("Total Rotations:", nt * n_intervals)

    for ki in 1:n_intervals
        accumulated_error = 0.0

        for j in 1:nt
            theta = 2 * dt * angles[j]
            pb = pbs[j]

            evolve!(Wt, pb, theta)

            coeff_clip!(Wt, thresh=1.0e-12)

            WWt = W * Wt
            e1 = expectation_value(WWt, ket)

            coeff_clip!(Wt, thresh=thresh)
            weight_clip!(Wt, max_weight)

            # recompute after pruning
            WWt = W * Wt
            e2 = expectation_value(WWt, ket)

            accumulated_error += e2 - e1
        end

        # final value
        WWt = W * Wt
        expval = expectation_value(WWt, ket) + accumulated_error

        idx = ki + 1

        rCtvals[idx] = real(expval)
        iCtvals[idx] = imag(expval)

        out["ReCt"][idx] = rCtvals[idx]
        out["ImCt"][idx] = iCtvals[idx]
        out["intervals"][idx] = ki * dt
        out["norms"][idx] = norm(Wt)

        # Checkpointing every 50 intervals 
        if checkfile !== nothing && (ki % 50 == 0)
            @save "$(checkfile)_checkpoint.jld2" out ki
        end
    end

    # final save
    if checkfile !== nothing
        @save "$(checkfile).jld2" out
    end

    tgrid = range(0.0, stop=n_intervals*dt, length=n_intervals+1)

    return rCtvals, iCtvals, collect(tgrid)
end

"""
Pauli-propagation time series with pruning controls.
Using First order Trotterization
THIS FUNCTION DO NOT INCLUDE THE PRUNING ERROR CORRECTION.
Arguments:
- generators, angles  : from UnitaryPruning.heisenberg_1D(...)
- o                   : PauliSum (used as both V and W by default)
- ket                 : state for estimator (your current approach)
- thresh              : magnitude threshold (|coeff|)
- n_intervals : number of intervals,  n_intervals =tot_time/dt
"""
function QSP_evolution_op_no_corr(ket, o::PauliSum{N,T}, 
                                  H::PauliSum{N,T}, n_intervals, dt;
                                  thresh::Float64=1e-3,
                                  max_weight=N,
                                  checkfile=nothing) where {N,T}
    
    Wt = deepcopy(o)
    W  = deepcopy(o)

    # Preallocations
    rCtvals = Vector{Float64}(undef, n_intervals + 1)
    iCtvals = Vector{Float64}(undef, n_intervals + 1)

    out = Dict(
        "ReCt" => Vector{Float64}(undef, n_intervals + 1),
        "ImCt" => Vector{Float64}(undef, n_intervals + 1),
        "intervals" => Vector{Float64}(undef, n_intervals + 1),
        "norms" => Vector{Float64}(undef, n_intervals + 1),
    )

    # t = 0
    WW = W * W
    expval0 = expectation_value(WW, ket)

    rCtvals[1] = real(expval0)
    iCtvals[1] = imag(expval0)

    out["ReCt"][1] = rCtvals[1]
    out["ImCt"][1] = iCtvals[1]
    out["intervals"][1] = 0.0
    out["norms"][1] = norm(Wt)

    # Precompute generators
    generators, angles = gens_from_H(H)
    nt = length(angles)

    # Precompute PauliBasis
    pbs = [PauliBasis(g) for g in generators]

    println("Total Rotations:", nt * n_intervals)

    for ki in 1:n_intervals
        accumulated_error = 0.0

        for j in 1:nt
            theta = 2 * dt * angles[j]
            pb = pbs[j]
            evolve!(Wt, pb, theta)
            coeff_clip!(Wt, thresh=thresh)
            weight_clip!(Wt, max_weight)
        end

        # final value
        WWt = W * Wt
        expval = expectation_value(WWt, ket) #+ accumulated_error  # NO CORRECTION
        idx = ki + 1

        rCtvals[idx] = real(expval)
        iCtvals[idx] = imag(expval)

        out["ReCt"][idx] = rCtvals[idx]
        out["ImCt"][idx] = iCtvals[idx]
        out["intervals"][idx] = ki * dt
        out["norms"][idx] = norm(Wt)

        # Checkpointing every 50 intervals 
        if checkfile !== nothing && (ki % 50 == 0)
            @save "$(checkfile)_nocorr_checkpoint.jld2" out ki
        end
    end

    # final save
    if checkfile !== nothing
        @save "$(checkfile)_nocorr.jld2" out
    end

    tgrid = range(0.0, stop=n_intervals*dt, length=n_intervals+1)

    return rCtvals, iCtvals, collect(tgrid)
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
function qdrift_propagator(ket, o::PauliSum{N,T}, H::PauliSum{N,T}, 
                           thresh::Real, tot_measurements::Int,
                           t::Real, eps::Real;
                           seed::Union{Int,Nothing}=nothing,
                           plot::Bool=true, print_selection::Bool=true) where {N,T}

    # Extract Pauli strings and coefficients
    ops, coeffs = gens_from_H(H)
    @printf("Hamiltonian has %.2f terms \n", length(coeffs))
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
    Nsamples = ceil(Int, 2 * (λ * t)^2 / eps)
    if Nsamples < 1
        Nsamples = 1
    end

    # Set RNG seed if provided
    if seed !== nothing
        Random.seed!(seed)
    end

    # i.i.d. sampling of indices 1..L with replacement according to probs
    sample_list = sample(1:L, Weights(probs), Nsamples; replace=true)

    # Per-step angle tau = t * λ / N
    τ = t * λ / Nsamples

    #
    #* * * QSP setup * * *
    #
    Wt = deepcopy(o)        # evolve W ≡ U*OU
    W  = deepcopy(o)        # initial operator 
    rCtvals = Vector{Float64}([])# vector for C(t) values storing
    iCtvals = Vector{Float64}([])

    # t = 0 
    WW = W*W
    expval0 = expectation_value(WW, ket)

    C0real = real(expval0)
    C0imag = imag(expval0)

    push!(rCtvals, C0real)
    push!(iCtvals, C0imag)

    # Number of sampling steps
    N_tau = length(sample_list)

    # Desired number of measurement points (coarse-grain the timeline for plotting)
    N_meas = min(tot_measurements, N_tau)   # do not request more measurements than steps

    # Build measurement indices evenly spaced between 1 and N_tau (inclusive)
    # Use LinRange then round to Int to avoid zero/step problems.
    meas_lst = unique(round.(Int, LinRange(1, N_tau, N_meas)))

    # Corresponding time grid:
    time_grid =  collect(range(0.0, stop=N_tau * dt/lamb, length=N_meas))

    sort!(meas_lst)             # ensure ascending order
    # optional: show how many actual measurement points we will have
    @printf("Will take %d measurements at indices (first 20): %s\n", length(meas_lst),
            string(meas_lst[1:min(end,20)]))

    # prepare counters
    meas_idx = 1
    next_meas = meas_lst[meas_idx]
    sample_idx = 1

    # Perform C(t) estimation using qDRIFT
    for idx in sample_list
        Pi = ops[idx]
        s = sign(coeffs[idx])              # sign(h_j)
        theta = 2 * τ * s                  #Angles -- keep if evolve! expects this convention
        pb = PauliBasis(Pi)
        evolve!(Wt, pb, theta)             # in-place evolve the operator Wt
        coeff_clip!(Wt, thresh=thresh)
        WWt = W * Wt                       # OTOC-like product

        # Check whether to measure at this step
        if sample_idx == next_meas
            # Measurements:
            expval = expectation_value(WWt, ket) # Contraction with reference ket
            Ctreal = real(expval)
            Ctimag = imag(expval)
            push!(rCtvals, Ctreal)
            push!(iCtvals, Ctimag)

            # advance to next measurement index (if any)
            meas_idx += 1
            if meas_idx <= length(meas_lst)
                next_meas = meas_lst[meas_idx]
            else
                # no more measurements; you can break if you like to save time
                # break
                next_meas = typemax(Int) # sentinel that will never be reached
            end
        end

        sample_idx += 1
    end

    # Print selection if requested
    if print_selection
        println("=== qDRIFT Selection (first 200 shown) ===")
        println(sample_list[1:min(end,200)])
        println("... (total samples = $Nsamples)")
    end

    # Plot distribution if requested
    if plot
        plt_bar = bar!(probs, legend=false, xlabel="Term index j", ylabel="p_j",
            title = @sprintf("qDRIFT sampling distribution (λ = %.6g)", λ))
        savefig(plt_bar, "qdrift_sampling_distribution.pdf")
    end

    return (sample=sample_list, probs=probs, time_grid=time_grid,
            λ=λ, τ=τ, N=Nsamples, RCt=rCtvals, ICt=iCtvals)
end
