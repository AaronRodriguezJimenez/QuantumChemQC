using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Statistics
using Printf
using Random
using Plots
using NPZ

coeff_clip!(ps; thresh=1e-16) = filter!(p -> abs(p.second) > thresh, ps)

function string_to_ket(bits::AbstractString)
    idx = 0
    for (i, ch) in enumerate(bits)
        if ch == '1'
            idx += 1 << (i - 1)
        end
    end
    return Ket{length(bits)}(idx), idx
end

# ------------------------------------------------------------
# Pauli propagation
# ------------------------------------------------------------

function evolve!(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N, T}
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

function extract_hamiltonian_coeffs_and_ops(H::PauliSum{N, T}) where {N, T}
    ops = PauliBasis{N}[]
    coeffs = Float64[]
    for (p, c) in H
        push!(ops, p)
        push!(coeffs, float(c))
    end
    return ops, coeffs
end

# ------------------------------------------------------------
# Weight diagnostics
# ------------------------------------------------------------

weight(p::PauliBasis) = count_ones(p.x | p.z)

function weight_profile(W::PauliSum{N, T}) where {N, T}
    hist = zeros(Float64, N + 1)
    total = 0.0

    for (P, c) in W
        k = weight(P)
        a2 = abs2(c)
        hist[k + 1] += a2
        total += a2
    end

    if total > 0
        hist ./= total
    end

    return hist
end

function weight_stats(w::AbstractVector)
    ks = 0:length(w)-1
    μ = sum(ks .* w)
    σ2 = sum(((ks .- μ).^2) .* w)
    return μ, σ2
end

# ------------------------------------------------------------
# DMD
# ------------------------------------------------------------

struct DMDResult
    A_ls::Matrix{Float64}
    A_tilde::Matrix{ComplexF64}
    modes::Matrix{ComplexF64}
    evals::Vector{ComplexF64}
    amplitudes::Vector{ComplexF64}
    singular_values::Vector{Float64}
    residual_rel::Float64
end

function fit_dmd(X::AbstractMatrix{<:Real}; r::Union{Nothing,Int}=nothing, tol=1e-10)
    X1 = X[:, 1:end-1]
    X2 = X[:, 2:end]

    A_ls = Matrix(X2 * pinv(X1))
    residual_rel = norm(X2 - A_ls * X1) / max(norm(X2), eps())

    F = svd(X1; full=false)
    U, s, V = F.U, F.S, F.V

    if r === nothing
        r = max(count(>(tol * s[1]), s), 1)
    else
        r = min(r, length(s))
    end

    Ur = U[:, 1:r]
    sr = s[1:r]
    Vr = V[:, 1:r]

    Sinv = Diagonal(1.0 ./ sr)
    A_tilde = Matrix(Ur' * X2 * Vr * Sinv)

    eig = eigen(A_tilde)
    λ = eig.values
    W = eig.vectors

    Φ = Matrix(X2 * Vr * Sinv * W)
    b = Φ \ complex.(X[:, 1])

    return DMDResult(A_ls, A_tilde, Φ, λ, b, s, residual_rel)
end

function print_dmd_summary(res::DMDResult; dt::Real=1.0, topk::Int=5)
    λ = res.evals
    growth = log.(abs.(λ)) ./ dt
    freq = angle.(λ) ./ dt
    idx = sortperm(abs.(λ), rev=true)

    println("---- DMD summary ----")
    println("relative LS residual = $(res.residual_rel)")
    println("numerical rank       = $(size(res.A_tilde, 1))")
    println("top singular values   = ", res.singular_values[1:min(end, topk)])

    println("\nDominant modes:")
    for j in 1:min(topk, length(idx))
        i = idx[j]
        println("  λ[$i] = $(λ[i])")
        println("      |λ| = $(abs(λ[i]))")
        println("      growth rate = $(growth[i])")
        println("      frequency    = $(freq[i])")
        println("      amplitude    = $(res.amplitudes[i])")
    end
end

function delay_embed(snapshots::Vector{<:AbstractVector}, q::Int)
    m = length(snapshots)
    d = length(snapshots[1])
    ncols = m - q + 1
    X = zeros(Float64, d * q, ncols)

    for k in 1:ncols
        for j in 1:q
            X[(j-1)*d + 1 : j*d, k] .= snapshots[k + j - 1]
        end
    end

    return X
end

function dmd_bootstrap(snapshots::Vector{<:AbstractVector}; B::Int=100, r::Union{Nothing,Int}=nothing, seed::Int=1)
    rng = MersenneTwister(seed)
    m = length(snapshots)
    evals = Vector{Vector{ComplexF64}}(undef, B)

    for b in 1:B
        lo = rand(rng, 1:m-3)
        hi = rand(rng, lo+2:m)
        Xb = hcat(snapshots[lo:hi]...)
        evals[b] = fit_dmd(Xb; r=r).evals
    end

    return evals
end

# ------------------------------------------------------------
# Main evolution loop
# ------------------------------------------------------------

function evolution_op(ket, o::PauliSum{N, T}, H::PauliSum{N, T}, n_intervals, dt;
    thresh::Float64=1e-3) where {N, T}

    O0 = deepcopy(o)
    Ot = deepcopy(o)

    corr_real = Float64[]
    corr_imag = Float64[]
    snapshots = Vector{Vector{Float64}}()

    push!(snapshots, weight_profile(Ot))

    c0 = expectation_value(O0 * Ot, ket)
    push!(corr_real, real(c0))
    push!(corr_imag, imag(c0))

    generators, angles = extract_hamiltonian_coeffs_and_ops(H)
    nt = length(angles)

    println("Total Rotations: ", nt * n_intervals)

    for _ in 1:n_intervals
        accumulated_error = 0.0 + 0.0im

        for j in 1:nt
            θ = 2 * dt * angles[j]
            evolve!(Ot, generators[j], θ)

            coeff_clip!(Ot; thresh=1e-12)
            before = expectation_value(O0 * Ot, ket)

            coeff_clip!(Ot; thresh=thresh)
            after = expectation_value(O0 * Ot, ket)

            accumulated_error += after - before
        end

        push!(snapshots, weight_profile(Ot))

        c = expectation_value(O0 * Ot, ket) + accumulated_error
        push!(corr_real, real(c))
        push!(corr_imag, imag(c))
    end

    tgrid = collect(range(0.0, stop=n_intervals * dt, length=length(corr_real)))
    return corr_real, corr_imag, tgrid, snapshots
end

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

function plot_weight_heatmap(w_snapshots; tgrid=nothing)
    W = hcat(w_snapshots...)
    nweights, nt = size(W)

    if tgrid === nothing
        tgrid = 0:nt-1
    end

    heatmap(
        tgrid,
        0:nweights-1,
        W,
        xlabel = "time",
        ylabel = "Pauli weight",
        title = "Pauli Weight Dynamics",
        #legend = false
    )
end

function plot_weight_stack(w_snapshots; tgrid=nothing)
    W = hcat(w_snapshots...)
    nweights, nt = size(W)

    if tgrid === nothing
        tgrid = 0:nt-1
    end

    plt = plot(
        xlabel = "time",
        ylabel = "weight fraction",
        title = "Pauli Weight Distribution",
        legend = :right
    )

    for k in 1:nweights
        plot!(plt, tgrid, W[k, :], label = "k=$(k-1)", lw=2)
    end

    return plt
end

function plot_dmd_modes_abs(res::DMDResult; nmodes=4)
    r = min(nmodes, size(res.modes, 2))
    ks = 0:(size(res.modes, 1)-1)

    p = plot(
        xlabel = "Pauli weight",
        ylabel = "|mode amplitude|",
        title = "Magnitude of DMD modes",
        legend = :right
    )

    for j in 1:r
        plot!(p, ks, abs.(res.modes[:, j]), label = "mode $j", lw=2)
    end

    return p
end

# ------------------------------------------------------------
# Basic utilities
# ------------------------------------------------------------

mean_weight(w::AbstractVector) = sum((0:length(w)-1) .* w)

function weight_entropy(w::AbstractVector; ϵ::Float64=1e-15)
    p = clamp.(w, ϵ, 1.0)
    return -sum(p .* log.(p))
end

function dominant_rank(s::AbstractVector{<:Real}; energy::Float64=0.95)
    tot = sum(abs2, s)
    tot <= 0 && return 0
    acc = 0.0
    for (i, σ) in enumerate(s)
        acc += abs2(σ)
        if acc / tot >= energy
            return i
        end
    end
    return length(s)
end

function l2_relative(a::AbstractVector, b::AbstractVector)
    return norm(a .- b) / max(norm(b), eps())
end

function time_series_metrics(snapshots::Vector{<:AbstractVector})
    nt = length(snapshots)
    μ = zeros(Float64, nt)
    σ2 = zeros(Float64, nt)
    ent = zeros(Float64, nt)

    for t in 1:nt
        μ[t], σ2[t] = weight_stats(snapshots[t])
        ent[t] = weight_entropy(snapshots[t])
    end

    return (mean_weight = μ, variance = σ2, entropy = ent)
end

# ------------------------------------------------------------
# Delay embedding
# ------------------------------------------------------------

"""
Takens-style delay embedding.
Each column is [x_t; x_{t+1}; ...; x_{t+q-1}].
"""
function delay_embed(snapshots::Vector{<:AbstractVector}, q::Int)
    m = length(snapshots)
    d = length(snapshots[1])
    @assert q >= 1 "q must be at least 1"
    @assert m >= q "need at least q snapshots"

    ncols = m - q + 1
    X = zeros(Float64, d * q, ncols)

    for k in 1:ncols
        for j in 1:q
            X[(j-1)*d + 1 : j*d, k] .= snapshots[k + j - 1]
        end
    end

    return X
end

embed_snapshots(snapshots::Vector{<:AbstractVector}; q::Int=1) = q == 1 ? hcat(snapshots...) : delay_embed(snapshots, q)

function embed_fit_summary(snapshots::Vector{<:AbstractVector}; q::Int=1, r::Union{Nothing,Int}=2)
    X = embed_snapshots(snapshots; q=q)
    res = fit_dmd(X; r=r)
    X1 = X[:, 1:end-1]
    X2 = X[:, 2:end]
    rel_fit = norm(X2 - res.A_ls * X1) / max(norm(X2), eps())
    rank95 = dominant_rank(res.singular_values; energy=0.95)
    rank99 = dominant_rank(res.singular_values; energy=0.99)
    return (
        X = X,
        res = res,
        rel_fit = rel_fit,
        rank95 = rank95,
        rank99 = rank99,
    )
end

# ------------------------------------------------------------
# Signal generators
# ------------------------------------------------------------

function synthetic_signal(m::Int; ω1::Float64=0.12, ω2::Float64=0.31)
    t = collect(0:m-1)
    x1 = sin.(ω1 .* t) .+ 0.3 .* cos.(ω2 .* t)
    x2 = cos.(ω1 .* t) .- 0.2 .* sin.(ω2 .* t)
    x3 = 0.5 .* sin.(0.5 .* ω1 .* t .+ 0.1)
    return [Float64[x1[i], x2[i], x3[i]] for i in 1:m]
end

function white_noise_signal(m::Int, d::Int; seed::Int=1)
    rng = MersenneTwister(seed)
    return [randn(rng, d) for _ in 1:m]
end

"""
 Correlated random signal (AR process)

Pure white noise is almost too easy. 
A more realistic null model is correlated stochastic dynamics
AR(1) formula: x_t = α * x_{t-1} + noise
"""
function ar1_signal(m::Int, d::Int; α::Float64=0.9, seed::Int=1)
    rng = MersenneTwister(seed)
    x = zeros(Float64, d)
    out = Vector{Vector{Float64}}(undef, m)
    for t in 1:m
        x = α .* x .+ 0.1 .* randn(rng, d)
        out[t] = copy(x)
    end
    return out
end

"""
 Take PP snapshots and shuffle time order
"""
function shuffled_snapshots(snapshots::Vector{<:AbstractVector}; seed::Int=1)
    rng = MersenneTwister(seed)
    idx = randperm(rng, length(snapshots))
    return snapshots[idx]
end

# ------------------------------------------------------------
# Benchmark container
# ------------------------------------------------------------

struct BenchmarkResult
    name::String
    q_values::Vector{Int}
    fit_error::Vector{Float64}
    rank95::Vector{Int}
    rank99::Vector{Int}
    sv_first::Vector{Float64}
    sv_second::Vector{Float64}
end

function benchmark_signal(name::String, snapshots::Vector{<:AbstractVector}; dt = 0.1, q_values=1:10, r::Union{Nothing,Int}=2)
    qv = collect(q_values)
    fit_error = Float64[]
    rank95 = Int[]
    rank99 = Int[]
    sv_first = Float64[]
    sv_second = Float64[]
    @printf("\nBenchmarking %s...\n", name)
    for q in qv
        if length(snapshots) < q + 1
            @printf("Skipping q=%d for %s (not enough snapshots)\n", q, name)
            push!(fit_error, NaN)
            push!(rank95, 0)
            push!(rank99, 0)
            push!(sv_first, NaN)
            push!(sv_second, NaN)
            continue
        end

        out = embed_fit_summary(snapshots; q=q, r=r)
        print_dmd_summary(out.res; dt=dt, topk=5) # Summary for dominant topk modes
        
        push!(fit_error, out.rel_fit)
        push!(rank95, out.rank95)
        push!(rank99, out.rank99)
        push!(sv_first, isempty(out.res.singular_values) ? NaN : out.res.singular_values[1])
        push!(sv_second, length(out.res.singular_values) >= 2 ? out.res.singular_values[2] : NaN)
    end

    return BenchmarkResult(name, qv, fit_error, rank95, rank99, sv_first, sv_second)
end

# ------------------------------------------------------------
# PP-specific wrappers
# ------------------------------------------------------------

function run_pp_snapshots(ket, o, H, n_intervals::Int, dt::Real; threshold::Float64=1e-10)
    rRES, iRES, tgrid, w_snapshots = evolution_op(ket, o, H, n_intervals, dt; thresh=threshold)
    return (
        rRES = rRES,
        iRES = iRES,
        tgrid = tgrid,
        w_snapshots = w_snapshots,
        metrics = time_series_metrics(w_snapshots),
    )
end

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

function plot_benchmark_metric(results::Vector{BenchmarkResult}; metric::Symbol=:fit_error)
    p = plot(
        xlabel = "embedding q",
        ylabel = string(metric),
        title = "Takens + DMD benchmark",
        legend = :best,
        xscale = :linear,
    )

    for r in results
        vals = getproperty(r, metric)
        plot!(p, r.q_values, vals, marker=:circle, lw=2, label=r.name)
    end

    return p
end

function plot_benchmark_ranks(results::Vector{BenchmarkResult})
    p1 = plot(
        xlabel = "embedding q",
        ylabel = "rank for 95% energy",
        title = "Effective rank (95%)",
        legend = :best,
    )
    p2 = plot(
        xlabel = "embedding q",
        ylabel = "rank for 99% energy",
        title = "Effective rank (99%)",
        legend = :best,
    )

    for r in results
        plot!(p1, r.q_values, r.rank95, marker=:circle, lw=2, label=r.name)
        plot!(p2, r.q_values, r.rank99, marker=:circle, lw=2, label=r.name)
    end

    return p1, p2
end

function plot_singular_values_vs_q(results::Vector{BenchmarkResult})
    p = plot(
        xlabel = "embedding q",
        ylabel = "leading singular value",
        title = "Leading singular values vs embedding dimension",
        legend = :best,
    )

    for r in results
        plot!(p, r.q_values, r.sv_first, marker=:circle, lw=2, label="$(r.name): σ₁")
        plot!(p, r.q_values, r.sv_second, marker=:square, lw=2, label="$(r.name): σ₂")
    end

    return p
end

"""
 Build effective Heisenberg Hamiltonian from a molecule.
 This Hamiltonian is constructed as a matrix H_eff, 
 containing  the exchange coupling localized spins
 H = sum_{i<j} J_ij S_i * S_j
 with J = Heff_final[i, j] / np.sqrt(Si * Sj)
 For spin 1/2 systems S_i * S_j = 1/4 * (XiXj + YiYj + ZiZj)

 Effective Hamiltonian constructed following:
 Phys. Chem. Lett. 2015, 6, 1982−1988 

 NOTE: J matrix is in cm^-1 units.

 If include_diagonal, we include the onsite energies as Z
 terms. Usually false for pure Heisenberg models.
"""
function build_H_eff(data_path)

    data = npzread(data_path)

    Heff = data["Heff"] # Effective Hamiltonian in local ground state basis
    J = Matrix{Float64}(data["Jmatrix"])
    N = size(J, 1)

    H = PauliSum(N, Float64)

    #Pure Heisenberg exchange Xterms
    for i in 1:N-1
        for j in i+1:N
            Jij = J[i, j]
            coeff = Jij / 4.0
            
            # Pauli terms
            H += coeff * Pauli(N, X=[i,j])
            H += coeff * Pauli(N, Y=[i,j])
            H += coeff * Pauli(N, Z=[i,j])
        end
    end
   
    @printf("J matrix:")
    display(J)

    @printf("Reconstructed Hamiltonian:")
    display(H)

    coeff_clip!(H)

    return H, J, Heff   
end

# ------------------------------------------------------------
# Main benchmark
# ------------------------------------------------------------

function main_benchmark()
    println("========================================")
    println("Takens + DMD comparison benchmark")
    println("========================================")

    q_values = 1:12
    r = 4

    # Structured synthetic control
    structured = synthetic_signal(250)

    # Null controls
    white = white_noise_signal(250, 3; seed=2)
    ar1 = ar1_signal(250, 3; α=0.9, seed=3)

    # Shuffled version of the structured signal
    shuffled = shuffled_snapshots(structured; seed=4)

    # PP data
    # Get Molecular Hamiltonian
    script_dir = dirname(@__FILE__) 
    data_path = joinpath(script_dir, "ni_cubane_excoups.npz")
    H, J_matrix, Heff = build_H_eff(data_path)
    
    O = Pauli(4, X=[1,3])
    O = PauliSum(O)
    ket, _ = QuantumChemQC.string_to_ket("0000")
    
    time_step = 0.5
    pp = run_pp_snapshots(ket, O, H, 250, time_step; threshold=1e-10)
    shuffled_pp = shuffled = shuffled_snapshots(pp.w_snapshots; seed=4)

    # Benchmark all cases
    cases = BenchmarkResult[]
    push!(cases, benchmark_signal("synthetic", structured; dt=1.0, q_values=q_values, r=r))
    push!(cases, benchmark_signal("white noise", white; dt=1.0, q_values=q_values, r=r))
    #push!(cases, benchmark_signal("AR(1)", ar1; dt=1.0, q_values=q_values, r=r))
    #push!(cases, benchmark_signal("shuffled synthetic", shuffled; dt=1.0, q_values=q_values, r=r))
    push!(cases, benchmark_signal("PP weight snapshots", pp.w_snapshots; dt=time_step, q_values=1:min(12, length(pp.w_snapshots)-1), r=r))
    push!(cases, benchmark_signal("PP shuffled snapshots",shuffled_pp; dt=time_step, q_values=1:min(12, length(pp.w_snapshots)-1), r=r))

    println("\n--- Summary ---")
    for c in cases
        println("Case: ", c.name)
        println("  q values    : ", c.q_values)
        println("  fit error   : ", c.fit_error)
        println("  rank95      : ", c.rank95)
        println("  rank99      : ", c.rank99)
    end

    # Plots
    p_fit = plot_benchmark_metric(cases; metric=:fit_error)
    p_rank95, p_rank99 = plot_benchmark_ranks(cases)
    p_sv = plot_singular_values_vs_q(cases)

    display(p_fit)
    display(p_rank95)
    display(p_rank99)
    display(p_sv)

    # PP-specific diagnostics
    display(plot_weight_heatmap(pp.w_snapshots; tgrid=pp.tgrid))

    return (; cases, pp)
end

function signals()
    println("========================================")
    println("PLOT SIGNALS")
    println("========================================")
    
    m = 250
    labels = ["x1" "x2" "x3"] # Horizontal matrix for series labeling

    # 1. Generate Signals (Vector of 3-element Vectors)
    structured_raw = synthetic_signal(m)
    white_raw      = white_noise_signal(m, 3; seed=2)
    ar1_raw        = ar1_signal(m, 3; α=0.9, seed=3)
    shuffled_raw   = shuffled_snapshots(structured_raw; seed=4)
    # PP data
    # Get Effective Hamiltonian (computed with pyQCTools)
    # Define data file path
    script_dir = dirname(@__FILE__) 
    data_path = joinpath(script_dir, "ni_cubane_excoups.npz")
    H, J_matrix, Heff = build_H_eff(data_path)
    
    
    O = Pauli(4, X=[1,3])
    O = PauliSum(O)
    ket, _ = QuantumChemQC.string_to_ket("0000")

    pp = run_pp_snapshots(ket, O, H, 250, 0.5; threshold=1e-10)

    # 2. Conversion helper 
    # Transforms Vector{Vector{Float64}} into a Matrix{Float64} of size (m, 3)
    prepare(s) = reduce(hcat, s)'

    # 3. Plotting
    # We plot all 3 components on the same subplot for each signal type
    p1 = plot(prepare(structured_raw), title="Structured Signal", label=labels)
    p2 = plot(prepare(white_raw),      title="White Noise",       label=false)
    p3 = plot(prepare(ar1_raw),        title="AR(1) Process",     label=false)
    p4 = plot(prepare(shuffled_raw),   title="Shuffled Snapshots", label=false)
    p5 = plot(prepare(pp.w_snapshots), title="PP Weight Snapshots", label=false)

    # Combine into a 4-row stack
    combined_plot = plot(p1, p2, p3, p4, p5,
        layout=(5, 1), 
        size=(900, 1300), 
        link=:x,           # Links the x-axis for easier scrolling/comparison
        margin=5Plots.mm,
        ylabel="Value")

    display(combined_plot)
    return combined_plot

end

# run 
results = main_benchmark()
signals()

#=
E_orig = eigen(Heff)
V_orig = E_orig.vectors
Λ_orig = E_orig.values
println("Original Heff Eigenvalues:")
display(Λ_orig)

Hmat = Matrix(H)
# Siagonalization
E = eigen(Hmat)
V = E.vectors
Λ = E.values

println("Reconstructed Eigenvalues:")
display(Λ)
=#   