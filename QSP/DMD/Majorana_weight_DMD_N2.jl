# PP + DMD
# Script designed to track the majorana weight distribution of 
# physically motivated operators evolved under the Heisenberg picture, 
# as a function of time.
#
# - Physics Layer: Hamiltonian, initial seed operators, ket/state.
# - Propagation layer: Hamiltonian evolution, Trotterization, Pauli propagation.
# - Diagnostics layer: Majorana-weight profiles, moments, entropy, norm loss.
# - DMD/Takens layer: DMD/Takens analysis of the snapshots, reconstruction, prediction.

using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Statistics
using Printf
using Random
using Plots
using NPZ

# ============================================================
# Configuration
# ============================================================

struct PPConfig
    n_steps::Int
    dt::Float64
    clip_rotation::Float64
    clip_step::Float64
    normalize_weights::Bool
end

struct DMDConfig
    q_values::Vector{Int}
    r::Union{Nothing,Int}
    train_fraction::Float64
end

struct SeedOperator{N,T}
    name::String
    op::PauliSum{N,T}
    sector::Symbol   # :odd, :even, :number_conserving, etc.
end

struct PPRunResult
    seed_name::String
    tgrid::Vector{Float64}
    snapshots::Vector{Vector{Float64}}
    raw_norms::Vector{Float64}
    corr_real::Vector{Float64}
    corr_imag::Vector{Float64}
    metrics::NamedTuple
end

weight(p::PauliBasis) = count_ones(p.x | p.z)

"""
Majorana weight of a Pauli string under the chosen Jordan-Wigner convention.

Important: this assumes the same qubit/orbital ordering used by the molecular
Hamiltonian construction.
"""
function weight_profile(W::PauliSum{N,T}; normalize::Bool=true) where {N,T}
    hist = zeros(Float64, 2N + 1)
    total = 0.0

    for (P, c) in W
        k = QuantumChemQC.majorana_weight(P)
        a2 = abs2(c)
        hist[k + 1] += a2
        total += a2
    end

    if normalize && total > 0
        hist ./= total
    end

    return hist, total
end

function normalized(w::AbstractVector)
    s = sum(w)
    s <= 0 && return copy(w)
    return w ./ s
end

function weight_stats(w::AbstractVector)
    p = normalized(w)
    ks = collect(0:length(p)-1)
    μ = sum(ks .* p)
    σ2 = sum(((ks .- μ).^2) .* p)
    return μ, σ2
end

function weight_entropy(w::AbstractVector; ϵ::Float64=1e-15)
    p = normalized(w)
    p = p[p .> ϵ]
    return -sum(p .* log.(p))
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

    return (
        mean_weight = μ,
        variance = σ2,
        entropy = ent,
    )
end

#
# Now snapshots are probability distributions over Majorana weight,
# while raw_norms separately tell whether clipping is changing the operator norm.

# =============================================================
# Propagation Layer 
# =============================================================
coeff_clip!(ps; thresh=1e-16) = filter!(p -> abs(p.second) > thresh, ps)

"""
 Evolution making clipping and norm tracking explicit.
"""
function heisenberg_weight_dynamics(
    ket,
    O0::PauliSum{N,T},
    H::PauliSum{N,T},
    cfg::PPConfig,
) where {N,T}

    Ot = deepcopy(O0)

    snapshots = Vector{Vector{Float64}}()
    raw_norms = Float64[]
    corr_real = Float64[]
    corr_imag = Float64[]

    w0, n0 = weight_profile(Ot; normalize=cfg.normalize_weights)
    push!(snapshots, w0)
    push!(raw_norms, n0)

    c0 = expectation_value(O0 * Ot, ket)
    push!(corr_real, real(c0))
    push!(corr_imag, imag(c0))

    generators, coeffs = QuantumChemQC.gens_from_H(H)

    @printf("Total Pauli rotations: %d\n", length(coeffs) * cfg.n_steps)

    for step in 1:cfg.n_steps
        accumulated_error = 0.0 + 0.0im

        for j in eachindex(coeffs)
            θ = 2 * cfg.dt * coeffs[j]

            evolve!(Ot, generators[j], θ)

            coeff_clip!(Ot; thresh=cfg.clip_rotation)
            before = expectation_value(O0 * Ot, ket)

            coeff_clip!(Ot; thresh=cfg.clip_step)
            after = expectation_value(O0 * Ot, ket)

            accumulated_error += after - before
        end

        w, norm2 = weight_profile(Ot; normalize=cfg.normalize_weights)
        push!(snapshots, w)
        push!(raw_norms, norm2)

        c = expectation_value(O0 * Ot, ket) + accumulated_error
        push!(corr_real, real(c))
        push!(corr_imag, imag(c))
    end

    tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

    return (
        tgrid = tgrid,
        snapshots = snapshots,
        raw_norms = raw_norms,
        corr_real = corr_real,
        corr_imag = corr_imag,
        metrics = time_series_metrics(snapshots),
    )
end

# ============================================================
# DMD/Takens Layer
# ============================================================
function train_test_dmd_error(
    snapshots::Vector{<:AbstractVector};
    q::Int,
    r::Union{Nothing,Int},
    train_fraction::Float64=0.7,
)
    X = QuantumChemQC.embed_snapshots(snapshots; q=q)
    n = size(X, 2)

    ntrain = clamp(floor(Int, train_fraction * n), q + 2, n - 1)

    Xtrain = X[:, 1:ntrain]
    Xtest  = X[:, ntrain:end]

    res = QuantumChemQC.fit_dmd(Xtrain; r=r)

    X1_test = Xtest[:, 1:end-1]
    X2_test = Xtest[:, 2:end]

    pred = res.A_ls * X1_test
    test_error = norm(X2_test - pred) / max(norm(X2_test), eps())

    return res, test_error
end

function dominant_rank(s::AbstractVector{<:Real}; energy::Float64=0.95)
    total = sum(abs2, s)
    total <= 0 && return 0

    acc = 0.0
    for (i, σ) in enumerate(s)
        acc += abs2(σ)
        if acc / total >= energy
            return i
        end
    end

    return length(s)
end

# ============================================================
# Load Hamiltonian and define operators, reference ket
# ============================================================
data_path = "/Users/admin/PycharmProjects/pyQCTools/DBF//N2_tensors_rhf_stable/1.1_tensors.npz" #Boys canonical localized
data = npzread(data_path)
H1 = data["h1e"]
Norbs = size(H1,1)    
H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=false)

# - - - Define Reference ket
#ket, _  = QuantumChemQC.string_to_ket("11111110001111111000") #Block
ket, _  = QuantumChemQC.string_to_ket("11111111111111000000") #Interleaved
println("Initial state Expectation Value:")
e1 = expectation_value(H,ket)
@printf("E(0) = %.6f\n", e1)

# - - - Define Excitation Operator
O = QuantumChemQC.homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
println("HOMO-LUMO excitation operator:")
display(O)

# ============================================================
# Run Pauli Propagation
# ============================================================
cfg = PPConfig(
    50,      # n_steps
    0.1,      # dt
    1e-12,    # clip_rotation
    1e-3,    # clip_step
    false,     # normalize_weights
)

out = heisenberg_weight_dynamics(ket, O, H, cfg)

# ============================================================
# 5. Basic diagnostics
# ============================================================

println("Initial norm² = ", out.raw_norms[1])
println("Final norm²   = ", out.raw_norms[end])
println("Relative norm² = ", out.raw_norms[end] / out.raw_norms[1])

println("Initial mean Majorana weight = ", out.metrics.mean_weight[1])
println("Final mean Majorana weight   = ", out.metrics.mean_weight[end])

println("Initial entropy = ", out.metrics.entropy[1])
println("Final entropy   = ", out.metrics.entropy[end])

# ============================================================
# 6. Plots
# ============================================================

W = hcat(out.snapshots...)

p1 = heatmap(
    out.tgrid,
    0:size(W, 1)-1,
    W,
    xlabel = "time",
    ylabel = "Majorana weight",
    title = "HOMO-LUMO Majorana-weight dynamics",
)

p2 = plot(
    out.tgrid,
    out.metrics.mean_weight,
    xlabel = "time",
    ylabel = "mean Majorana weight",
    title = "Mean Majorana weight",
    label = false,
    lw = 2,
)

p3 = plot(
    out.tgrid,
    out.metrics.entropy,
    xlabel = "time",
    ylabel = "entropy",
    title = "Majorana-weight entropy",
    label = false,
    lw = 2,
)

p4 = plot(
    out.tgrid,
    out.raw_norms ./ out.raw_norms[1],
    xlabel = "time",
    ylabel = "relative norm²",
    title = "Coefficient norm retention",
    label = false,
    lw = 2,
)

display(p1)
display(p2)
display(p3)
display(p4)