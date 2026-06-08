#=
  We use Pauli Propagation to study the operator-size distribution of
  physically motivated fermionic operators. 
  In particular, we look at the Majorana-weight dynamics, where we track
  the L2 norm of the operators with a given Majorana weight.
  The growth, saturations and basis dependence of this distributions,
  are interpretes as quantum-information diagnostics of interaction
  induced operato complexity.
  This idea is based from the work from Joven and Bastidas:
  https://arxiv.org/abs/2405.12289
=#
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

# ============================================================
# Pauli Propagation and Metrics
# ============================================================
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

            QuantumChemQC.coeff_clip!(Ot; thresh=cfg.clip_rotation)
            before = expectation_value(O0 * Ot, ket)

            QuantumChemQC.coeff_clip!(Ot; thresh=cfg.clip_step)
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
# Load Hamiltonian and define operators, reference ket
# ============================================================
#data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors/h4_RHF_H2H2_R_0p7414_tensors.npz"
data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors/h4_RHF_H2H2_R_6p0_tensors.npz"
data = npzread(data_path)
H1 = data["h1e"]
Norbs = size(H1,1)    
println("Number of orbitals: $Norbs")
H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=false)

# - - - Define Reference ket
ket, _  = QuantumChemQC.string_to_ket("11110000") #Interleaved
println("Initial state Expectation Value:")
e1 = expectation_value(H,ket)
@printf("E(0) = %.6f\n", e1)

# - - - Define Initial Operator
#O = QuantumChemQC.homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
#println("HOMO-LUMO excitation operator:")
#display(O)
O = Pauli(2*Norbs, Z=[4])
O = PauliSum(O)
# ============================================================
# Run Pauli Propagation
# ============================================================
cfg = PPConfig(
    50,      # n_steps
    0.1,      # dt
    1e-12,    # clip_rotation
    1e-12,    # clip_step
    false,     # normalize_weights
)

out = heisenberg_weight_dynamics(ket, O, H, cfg)

# ============================================================
# Basic diagnostics
# ============================================================

println("Initial norm² = ", out.raw_norms[1])
println("Final norm²   = ", out.raw_norms[end])
println("Relative norm² = ", out.raw_norms[end] / out.raw_norms[1])

println("Initial mean Majorana weight = ", out.metrics.mean_weight[1])
println("Final mean Majorana weight   = ", out.metrics.mean_weight[end])

println("Initial entropy = ", out.metrics.entropy[1])
println("Final entropy   = ", out.metrics.entropy[end])

# ============================================================
# Plots
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
