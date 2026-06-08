#=
  Simlar as weight_dynamics.jl, but here we test for the full PES
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
# Distance scan over H2-H2 separations
# ============================================================

const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors"

const CASES = [
    (0.7414, joinpath(BASE, "h4_RHF_H2H2_R_0p7414_tensors.npz")),
    (1.0,    joinpath(BASE, "h4_RHF_H2H2_R_1p0_tensors.npz")),
    (1.5,    joinpath(BASE, "h4_RHF_H2H2_R_1p5_tensors.npz")),
    (2.0,    joinpath(BASE, "h4_RHF_H2H2_R_2p0_tensors.npz")),
    (2.5,    joinpath(BASE, "h4_RHF_H2H2_R_2p5_tensors.npz")),
    (3.0,    joinpath(BASE, "h4_RHF_H2H2_R_3p0_tensors.npz")),
    (4.0,    joinpath(BASE, "h4_RHF_H2H2_R_4p0_tensors.npz")),
    (5.0,    joinpath(BASE, "h4_RHF_H2H2_R_5p0_tensors.npz")),
    (6.0,    joinpath(BASE, "h4_RHF_H2H2_R_6p0_tensors.npz")),
]


function run_single_distance(
    R::Float64,
    data_path::String,
    cfg::PPConfig;
    ket_string::String = "11110000",
    z_qubit::Int = 4,
    block::Bool = false,
    NOI::Bool = false,
)
    data = npzread(data_path)
    H1 = data["h1e"]
    Norbs = size(H1, 1)
    Nqubits = 2 * Norbs

    println()
    println("="^80)
    @printf("Running R = %.6f, Nqubits = %d\n", R, Nqubits)
    println("Path: ", data_path)

    H = QuantumChemQC.molecular_hamiltonian(
        Norbs,
        data_path,
        NOI = NOI,
        block = block,
    )

    ket, _ = QuantumChemQC.string_to_ket(ket_string)

    E0 = expectation_value(H, ket)
    @printf("Reference energy <ket|H|ket> = %.10f\n", real(E0))


    O = PauliSum(Pauli(Nqubits, Z = [z_qubit]))

    w_init, n_init = weight_profile(O; normalize = true)
    println("Initial normalized Majorana-weight profile:")
    display(w_init)

    out = heisenberg_weight_dynamics(ket, O, H, cfg)

    @printf("Initial norm² = %.12e\n", out.raw_norms[1])
    @printf("Final norm²   = %.12e\n", out.raw_norms[end])
    @printf("Relative norm² = %.12e\n", out.raw_norms[end] / out.raw_norms[1])
    @printf("Initial mean weight = %.8f\n", out.metrics.mean_weight[1])
    @printf("Final mean weight   = %.8f\n", out.metrics.mean_weight[end])
    @printf("Initial entropy = %.8f\n", out.metrics.entropy[1])
    @printf("Final entropy   = %.8f\n", out.metrics.entropy[end])

    return (
        R = R,
        data_path = data_path,
        Norbs = Norbs,
        Nqubits = Nqubits,
        E0 = real(E0),
        out = out,
    )
end


function run_distance_scan(
    cases,
    cfg::PPConfig;
    ket_string::String = "11110000",
    z_qubit::Int = 4,
    block::Bool = false,
    NOI::Bool = false,
)
    scan = Dict{Float64,Any}()

    for (R, path) in cases
        scan[R] = run_single_distance(
            R,
            path,
            cfg;
            ket_string = ket_string,
            z_qubit = z_qubit,
            block = block,
            NOI = NOI,
        )
    end

    return scan
end

# ============================================================
# Comparative plots across H2-H2 separations
# ============================================================

function sorted_Rs(scan)
    return sort(collect(keys(scan)))
end


function normalized_snapshots_matrix(out)
    W = hcat(out.snapshots...)

    for j in axes(W, 2)
        s = sum(W[:, j])
        if s > 0
            W[:, j] ./= s
        end
    end

    return W
end


function plot_distance_scan_summary(scan)
    Rs = sorted_Rs(scan)

    p_mean = plot(
        xlabel = "time",
        ylabel = "mean Majorana weight",
        title = "Mean Majorana weight vs H2-H2 separation",
        lw = 2,
    )

    p_entropy = plot(
        xlabel = "time",
        ylabel = "weight entropy",
        title = "Majorana-weight entropy",
        lw = 2,
    )

    p_norm = plot(
        xlabel = "time",
        ylabel = "relative norm²",
        title = "Coefficient norm retention",
        lw = 2,
        yformatter = y -> @sprintf("%.5f", y) # Shows 2 decimal places
    )

    for R in Rs
        out = scan[R].out
        label = @sprintf("R=%.4g", R)

        plot!(
            p_mean,
            out.tgrid,
            out.metrics.mean_weight;
            label = label,
            lw = 2,
        )

        plot!(
            p_entropy,
            out.tgrid,
            out.metrics.entropy;
            label = label,
            lw = 2,
        )

        plot!(
            p_norm,
            out.tgrid,
            out.raw_norms ./ out.raw_norms[1];
            label = label,
            lw = 2,
        )
    end

    final_mean = [scan[R].out.metrics.mean_weight[end] for R in Rs]
    final_entropy = [scan[R].out.metrics.entropy[end] for R in Rs]

    p_final = plot(
        Rs,
        final_mean;
        marker = :circle,
        xlabel = "H2-H2 separation R",
        ylabel = "final mean Majorana weight",
        title = "Final operator size vs separation",
        label = "mean weight",
        lw = 2,
    )

    plot!(
        p_final,
        Rs,
        final_entropy;
        marker = :square,
        ylabel = "final diagnostic value",
        label = "entropy",
        lw = 2,
    )

    return plot(
        p_mean,
        p_entropy,
        p_norm,
        p_final;
        layout = (2, 2),
        size = (1100, 800),
    )
end

function plot_selected_heatmaps(scan; selected_Rs = [0.7414, 2.0, 6.0])
    plots = []

    for R in selected_Rs
        out = scan[R].out
        W = normalized_snapshots_matrix(out)

        p = heatmap(
            out.tgrid,
            0:size(W, 1)-1,
            W;
            xlabel = "time",
            ylabel = "Majorana weight",
            title = @sprintf("R = %.4g", R),
        )

        push!(plots, p)
    end

    return plot(
        plots...;
        layout = (length(plots), 1),
        size = (850, 300 * length(plots)),
    )
end

function distance_scan_table(scan)
    Rs = sorted_Rs(scan)

    println()
    println("R        E0              final_mu      max_mu        final_S       max_S        norm_ret")
    println("-"^95)

    for R in Rs
        out = scan[R].out

        final_mu = out.metrics.mean_weight[end]
        max_mu = maximum(out.metrics.mean_weight)
        final_S = out.metrics.entropy[end]
        max_S = maximum(out.metrics.entropy)
        norm_ret = out.raw_norms[end] / out.raw_norms[1]

        @printf(
            "%-8.4f %-15.8f %-13.6f %-13.6f %-13.6f %-13.6f %-13.8f\n",
            R,
            scan[R].E0,
            final_mu,
            max_mu,
            final_S,
            max_S,
            norm_ret,
        )
    end
end

cfg = PPConfig(
    50,       # n_steps
    0.1,     # dt
    1e-12,   # clip_rotation
    1e-12,   # clip_step
    false,   # normalize_weights
)

scan = run_distance_scan(
    CASES,
    cfg;
    ket_string = "11110000",
    z_qubit = 4,
    block = false,
    NOI = false,
)


p_summary = plot_distance_scan_summary(scan)
display(p_summary)
savefig(p_summary, "majorana_distance_scan_summary.png")

p_heat = plot_selected_heatmaps(scan; selected_Rs = [0.7414, 2.0, 6.0])
display(p_heat)
savefig(p_heat, "majorana_weight_heatmaps_selected_R.png")

distance_scan_table(scan)