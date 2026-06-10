#
# Here we track the L2 norm of the majorana weight for H2
#
using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Statistics
using Printf
using Random
using Plots
using Plots.PlotMeasures
using NPZ

# ============================================================
# Configuration
# ============================================================
struct PPConfig
    n_steps::Int
    dt::Float64
    coeff_threshold::Float64
    normalize_weights::Bool
end

"""
 Builds p_k(t) = sum_{P: majorana_weight(P) = k} |c_P(t)|², 
 the distribution of the majorana weight
"""
function weight_profile(W::PauliSum{N,T}; normalize::Bool=true) where {N,T}
    hist = zeros(Float64, 2*N + 1)
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

"""
 Entropy of the Majorana-weight distribution.
 How broadly the operator spreads accross Majorana weights.
 Not the entropy of the full propagated operator over
 individual Pauli strings.
"""
function weight_entropy(w::AbstractVector; ϵ::Float64=1e-15)
    p = normalized(w)
    p = p[p .> ϵ]
    return -sum(p .* log.(p))
end

function time_series_metrics(snapshots::Vector{<:AbstractVector})
    nt = length(snapshots)
    ent = zeros(Float64, nt)

    for t in 1:nt
        ent[t] = weight_entropy(snapshots[t])
    end

    return (
        entropy = ent,
    )
end

# ============================================================
# Pauli Propagation and Metrics
# ============================================================
"""
 Pauli Propagation, tracking the majorana weight distribution in this case,
 we do not track/compute the correlation function and in previous versions. 
  cfg :: is a structure containing the parameters for the Pauli Propagation.
  This new version needs to be tested to check performance (time).
"""
function weight_dynamics(O0::PauliSum{N,T}, H::PauliSum{N,T}, cfg::PPConfig) where {N,T}

    Ot = deepcopy(O0)
    snapshots = Vector{Vector{Float64}}()
    raw_norms = Float64[]

    w0, n0 = weight_profile(Ot; normalize=cfg.normalize_weights)
    push!(snapshots, w0)
    push!(raw_norms, n0)

    generators, coeffs = QuantumChemQC.gens_from_H(H)

    @printf("Total Pauli rotations: %d\n", length(coeffs) * cfg.n_steps)

    for step in 1:cfg.n_steps
        for j in eachindex(coeffs)
            θ = 2 * cfg.dt * coeffs[j]

            evolve!(Ot, generators[j], θ)
            QuantumChemQC.coeff_clip!(Ot; thresh=cfg.coeff_threshold)
        end

        w, norm2 = weight_profile(Ot; normalize=cfg.normalize_weights)
        push!(snapshots, w)
        push!(raw_norms, norm2)
    end

    tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

    return (
        tgrid = tgrid,
        snapshots = snapshots,
        raw_norms = raw_norms,
        metrics = time_series_metrics(snapshots),
    )
end

# ============================================================
# Running functions
# ============================================================
"""
 Run a single point in the PES 
 This function is intended to propagate a Z_i operator,
 whose i-th qubit acts on z_qubit.
"""
function run_single_distance(
    R::Float64,
    data_path::String,
    cfg::PPConfig;
    z_qubit::Int = 1,
    block::Bool = false,
    NOI::Bool = false,
    UNRESTRICTED::Bool = false,
)
    data = npzread(data_path)
    if UNRESTRICTED
        @printf("Using UNRESTRICTED data for R = %.6f\n", R)
        H1 = data["h1_a"]
    else
        H1 = data["h1e"]
    end
    Norbs = size(H1, 1)
    Nqubits = 2 * Norbs

    println()
    println("="^80)
    @printf("Running R = %.6f, Nqubits = %d\n", R, Nqubits)
    println("Path: ", data_path)

    if UNRESTRICTED
        println("Building UNRESTRICTED Hamiltonian")
        H = QuantumChemQC.molecular_hamiltonian_uhf(
            Norbs, 
            data_path,
            NOI=NOI,
            block=block)
    else
        println("Building RESTRICTED Hamiltonian")
        H = QuantumChemQC.molecular_hamiltonian(
        Norbs,
        data_path,
        NOI = NOI,
        block = block)
    end
    

    O = PauliSum(Pauli(Nqubits, Z = [z_qubit]))

    w_init, _ = weight_profile(O; normalize = true)
    println("Initial normalized Majorana-weight profile:")
    display(w_init)

    out = weight_dynamics(O, H, cfg)

    @printf("Initial norm² = %.12e\n", out.raw_norms[1])
    @printf("Final norm²   = %.12e\n", out.raw_norms[end])
    @printf("Relative norm² = %.12e\n", out.raw_norms[end] / out.raw_norms[1])
    @printf("Initial entropy = %.8f\n", out.metrics.entropy[1])
    @printf("Final entropy   = %.8f\n", out.metrics.entropy[end])

    return (
        R = R,
        data_path = data_path,
        Norbs = Norbs,
        Nqubits = Nqubits,
        out = out,
    )
end

"""
 Run a scan over multiple distances, given a list of (R, data_path) pairs.
"""
function run_distance_scan(cases, cfg::PPConfig;
             z_qubit::Int = 1, block::Bool = false, 
             NOI::Bool = false,
             UNRESTRICTED::Bool = false,)
    
             scan = Dict{Float64,Any}()

    for (R, path) in cases
        scan[R] = run_single_distance(
            R,
            path,
            cfg;
            z_qubit = z_qubit,
            block = block,
            NOI = NOI,
            UNRESTRICTED = UNRESTRICTED,
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

    p_entropy = plot(
        xlabel = "time",
        ylabel = "weight entropy",
        title = "Majorana-weight distribution entropy",
        lw = 2,
    )


    for R in Rs
        out = scan[R].out
        label = @sprintf("R=%.4g", R)

        plot!(
            p_entropy,
            out.tgrid,
            out.metrics.entropy;
            label = label,
            lw = 2,
        )
    end

    return plot(
        p_entropy,
        #layout = (2, 2),
        #size = (1100, 800),
    )
end

function plot_selected_heatmaps(scan; selected_Rs = [0.7414])
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
            titlelocation = :center,
            left_margin = 15mm,
            right_margin = 15mm,
            top_margin = 5mm,
            bottom_margin = 8mm,
        )

        push!(plots, p)
    end

    return plot(
        plots...;
        layout = (length(plots), 1),
        size = (850, 300 * length(plots)),
    )
end

function results_scan_table(scan)
    Rs = sorted_Rs(scan)

    println()
    println("R        final_S       max_S        norm_ret")
    println("-"^95)

    for R in Rs
        out = scan[R].out
        final_S = out.metrics.entropy[end]
        max_S = maximum(out.metrics.entropy)
        norm_ret = out.raw_norms[end] / out.raw_norms[1]

        @printf(
            "%-8.4f %-13.6f %-13.6f %-13.8f\n",
            R,
            final_S,
            max_S,
            norm_ret,
        )
    end
end

# ============================================================
# DATA INPUT AND PARAMETERS
# ============================================================
# D4 Dissociation model
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors" # Canonical RHF
const CASES = [
    (0.50, joinpath(BASE, "h4_RHF_D4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_RHF_D4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_RHF_D4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_RHF_D4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_RHF_D4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_RHF_D4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_RHF_D4_R_3_tensors.npz"))
]=#
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/UHF-tensors" # UHF
const CASES = [
    (0.50, joinpath(BASE, "h4_UHF_D4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_UHF_D4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_UHF_D4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_UHF_D4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_UHF_D4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_UHF_D4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_UHF_D4_R_3_tensors.npz"))
]=#
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/NO-CCSD_tensors" # CCSD natural orbitals
const CASES = [
    (0.50, joinpath(BASE, "h4_CCSDNO_D4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_CCSDNO_D4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_CCSDNO_D4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_CCSDNO_D4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_CCSDNO_D4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_CCSDNO_D4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_CCSDNO_D4_R_3_tensors.npz"))
]=#

# LIN4 Dissociation model
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors" # Canonical RHF
const CASES = [
    (0.50, joinpath(BASE, "h4_RHF_LIN4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_RHF_LIN4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_RHF_LIN4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_RHF_LIN4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_RHF_LIN4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_RHF_LIN4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_RHF_LIN4_R_3_tensors.npz"))
]=#
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/UHF-tensors" # UHF
const CASES = [
    (0.50, joinpath(BASE, "h4_UHF_LIN4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_UHF_LIN4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_UHF_LIN4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_UHF_LIN4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_UHF_LIN4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_UHF_LIN4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_UHF_LIN4_R_3_tensors.npz"))
]=#
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/NO-CCSD_tensors" # CCSD natural orbitals
const CASES = [
    (0.50, joinpath(BASE, "h4_CCSDNO_LIN4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_CCSDNO_LIN4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_CCSDNO_LIN4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_CCSDNO_LIN4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_CCSDNO_LIN4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_CCSDNO_LIN4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_CCSDNO_LIN4_R_3_tensors.npz"))
]=#

# S4 Dissociation model
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors" # Canonical RHF
const CASES = [
    (0.50, joinpath(BASE, "h4_RHF_S4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_RHF_S4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_RHF_S4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_RHF_S4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_RHF_S4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_RHF_S4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_RHF_S4_R_3_tensors.npz"))
]=#
#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/UHF-tensors" # UHF
const CASES = [
    (0.50, joinpath(BASE, "h4_UHF_S4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_UHF_S4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_UHF_S4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_UHF_S4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_UHF_S4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_UHF_S4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_UHF_S4_R_3_tensors.npz"))
]=#
const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/NO-CCSD_tensors" # CCSD natural orbitals
const CASES = [
    (0.50, joinpath(BASE, "h4_CCSDNO_S4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_CCSDNO_S4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_CCSDNO_S4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_CCSDNO_S4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_CCSDNO_S4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_CCSDNO_S4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_CCSDNO_S4_R_3_tensors.npz"))
]

cfg = PPConfig(
    500,       # n_steps
    0.1,     # dt
    1e-12,   # coeff_threshold
    false,   # normalize_weights
)

scan = run_distance_scan(
    CASES,
    cfg;
    z_qubit = 4,
    block = false,
    NOI = false,
    UNRESTRICTED = false,
)

p_summary = plot_distance_scan_summary(scan)
display(p_summary)
savefig(p_summary, "majorana_distance_scan_summary.png")

p_heat = plot_selected_heatmaps(scan; selected_Rs = [0.50, 0.7414, 1.50, 3.00])
display(p_heat)
savefig(p_heat, "majorana_weight_heatmaps_selected_R.png")

results_scan_table(scan)