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

"""
Entropy and effective support of the full propagated operator

W(t) = sum_P c_P(t) P

This is different from the entropy of the Majorana-weight histogram.
It measures spreading over individual Pauli strings.
"""
function operator_entropy(W::PauliSum; ϵ::Float64 = 1e-15)
coeff_weights = Float64[]
total = 0.0

for (P, c) in W
a2 = abs2(c)
if a2 > 0
push!(coeff_weights, a2)
total += a2
end
end

if total <= 0
return (
entropy = 0.0,
neff_shannon = 0.0,
neff_ipr = 0.0,
total = 0.0,
)
end

p = coeff_weights ./ total
p = p[p .> ϵ]

S = -sum(p .* log.(p))
neff_shannon = exp(S)
neff_ipr = 1.0 / sum(p .^ 2)

return (
entropy = S,   #Full Pauli-string entropy S_op
neff_shannon = neff_shannon, #Effective number of active Pauli strings
neff_ipr = neff_ipr,  #Inverse participation ratio
total = total,
)
end

"""
Moments and entropy of the Majorana-weight distribution.

Input can be normalized or unnormalized.
The function normalizes internally.
"""
function weight_moments(hist::AbstractVector; ϵ::Float64 = 1e-15)
    p = normalized(hist)
    K = length(p)
    ks = collect(0:K-1)

    mean_w = sum(ks .* p)
    var_w = sum(((ks .- mean_w) .^ 2) .* p)

    p_nonzero = p[p .> ϵ]
    S_weight = -sum(p_nonzero .* log.(p_nonzero))
    neff_weight = exp(S_weight)

    return (
    mean = mean_w,
    variance = var_w,
    std = sqrt(var_w),
    entropy = S_weight,
    neff = neff_weight,
    )
end

"""
Total probability in Majorana weights k >= kmin.

For H2, useful choices are kmin = 4 and kmin = 6.
"""
function high_weight_leakage(hist::AbstractVector, kmin::Int)
    p = normalized(hist)

    if kmin < 0
        error("kmin must be nonnegative")
    end

    if kmin >= length(p)
        return 0.0
    end

    return sum(p[(kmin + 1):end])
end

# ============================================================
# Weight Index
# ============================================================
function trapz(times::AbstractVector, values::AbstractVector)
    length(times) == length(values) || error("times and values must have the same length")
    length(times) < 2 && return 0.0

    acc = 0.0

    for i in 1:(length(times)-1)
        dt = times[i+1] - times[i]
        acc += 0.5 * dt * (values[i+1] + values[i])
    end

    return acc
end


"""
Time-averaged growth index:

    Γ_X = (1/T) ∫ [X(t) - X(0)] dt

If normalize=true, divide by the maximum absolute deviation.
"""
function growth_index(
    times::AbstractVector,
    values::AbstractVector;
    subtract_initial::Bool = true,
    normalize::Bool = false,
)
    length(times) == length(values) || error("times and values must have the same length")

    T = times[end] - times[1]
    T <= 0 && error("time interval must be positive")

    x = collect(values)

    if subtract_initial
        x .-= x[1]
    end

    if normalize
        scale = maximum(abs.(x))
        if scale > 0
            x ./= scale
        end
    end

    return trapz(times, x) / T
end


"""
Maximum growth over the time window:

    Γ_X^max = max_t [X(t) - X(0)]
"""
function max_growth(values::AbstractVector)
    x0 = values[1]
    return maximum(values .- x0)
end

"""
Compute all growth indices for the various metrics in the output of weight_dynamics.
"""
function compute_gamma_indices(out; normalize::Bool = false)
    t = out.tgrid
    m = out.metrics

    return (
        Γ_operator_entropy = growth_index(t, m.operator_entropy; normalize = normalize),
        Γ_operator_log_neff = growth_index(t, log.(m.operator_neff_shannon .+ eps()); normalize = normalize),
        Γ_operator_neff_ipr = growth_index(t, m.operator_neff_ipr; normalize = normalize),

        Γ_majorana_weight_mean = growth_index(t, m.majorana_weight_mean; normalize = normalize),
        Γ_majorana_weight_entropy = growth_index(t, m.majorana_weight_entropy; normalize = normalize),
        Γ_majorana_weight_neff = growth_index(t, m.majorana_weight_neff; normalize = normalize),

        Γ_leakage_ge_4 = growth_index(t, m.leakage_ge_4; subtract_initial = false, normalize = normalize),
        Γ_leakage_ge_6 = growth_index(t, m.leakage_ge_6; subtract_initial = false, normalize = normalize),

        max_operator_entropy = max_growth(m.operator_entropy),
        max_majorana_weight_entropy = max_growth(m.majorana_weight_entropy),
        max_majorana_weight_mean = max_growth(m.majorana_weight_mean),
        max_leakage_ge_4 = maximum(m.leakage_ge_4),
        max_leakage_ge_6 = maximum(m.leakage_ge_6),
    )
end

function collect_gamma_scan(scan; normalize_gamma::Bool = false)
    Rs = sorted_Rs(scan)

    Γ_Sop = Float64[]
    Γ_Sw = Float64[]
    Γ_wmean = Float64[]
    Γ_P4 = Float64[]
    Γ_P6 = Float64[]

    max_P4 = Float64[]
    max_P6 = Float64[]

    for R in Rs
        out = scan[R].out
        γ = compute_gamma_indices(out; normalize = normalize_gamma)

        push!(Γ_Sop, γ.Γ_operator_entropy)
        push!(Γ_Sw, γ.Γ_majorana_weight_entropy)
        push!(Γ_wmean, γ.Γ_majorana_weight_mean)
        push!(Γ_P4, γ.Γ_leakage_ge_4)
        push!(Γ_P6, γ.Γ_leakage_ge_6)

        push!(max_P4, γ.max_leakage_ge_4)
        push!(max_P6, γ.max_leakage_ge_6)
    end

    return (
        R = Rs,
        Γ_operator_entropy = Γ_Sop,
        Γ_majorana_weight_entropy = Γ_Sw,
        Γ_majorana_weight_mean = Γ_wmean,
        Γ_leakage_ge_4 = Γ_P4,
        Γ_leakage_ge_6 = Γ_P6,
        max_leakage_ge_4 = max_P4,
        max_leakage_ge_6 = max_P6,
    )
end


function plot_gamma_scan(scan; normalize_gamma::Bool = false)
    g = collect_gamma_scan(scan; normalize_gamma = normalize_gamma)

    p = plot(
        xlabel = "R",
        ylabel = "growth index",
        title = "Operator-growth indices vs bond distance",
        lw = 2,
        marker = :circle,
    )

    plot!(p, g.R, g.Γ_operator_entropy; label = "Γ S_op", lw = 2, marker = :circle)
    plot!(p, g.R, g.Γ_majorana_weight_entropy; label = "Γ S_w", lw = 2, marker = :circle)
    plot!(p, g.R, g.Γ_majorana_weight_mean; label = "Γ <w_M>", lw = 2, marker = :circle)
    plot!(p, g.R, g.Γ_leakage_ge_4; label = "Γ P>=4", lw = 2, marker = :circle)
    plot!(p, g.R, g.Γ_leakage_ge_6; label = "Γ P>=6", lw = 2, marker = :circle)

    return p
end

# ============================================================
# Pauli Propagation and Metrics
# ============================================================
"""
Compute all operator-growth metrics from one propagated operator snapshot.
"""
function snapshot_metrics(W::PauliSum{N,T}; ϵ::Float64 = 1e-15) where {N,T}
    hist, norm2 = weight_profile(W; normalize = false)

    wm = weight_moments(hist; ϵ = ϵ)
    op = operator_entropy(W; ϵ = ϵ)

    return (
        norm2 = norm2,

        operator_entropy = op.entropy,
        operator_neff_shannon = op.neff_shannon,
        operator_neff_ipr = op.neff_ipr,

        majorana_weight_mean = wm.mean,
        majorana_weight_variance = wm.variance,
        majorana_weight_std = wm.std,
        majorana_weight_entropy = wm.entropy,
        majorana_weight_neff = wm.neff,

        leakage_ge_4 = high_weight_leakage(hist, 4),
        leakage_ge_6 = high_weight_leakage(hist, 6),

        majorana_weight_profile = normalized(hist),
        majorana_weight_profile_raw = hist,
    )
end

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

"""
Pauli propagation while tracking full operator entropy,
Majorana-weight distribution, Majorana-weight moments,
and high-weight leakage.
"""
function weight_dynamics_v2(
    O0::PauliSum{N,T},
    H::PauliSum{N,T},
    cfg::PPConfig;
    ϵ::Float64 = 1e-15,
) where {N,T}

    Ot = deepcopy(O0)

    weight_profiles = Vector{Vector{Float64}}()
    raw_weight_profiles = Vector{Vector{Float64}}()

    norm2 = Float64[]

    operator_entropy_series = Float64[]
    operator_neff_shannon_series = Float64[]
    operator_neff_ipr_series = Float64[]

    majorana_weight_mean_series = Float64[]
    majorana_weight_variance_series = Float64[]
    majorana_weight_entropy_series = Float64[]
    majorana_weight_neff_series = Float64[]

    leakage_ge_4_series = Float64[]
    leakage_ge_6_series = Float64[]

    # Initial snapshot
    m0 = snapshot_metrics(Ot; ϵ = ϵ)

    push!(weight_profiles, m0.majorana_weight_profile)
    push!(raw_weight_profiles, m0.majorana_weight_profile_raw)

    push!(norm2, m0.norm2)

    push!(operator_entropy_series, m0.operator_entropy)
    push!(operator_neff_shannon_series, m0.operator_neff_shannon)
    push!(operator_neff_ipr_series, m0.operator_neff_ipr)

    push!(majorana_weight_mean_series, m0.majorana_weight_mean)
    push!(majorana_weight_variance_series, m0.majorana_weight_variance)
    push!(majorana_weight_entropy_series, m0.majorana_weight_entropy)
    push!(majorana_weight_neff_series, m0.majorana_weight_neff)

    push!(leakage_ge_4_series, m0.leakage_ge_4)
    push!(leakage_ge_6_series, m0.leakage_ge_6)

    generators, coeffs = QuantumChemQC.gens_from_H(H)

    @printf("Total Pauli rotations: %d\n", length(coeffs) * cfg.n_steps)

    for step in 1:cfg.n_steps
        for j in eachindex(coeffs)
            θ = 2 * cfg.dt * coeffs[j]

            evolve!(Ot, generators[j], θ)
            QuantumChemQC.coeff_clip!(Ot; thresh = cfg.coeff_threshold)
        end

        m = snapshot_metrics(Ot; ϵ = ϵ)

        push!(weight_profiles, m.majorana_weight_profile)
        push!(raw_weight_profiles, m.majorana_weight_profile_raw)

        push!(norm2, m.norm2)

        push!(operator_entropy_series, m.operator_entropy)
        push!(operator_neff_shannon_series, m.operator_neff_shannon)
        push!(operator_neff_ipr_series, m.operator_neff_ipr)

        push!(majorana_weight_mean_series, m.majorana_weight_mean)
        push!(majorana_weight_variance_series, m.majorana_weight_variance)
        push!(majorana_weight_entropy_series, m.majorana_weight_entropy)
        push!(majorana_weight_neff_series, m.majorana_weight_neff)

        push!(leakage_ge_4_series, m.leakage_ge_4)
        push!(leakage_ge_6_series, m.leakage_ge_6)
    end

    tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

    return (
        tgrid = tgrid,

        snapshots = weight_profiles,
        raw_snapshots = raw_weight_profiles,
        raw_norms = norm2,

        metrics = (
            operator_entropy = operator_entropy_series,
            operator_neff_shannon = operator_neff_shannon_series,
            operator_neff_ipr = operator_neff_ipr_series,

            majorana_weight_mean = majorana_weight_mean_series,
            majorana_weight_variance = majorana_weight_variance_series,
            majorana_weight_entropy = majorana_weight_entropy_series,
            majorana_weight_neff = majorana_weight_neff_series,

            leakage_ge_4 = leakage_ge_4_series,
            leakage_ge_6 = leakage_ge_6_series,
        ),
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

    #out = weight_dynamics(O, H, cfg)
    out = weight_dynamics_v2(O, H, cfg)

    @printf("Initial norm² = %.12e\n", out.raw_norms[1])
    @printf("Final norm²   = %.12e\n", out.raw_norms[end])
    @printf("Relative norm² = %.12e\n", out.raw_norms[end] / out.raw_norms[1])

    @printf("Initial S_weight = %.8f\n", out.metrics.majorana_weight_entropy[1])
    @printf("Final S_weight   = %.8f\n", out.metrics.majorana_weight_entropy[end])

    @printf("Initial S_op = %.8f\n", out.metrics.operator_entropy[1])
    @printf("Final S_op   = %.8f\n", out.metrics.operator_entropy[end])

    @printf("Final <w_M> = %.8f\n", out.metrics.majorana_weight_mean[end])
    @printf("Final P_ge_4 = %.8f\n", out.metrics.leakage_ge_4[end])
    @printf("Final P_ge_6 = %.8f\n", out.metrics.leakage_ge_6[end])

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

function results_scan_table_v2(scan; normalize_gamma::Bool = false)
    Rs = sorted_Rs(scan)

    println()
    println("R        Γ_Sop       Γ_Sw        Γ_<wM>      Γ_P>=4      Γ_P>=6      max_P>=4    norm_ret")
    println("-"^115)

    for R in Rs
        out = scan[R].out
        γ = compute_gamma_indices(out; normalize = normalize_gamma)

        norm_ret = out.raw_norms[end] / out.raw_norms[1]

        @printf(
            "%-8.4f %-11.6f %-11.6f %-11.6f %-11.6f %-11.6f %-11.6f %-13.8f\n",
            R,
            γ.Γ_operator_entropy,
            γ.Γ_majorana_weight_entropy,
            γ.Γ_majorana_weight_mean,
            γ.Γ_leakage_ge_4,
            γ.Γ_leakage_ge_6,
            γ.max_leakage_ge_4,
            norm_ret,
        )
    end
end

function plot_metric_across_R(scan, metric_name::Symbol; ylabel::String = string(metric_name))
    Rs = sorted_Rs(scan)

    p = plot(
        xlabel = "time",
        ylabel = ylabel,
        title = string(metric_name),
        lw = 2,
    )

    for R in Rs
        out = scan[R].out
        y = getproperty(out.metrics, metric_name)
        label = @sprintf("R=%.4g", R)

        plot!(
            p,
            out.tgrid,
            y;
            label = label,
            lw = 2,
        )
    end

    return p
end

# ============================================================
# DATA INPUT AND PARAMETERS
# ============================================================

const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/tensors" # Canonical RHF
const CASES = [
    (0.50, joinpath(BASE, "RHF_H2_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "RHF_H2_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "RHF_H2_R_1p0_tensors.npz"))
    (1.50, joinpath(BASE, "RHF_H2_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "RHF_H2_R_2p0_tensors.npz"))
    (2.50, joinpath(BASE, "RHF_H2_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "RHF_H2_R_3p0_tensors.npz"))
]

#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/uhf_tensors" # Canonical UHF
const CASES = [
    (0.50, joinpath(BASE, "UHF_H2_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "UHF_H2_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "UHF_H2_R_1p0_tensors.npz"))
    (1.50, joinpath(BASE, "UHF_H2_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "UHF_H2_R_2p0_tensors.npz"))
    (2.50, joinpath(BASE, "UHF_H2_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "UHF_H2_R_3p0_tensors.npz"))
]=#

#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/nuhf_tensors" # NATURAL UHF
const CASES = [
    (0.50, joinpath(BASE, "Natural_UHF_H2_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "Natural_UHF_H2_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "Natural_UHF_H2_R_1p0_tensors.npz"))
    (1.50, joinpath(BASE, "Natural_UHF_H2_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "Natural_UHF_H2_R_2p0_tensors.npz"))
    (2.50, joinpath(BASE, "Natural_UHF_H2_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "Natural_UHF_H2_R_3p0_tensors.npz"))

]=#

#=const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/nuccsd_tensors" # Natural UCCSD
const CASES = [
    (0.50, joinpath(BASE, "Natural_UCCSD_H2_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "Natural_UCCSD_H2_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "Natural_UCCSD_H2_R_1p0_tensors.npz"))
    (1.50, joinpath(BASE, "Natural_UCCSD_H2_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "Natural_UCCSD_H2_R_2p0_tensors.npz"))
    (2.50, joinpath(BASE, "Natural_UCCSD_H2_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "Natural_UCCSD_H2_R_3p0_tensors.npz"))

]=#

cfg = PPConfig(
    500,       # n_steps
    0.1,     # dt
    1e-12,   # coeff_threshold
    false,   # normalize_weights
)

scan = run_distance_scan(
    CASES,
    cfg;
    z_qubit = 3,
    block = false,
    NOI = false,
    UNRESTRICTED = false,
)

#p_summary = plot_distance_scan_summary(scan)
#display(p_summary)
#savefig(p_summary, "majorana_distance_scan_summary.png")
#results_scan_table(scan)

p_heat = plot_selected_heatmaps(scan; selected_Rs = [0.50, 0.7414, 1.50, 3.00])
display(p_heat)
savefig(p_heat, "majorana_weight_heatmaps_selected_R.png")

p_Sop = plot_metric_across_R(scan, :operator_entropy; ylabel = "operator entropy")
display(p_Sop)
savefig(p_Sop, "operator_entropy_scan.png")

p_Sw = plot_metric_across_R(scan, :majorana_weight_entropy; ylabel = "Majorana-weight entropy")
display(p_Sw)
savefig(p_Sw, "majorana_weight_entropy_scan.png")

p_w = plot_metric_across_R(scan, :majorana_weight_mean; ylabel = "<w_M>")
display(p_w)
savefig(p_w, "majorana_weight_mean_scan.png")

p_leak4 = plot_metric_across_R(scan, :leakage_ge_4; ylabel = "P(w_M >= 4)")
display(p_leak4)
savefig(p_leak4, "majorana_leakage_ge4_scan.png")

p_gamma = plot_gamma_scan(scan)
display(p_gamma)
savefig(p_gamma, "gamma_indices_vs_R.png")

results_scan_table_v2(scan)