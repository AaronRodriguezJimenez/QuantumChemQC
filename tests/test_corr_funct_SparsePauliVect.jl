using PauliOperators
using QuantumChemQC
using LinearAlgebra
using Printf
using NPZ
using Plots
using Test


# ============================================================
# 1. SETTINGS
# ============================================================

function run_H()
    data_path = joinpath(@__DIR__, "h2-RHF_test_ham.npz")
    data = npzread(data_path)

    H0 = data["hc"][1]
    H1 = data["h1e"]
    H2 = data["h2e"]

    println("H0 shape: ", size(H0))
    println("H0 type:  ", typeof(H0))

    println("H1 shape: ", size(H1))
    println("H1 type:  ", typeof(H1))

    println("H2 shape: ", size(H2))
    println("H2 type:  ", typeof(H2))

    Norbs = size(H1, 1)

    println("\nConstructing molecular Hamiltonian...")

    @time H = QuantumChemQC.molecular_hamiltonian(
        Norbs,
        data_path;
        NOI = false,
        block = true,
    )

    return H
end


# ------------------------------------------------------------
# Initial observable
# ------------------------------------------------------------

O_dict = PauliSum(Pauli(4, X=[1, 3]))
O_dict += Pauli(4, X=[2, 3])


# ------------------------------------------------------------
# Hamiltonian and reference state
# ------------------------------------------------------------

H = run_H()

ket, _ = QuantumChemQC.string_to_ket("0000")


# ------------------------------------------------------------
# Evolution settings
# ------------------------------------------------------------

n_intervals = 200
total_time = 50.0
dt = total_time / n_intervals

evol_thresh = 0.01
max_weight = 4

# First-order Trotter decomposition
trotter_order = 1
n_trotter = 1

# Keep window=1 for direct PauliSum/SparsePauliVector comparison.
#
# After validating the results, test values such as 4, 8, or 16.
sparse_window = 1


# ------------------------------------------------------------
# Sparse representation of the initial observable
#
# Float64 is appropriate here because O is Hermitian and its
# Pauli-basis coefficients are real.
# ------------------------------------------------------------

O_sparse = SparsePauliVector(
    O_dict;
    T = Float64,
    capacity_factor = 4.0,
    append_factor = 2.0,
    min_capacity = 256,
)


# ------------------------------------------------------------
# Initial conversion tests
# ------------------------------------------------------------

@test PauliSum(O_sparse) == O_dict
@test PauliOperators.check_spv(O_sparse)

println()
println("Hamiltonian terms:         ", length(H))
println("Initial observable terms:  ", length(O_dict))
println("Number of intervals:       ", n_intervals)
println("Total time:                ", total_time)
println("Time step:                 ", dt)
println("Coefficient threshold:     ", evol_thresh)
println("Maximum Pauli weight:      ", max_weight)
println("Sparse merge window:       ", sparse_window)


# ============================================================
# 2. CALCULATION FUNCTIONS AND HELPERS
# ============================================================

"""
Compute

    C = <ket| left * right |ket>

without explicitly constructing the product `left * right`.

This works with both PauliSum and SparsePauliVector.
"""
function correlation_value(left, right, ket)
    value = zero(ComplexF64)

    for (p, cp) in left
        for (q, cq) in right
            value += cp * cq * expectation_value(p * q, ket)
        end
    end

    return value
end


# ------------------------------------------------------------
# Complex correlation correction
#
# PauliOperators accumulates:
#
#     after_truncation - before_truncation
#
# Therefore, the corrected signal is:
#
#     current_signal - accumulated_change
# ------------------------------------------------------------

mutable struct CorrelationCorrection{N,L,K} <: CorrectionAccumulator
    left_operator::L
    ket::K
    accumulated_change::ComplexF64
end


function CorrelationCorrection(
    left_operator::Union{PauliSum{N},SparsePauliVector{N}},
    ket::Ket{N},
) where {N}

    return CorrelationCorrection{
        N,
        typeof(left_operator),
        typeof(ket),
    }(
        left_operator,
        ket,
        zero(ComplexF64),
    )
end


"""
Measure the complex correlation before and after truncation.
"""
function PauliOperators._measure(
    evolving_operator::Union{PauliSum{N},SparsePauliVector{N}},
    correction::CorrelationCorrection{N},
) where {N}

    return correlation_value(
        correction.left_operator,
        evolving_operator,
        correction.ket,
    )
end


"""
Accumulate the truncation-induced change:

    after - before
"""
function PauliOperators._accumulate!(
    correction::CorrelationCorrection,
    before,
    after,
)
    correction.accumulated_change += after - before
    return nothing
end


# ------------------------------------------------------------
# Dictionary evolution for one complete Trotter interval
# ------------------------------------------------------------

function evolve_trotter_interval!(
    Wt::PauliSum,
    generators,
    angles,
    truncation,
    correction;
    window::Int = 1,
)
    window == 1 || throw(
        ArgumentError(
            "PauliSum uses truncation after every rotation; " *
            "set window=1.",
        ),
    )

    @inbounds for i in eachindex(angles)
        evolve!(Wt, generators[i], angles[i])
        truncate!(Wt, truncation, correction)
    end

    return Wt
end


# ------------------------------------------------------------
# SparsePauliVector evolution for one complete Trotter interval
# ------------------------------------------------------------

function evolve_trotter_interval!(
    Wt::SparsePauliVector,
    generators,
    angles,
    truncation,
    correction;
    window::Int = 1,
)
    evolve!(
        Wt,
        generators,
        angles;
        window = window,
        truncation = truncation,
        correction = correction,
    )

    return Wt
end


"""
Compute

    C(t) = <ket| O(0) O(t) |ket>

using either PauliSum or SparsePauliVector.

The Trotter angles returned by `trotterize` are passed directly to
`evolve!`. They must not be multiplied by `2dt` again.
"""
function QSP_evolution_op_v3(
    ket::Ket{N},
    initial_operator::Union{PauliSum{N},SparsePauliVector{N}},
    H,
    n_intervals::Int,
    dt::Real;
    thresh::Real = 1e-3,
    max_weight::Int = N,
    window::Int = 1,
    trotter_order::Int = 1,
    n_trotter::Int = 1,
) where {N}

    n_intervals > 0 ||
        throw(ArgumentError("n_intervals must be positive"))

    dt > 0 ||
        throw(ArgumentError("dt must be positive"))

    window > 0 ||
        throw(ArgumentError("window must be positive"))

    # Fixed left operator and evolving operator
    W = copy(initial_operator)
    Wt = copy(initial_operator)

    # --------------------------------------------------------
    # Construct one Trotter interval
    #
    # Returned angles are already the rotation angles expected
    # by evolve!.
    # --------------------------------------------------------

    generators, angles = trotterize(
        H,
        dt;
        n_trotter = n_trotter,
        order = trotter_order,
    )

    nrotations = length(angles)

    println()
    println("Representation:          ", typeof(initial_operator))
    println("Rotations per interval:  ", nrotations)
    println("Total rotations:         ", nrotations * n_intervals)

    # --------------------------------------------------------
    # Combined pruning rule
    # --------------------------------------------------------

    truncation = CompositeTruncation(
        CoeffTruncation(Float64(thresh)),
        WeightTruncation(max_weight),
    )

    correction = CorrelationCorrection(W, ket)

    # --------------------------------------------------------
    # Preallocated output
    # --------------------------------------------------------

    rCtvals = Vector{Float64}(undef, n_intervals + 1)
    iCtvals = Vector{Float64}(undef, n_intervals + 1)

    raw_rCtvals = Vector{Float64}(undef, n_intervals + 1)
    raw_iCtvals = Vector{Float64}(undef, n_intervals + 1)

    norms = Vector{Float64}(undef, n_intervals + 1)
    term_counts = Vector{Int}(undef, n_intervals + 1)
    correction_values =
        Vector{ComplexF64}(undef, n_intervals + 1)

    tgrid = collect(
        range(
            0.0;
            stop = n_intervals * dt,
            length = n_intervals + 1,
        ),
    )

    # --------------------------------------------------------
    # t = 0
    # --------------------------------------------------------

    C0 = correlation_value(W, Wt, ket)

    rCtvals[1] = real(C0)
    iCtvals[1] = imag(C0)

    raw_rCtvals[1] = real(C0)
    raw_iCtvals[1] = imag(C0)

    norms[1] = norm(Wt)
    term_counts[1] = length(Wt)
    correction_values[1] = zero(ComplexF64)

    # --------------------------------------------------------
    # Time propagation
    # --------------------------------------------------------

    for interval in 1:n_intervals
        evolve_trotter_interval!(
            Wt,
            generators,
            angles,
            truncation,
            correction;
            window = window,
        )

        raw_correlation = correlation_value(W, Wt, ket)

        # PauliOperators stores after - before, so subtract it.
        corrected_correlation =
            raw_correlation - correction.accumulated_change

        idx = interval + 1

        raw_rCtvals[idx] = real(raw_correlation)
        raw_iCtvals[idx] = imag(raw_correlation)

        rCtvals[idx] = real(corrected_correlation)
        iCtvals[idx] = imag(corrected_correlation)

        norms[idx] = norm(Wt)
        term_counts[idx] = length(Wt)

        correction_values[idx] =
            correction.accumulated_change
    end

    return (
        ReCt = rCtvals,
        ImCt = iCtvals,
        RawReCt = raw_rCtvals,
        RawImCt = raw_iCtvals,
        intervals = tgrid,
        norms = norms,
        term_counts = term_counts,
        corrections = correction_values,
        final_operator = Wt,
        rotations_per_interval = nrotations,
    )
end


# ------------------------------------------------------------
# Reference energy
# ------------------------------------------------------------

function compute_ref_expval(H, ket)
    return expectation_value(H, ket)
end


# ============================================================
# 3. TEST, COMPARISON, AND COMPUTATIONAL COST
# ============================================================

println()
println("============================================================")
println("COMPILATION WARM-UP")
println("============================================================")

# Do not include Julia compilation in the benchmark.
QSP_evolution_op_v3(
    ket,
    O_dict,
    H,
    2,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
    window = 1,
    trotter_order = trotter_order,
    n_trotter = n_trotter,
)

QSP_evolution_op_v3(
    ket,
    O_sparse,
    H,
    2,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
    window = sparse_window,
    trotter_order = trotter_order,
    n_trotter = n_trotter,
)


# ------------------------------------------------------------
# Dictionary calculation
# ------------------------------------------------------------

println()
println("============================================================")
println("PAULISUM CALCULATION")
println("============================================================")

GC.gc()

dict_stats = @timed QSP_evolution_op_v3(
    ket,
    O_dict,
    H,
    n_intervals,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
    window = 1,
    trotter_order = trotter_order,
    n_trotter = n_trotter,
)

dict_result = dict_stats.value


# ------------------------------------------------------------
# Sparse calculation
# ------------------------------------------------------------

println()
println("============================================================")
println("SPARSEPAULIVECTOR CALCULATION")
println("============================================================")

GC.gc()

sparse_stats = @timed QSP_evolution_op_v3(
    ket,
    O_sparse,
    H,
    n_intervals,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
    window = sparse_window,
    trotter_order = trotter_order,
    n_trotter = n_trotter,
)

sparse_result = sparse_stats.value


# ------------------------------------------------------------
# Numerical comparison
# ------------------------------------------------------------

delta_re = sparse_result.ReCt .- dict_result.ReCt
delta_im = sparse_result.ImCt .- dict_result.ImCt
delta_norm = sparse_result.norms .- dict_result.norms

max_re_error = maximum(abs.(delta_re))
max_im_error = maximum(abs.(delta_im))
max_norm_error = maximum(abs.(delta_norm))

println()
println("============================================================")
println("NUMERICAL COMPARISON")
println("============================================================")

@printf("Maximum Re(C(t)) error:   %.8e\n", max_re_error)
@printf("Maximum Im(C(t)) error:   %.8e\n", max_im_error)
@printf("Maximum norm error:       %.8e\n", max_norm_error)

println(
    "Final PauliSum terms:     ",
    length(dict_result.final_operator),
)

println(
    "Final sparse terms:       ",
    length(sparse_result.final_operator),
)

println(
    "Maximum PauliSum terms:   ",
    maximum(dict_result.term_counts),
)

println(
    "Maximum sparse terms:     ",
    maximum(sparse_result.term_counts),
)


# Exact comparison is expected only for window=1.
if sparse_window == 1
    @test dict_result.intervals == sparse_result.intervals

    @test isapprox(
        dict_result.ReCt,
        sparse_result.ReCt;
        atol = 1e-8,
        rtol = 1e-8,
    )

    @test isapprox(
        dict_result.ImCt,
        sparse_result.ImCt;
        atol = 1e-8,
        rtol = 1e-8,
    )

    @test isapprox(
        dict_result.norms,
        sparse_result.norms;
        atol = 1e-8,
        rtol = 1e-8,
    )

    @test dict_result.term_counts ==
          sparse_result.term_counts

    @test PauliOperators.check_spv(
        sparse_result.final_operator,
    )

    println("Numerical parity tests passed.")
else
    println(
        "Sparse window is greater than one; " *
        "the truncation cadence differs from PauliSum.",
    )
end


# ------------------------------------------------------------
# Performance comparison
# ------------------------------------------------------------

dict_memory_mib = dict_stats.bytes / 2.0^20
sparse_memory_mib = sparse_stats.bytes / 2.0^20

speedup = dict_stats.time / sparse_stats.time

allocation_reduction =
    dict_stats.bytes / max(sparse_stats.bytes, 1)

println()
println("============================================================")
println("COMPUTATIONAL COST")
println("============================================================")

@printf(
    "%-24s %12s %18s %14s\n",
    "Representation",
    "Time (s)",
    "Allocated (MiB)",
    "GC time (s)",
)

@printf(
    "%-24s %12.6f %18.3f %14.6f\n",
    "PauliSum",
    dict_stats.time,
    dict_memory_mib,
    dict_stats.gctime,
)

@printf(
    "%-24s %12.6f %18.3f %14.6f\n",
    "SparsePauliVector",
    sparse_stats.time,
    sparse_memory_mib,
    sparse_stats.gctime,
)

println()

@printf(
    "Sparse speedup:             %.3fx\n",
    speedup,
)

@printf(
    "Allocation reduction:       %.3fx\n",
    allocation_reduction,
)


# ============================================================
# RESULTS AND PLOTS
# ============================================================

rRES = sparse_result.ReCt
iRES = sparse_result.ImCt
tgrid = sparse_result.intervals

println()
println("* * * * Number of snapshots collected: ", length(rRES))


# ------------------------------------------------------------
# Print correlation function
# ------------------------------------------------------------

println()
println("- - - Correlation function - - -")

@printf(
    "%10s %18s %18s %18s\n",
    "time",
    "Re(C(t))",
    "Im(C(t))",
    "|C(t)|",
)

for i in eachindex(tgrid)
    Cabs = hypot(rRES[i], iRES[i])

    @printf(
        "%10.4f %18.10f %18.10f %18.10f\n",
        tgrid[i],
        rRES[i],
        iRES[i],
        Cabs,
    )
end


# ------------------------------------------------------------
# PauliSum versus SparsePauliVector plot
# ------------------------------------------------------------

comparison_plot = plot(
    dict_result.intervals,
    dict_result.ReCt;
    lw = 2,
    label = "Re(C(t)) PauliSum",
    xlabel = "Time",
    ylabel = "< O(0) O(t) >",
)

plot!(
    comparison_plot,
    sparse_result.intervals,
    sparse_result.ReCt;
    lw = 2,
    linestyle = :dash,
    label = "Re(C(t)) SparsePauliVector",
)

plot!(
    comparison_plot,
    dict_result.intervals,
    dict_result.ImCt;
    lw = 2,
    label = "Im(C(t)) PauliSum",
)

plot!(
    comparison_plot,
    sparse_result.intervals,
    sparse_result.ImCt;
    lw = 2,
    linestyle = :dash,
    label = "Im(C(t)) SparsePauliVector",
)

savefig(
    comparison_plot,
    joinpath(
        @__DIR__,
        "QSP_H2_PauliSum_vs_Sparse_$evol_thresh.pdf",
    ),
)


# ------------------------------------------------------------
# Numerical-error plot
# ------------------------------------------------------------

error_plot = plot(
    dict_result.intervals,
    abs.(delta_re);
    lw = 2,
    label = "|Δ Re(C(t))|",
    xlabel = "Time",
    ylabel = "Absolute difference",
)

plot!(
    error_plot,
    dict_result.intervals,
    abs.(delta_im);
    lw = 2,
    label = "|Δ Im(C(t))|",
)

savefig(
    error_plot,
    joinpath(
        @__DIR__,
        "QSP_H2_PauliSum_vs_Sparse_error_$evol_thresh.pdf",
    ),
)


# ------------------------------------------------------------
# Sparse correlation plot
# ------------------------------------------------------------

correlation_plot = plot(
    tgrid,
    rRES;
    lw = 2,
    label = "Re(C(t)), threshold=$evol_thresh",
    xlabel = "Time",
    ylabel = "< O(0) O(t) >",
)

plot!(
    correlation_plot,
    tgrid,
    iRES;
    lw = 2,
    label = "Im(C(t)), threshold=$evol_thresh",
)

savefig(
    correlation_plot,
    joinpath(
        @__DIR__,
        "QSP_H2_deterministic_sparse_$evol_thresh.pdf",
    ),
)


# ------------------------------------------------------------
# Operator-population plot
# ------------------------------------------------------------

population_plot = plot(
    dict_result.intervals,
    dict_result.term_counts;
    lw = 2,
    label = "PauliSum",
    xlabel = "Time",
    ylabel = "Number of Pauli terms",
)

plot!(
    population_plot,
    sparse_result.intervals,
    sparse_result.term_counts;
    lw = 2,
    linestyle = :dash,
    label = "SparsePauliVector",
)

savefig(
    population_plot,
    joinpath(
        @__DIR__,
        "QSP_H2_operator_population_$evol_thresh.pdf",
    ),
)


# ============================================================
# PHASE CORRECTION
# ============================================================

println()
println("Initial state:")
display(ket)

ref_expval = compute_ref_expval(H, ket)

println()
println("Reference energy: ", ref_expval)

# This phase correction is appropriate only when the reference
# state is an eigenstate of H.
Ek = real(ref_expval)

signal = rRES .+ 1im .* iRES
phase = exp.(1im .* Ek .* tgrid)

F = phase .* signal


# ------------------------------------------------------------
# Corrected signal plot
# ------------------------------------------------------------

phase_plot = plot(
    tgrid,
    real.(F);
    lw = 2,
    label = "Re(F(t)), threshold=$evol_thresh",
    xlabel = "Time",
    ylabel = "exp(i E₀ t) < O(0) O(t) >",
)

plot!(
    phase_plot,
    tgrid,
    imag.(F);
    lw = 2,
    label = "Im(F(t)), threshold=$evol_thresh",
)

savefig(
    phase_plot,
    joinpath(
        @__DIR__,
        "QSP_H2_deterministic_F_sparse_$evol_thresh.pdf",
    ),
)


# ------------------------------------------------------------
# Print phase-corrected signal
# ------------------------------------------------------------

println()
println("- - - Phase-corrected signal - - -")

@printf(
    "%10s %18s %18s\n",
    "time",
    "Re(F(t))",
    "Im(F(t))",
)

for i in eachindex(tgrid)
    @printf(
        "%10.4f %18.10f %18.10f\n",
        tgrid[i],
        real(F[i]),
        imag(F[i]),
    )
end