# SparsePauliVector is a sparse representation of an N-qubit operator
# written as a sum of Pauli strings. 
# It is intended as a high-performance, allocation-conscious alternative
# to the Dict-backed PauliSum. Instead of stoing each pauli as a dictionary key
# it stores packed Pauli keys and coefficients in stored parallel arrays.

# Here we test the SparsePauliVector type and its methods, and compare with PauliSum
using QuantumChemQC
using PauliOperators
using Printf
using NPZ
using Plots
using Test

# ============================================================
# 1. SETTINGS
# ============================================================
function run_H()
    # Get Molecular Hamiltonian
    data_path =  joinpath(@__DIR__, "h2-RHF_test_ham.npz")
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
    
    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=true)

    return H
end

# Initial operator, ket and Hamiltonian
O = Pauli(4, X=[1,3])
O = PauliSum(O)
O += Pauli(4, X=[2,3])
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("0000")
#display(H)

# Time evolution parameters
n_intervals = 200
t = 50.0
dt = t/n_intervals

#Pruning parameters
evol_thresh = 0.01
max_weight = 4 

println()
println("Hamiltonian terms: ", length(H))
println("Initial operator terms: ", length(O))
println("Time step dt: ", dt)
println("Number of intervals: ", n_intervals)
println("Pruning threshold: ", evol_thresh)

# ------------------------------------------------------------
# Create the SparsePauliVector versions
# ------------------------------------------------------------

O_sparse = SparsePauliVector(O)

H_sparse = SparsePauliVector(H)

# ------------------------------------------------------------
# Generate the Hamiltonian rotations once
#
# Both calculations must use the same generator order because
# first-order Trotterization depends on rotation ordering.
# ------------------------------------------------------------>

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


generators, angles = gens_from_H(H)
pbs = PauliBasis.(generators)

println("Rotations per interval: ", length(angles))
println("Total rotations: ", length(angles) * n_intervals)

# ============================================================
# 2. CALCULATION FUNCTIONS AND HELPERS
# ============================================================
"""
Perform Pauli propagation using either:

    PauliSum
    SparsePauliVector

The Hamiltonian generators and angles are constructed before calling
this function, ensuring that both representations use exactly the same
Trotter sequence.
"""
function QSP_evolution_op(
    ket,
    o,
    pbs,
    angles,
    n_intervals,
    dt;
    thresh::Float64 = 1e-3,
    max_weight::Int = 4,
)
    # copy works for both PauliSum and SparsePauliVector
    Wt = copy(o)
    W = copy(o)

    rCtvals = Vector{Float64}(undef, n_intervals + 1)
    iCtvals = Vector{Float64}(undef, n_intervals + 1)
    term_counts = Vector{Int}(undef, n_intervals + 1)

    # --------------------------------------------------------
    # t = 0
    # --------------------------------------------------------

    WW = W * W
    expval0 = expectation_value(WW, ket)

    rCtvals[1] = real(expval0)
    iCtvals[1] = imag(expval0)
    term_counts[1] = length(Wt)

    nt = length(angles)

    @assert length(pbs) == nt

    # --------------------------------------------------------
    # Time propagation
    # --------------------------------------------------------

    for ki in 1:n_intervals
        accumulated_error = zero(expval0)

        for j in 1:nt
            theta = 2 * dt * angles[j]
            pb = pbs[j]

            # Heisenberg evolution
            evolve!(Wt, pb, theta)

            # Remove numerical noise
            PauliOperators.coeff_clip!(Wt, 1e-12)

            # Observable before pruning
            WWt = W * Wt
            e1 = expectation_value(WWt, ket)

            # Apply pruning
            PauliOperators.coeff_clip!(Wt, thresh)
            PauliOperators.weight_clip!(Wt, max_weight)

            # Observable after pruning
            WWt = W * Wt
            e2 = expectation_value(WWt, ket)

            # Same pruning-error correction used in the original code
            accumulated_error += e2 - e1
        end

        # Final value for the current interval
        WWt = W * Wt
        expval = expectation_value(WWt, ket) + accumulated_error

        idx = ki + 1

        rCtvals[idx] = real(expval)
        iCtvals[idx] = imag(expval)
        term_counts[idx] = length(Wt)
    end

    tgrid = collect(
        range(
            0.0;
            stop = n_intervals * dt,
            length = n_intervals + 1,
        ),
    )

    return (
        ReCt = rCtvals,
        ImCt = iCtvals,
        intervals = tgrid,
        term_counts = term_counts,
        final_operator = Wt,
    )
end


"""
Run a calculation and measure elapsed time and allocated memory.

This should only be called after compilation warm-up.
"""
function measure_calculation(f)
    GC.gc()
    return @timed f()
end

# ============================================================
# 3. TEST AND PERFORMANCE COMPARISON
# ============================================================

println()
println("Compiling PauliSum calculation...")

# Short warm-up run so compilation is not included in the benchmark
QSP_evolution_op(
    ket,
    O,
    pbs,
    angles,
    2,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
)

println("Compiling SparsePauliVector calculation...")

QSP_evolution_op(
    ket,
    O_sparse,
    pbs,
    angles,
    2,
    dt;
    thresh = evol_thresh,
    max_weight = max_weight,
)


# ------------------------------------------------------------
# PauliSum calculation
# ------------------------------------------------------------

println()
println("Running PauliSum calculation...")

paulisum_stats = measure_calculation() do
    QSP_evolution_op(
        ket,
        O,
        pbs,
        angles,
        n_intervals,
        dt;
        thresh = evol_thresh,
        max_weight = max_weight,
    )
end

result_paulisum = paulisum_stats.value


# ------------------------------------------------------------
# SparsePauliVector calculation
# ------------------------------------------------------------

println("Running SparsePauliVector calculation...")

sparse_stats = measure_calculation() do
    QSP_evolution_op(
        ket,
        O_sparse,
        pbs,
        angles,
        n_intervals,
        dt;
        thresh = evol_thresh,
        max_weight = max_weight,
    )
end

result_sparse = sparse_stats.value


# ------------------------------------------------------------
# Numerical comparison
# ------------------------------------------------------------

delta_re = result_sparse.ReCt .- result_paulisum.ReCt
delta_im = result_sparse.ImCt .- result_paulisum.ImCt

max_re_error = maximum(abs.(delta_re))
max_im_error = maximum(abs.(delta_im))


@test result_sparse.intervals == result_paulisum.intervals

@test isapprox(
    result_sparse.ReCt,
    result_paulisum.ReCt;
    atol = 1e-8,
    rtol = 1e-8,
)

@test isapprox(
    result_sparse.ImCt,
    result_paulisum.ImCt;
    atol = 1e-8,
    rtol = 1e-8,
)

@test PauliOperators.check_spv(result_sparse.final_operator)


# ------------------------------------------------------------
# Performance comparison
# ------------------------------------------------------------

paulisum_memory_mb = paulisum_stats.bytes / 2.0^20
sparse_memory_mb = sparse_stats.bytes / 2.0^20

speedup = paulisum_stats.time / sparse_stats.time

allocation_ratio = if sparse_stats.bytes > 0
    paulisum_stats.bytes / sparse_stats.bytes
else
    Inf
end


println()
println("============================================================")
println("NUMERICAL COMPARISON")
println("============================================================")

@printf(
    "Maximum error in Re(C(t)):      %.6e\n",
    max_re_error,
)

@printf(
    "Maximum error in Im(C(t)):      %.6e\n",
    max_im_error,
)

println(
    "Final PauliSum term count:       ",
    length(result_paulisum.final_operator),
)

println(
    "Final SparsePauliVector terms:   ",
    length(result_sparse.final_operator),
)


println()
println("============================================================")
println("COMPUTATIONAL COST")
println("============================================================")

@printf(
    "%-24s %12s %18s %12s\n",
    "Representation",
    "Time (s)",
    "Allocated (MiB)",
    "GC time (s)",
)

@printf(
    "%-24s %12.6f %18.3f %12.6f\n",
    "PauliSum",
    paulisum_stats.time,
    paulisum_memory_mb,
    paulisum_stats.gctime,
)

@printf(
    "%-24s %12.6f %18.3f %12.6f\n",
    "SparsePauliVector",
    sparse_stats.time,
    sparse_memory_mb,
    sparse_stats.gctime,
)

println()

@printf(
    "SparsePauliVector speedup:       %.3fx\n",
    speedup,
)

@printf(
    "Allocation reduction factor:    %.3fx\n",
    allocation_ratio,
)

# ------------------------------------------------------------
# Plot the two calculations
# ------------------------------------------------------------

comparison_plot = plot(
    result_paulisum.intervals,
    result_paulisum.ReCt;
    lw = 2,
    label = "Re(C(t)) PauliSum",
    xlabel = "Time",
    ylabel = "< O(0) O(t) >",
)

plot!(
    comparison_plot,
    result_sparse.intervals,
    result_sparse.ReCt;
    lw = 2,
    linestyle = :dash,
    label = "Re(C(t)) SparsePauliVector",
)

plot!(
    comparison_plot,
    result_paulisum.intervals,
    result_paulisum.ImCt;
    lw = 2,
    label = "Im(C(t)) PauliSum",
)

plot!(
    comparison_plot,
    result_sparse.intervals,
    result_sparse.ImCt;
    lw = 2,
    linestyle = :dash,
    label = "Im(C(t)) SparsePauliVector",
)

savefig(
    comparison_plot,
    joinpath(
        @__DIR__,
        "PauliSum_vs_SparsePauliVector.pdf",
    ),
)

difference_plot = plot(
    result_paulisum.intervals,
    abs.(delta_re);
    lw = 2,
    label = "|Δ Re(C(t))|",
    xlabel = "Time",
    ylabel = "Absolute difference",
)

plot!(
    difference_plot,
    result_paulisum.intervals,
    abs.(delta_im);
    lw = 2,
    label = "|Δ Im(C(t))|",
)

savefig(
    difference_plot,
    joinpath(
        @__DIR__,
        "PauliSum_vs_SparsePauliVector_error.pdf",
    ),
)