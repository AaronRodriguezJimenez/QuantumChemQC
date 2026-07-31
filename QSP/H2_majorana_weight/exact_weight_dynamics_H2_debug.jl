#
# Majorana-weight dynamics for H2 using a matrix-free Arnoldi-Liouville formulation.
#
# This script intentionally defines NO new top-level structs and NO module.
# It is safe to re-include in a VS Code/REPL session because functions can be
# redefined, while Julia structs cannot.
#
# Main propagation path:
#     vt <- exp(dt * L) * vt
# where L(O) = i[H, O].
# The exponential-times-vector is approximated by a non-Hermitian Arnoldi Krylov
# method without forming exp(dt * L) and without using ExponentialUtilities.
#

using QuantumChemQC
using PauliOperators
using LinearAlgebra
using SparseArrays
using Statistics
using Printf
using Random
using Plots
using Plots.PlotMeasures
using NPZ

# ============================================================
# Configuration
# ============================================================

"""
    make_config(; kwargs...)

Configuration as a NamedTuple, not a struct. This avoids invalid redefinition
errors when re-including this script in the same Julia session.

Fields:
    n_steps::Int
    dt::Float64
    coeff_threshold::Float64
    normalize_weights::Bool
    propagation::Symbol
        :arnoldi_matrixfree  recommended larger-system path
        :arnoldi_sparse      builds sparse L, then Arnoldi expv
        :dense_exp           small-system reference only
    krylov_m::Int
        Arnoldi subspace dimension. Increase this if norm drift is too large.
    krylov_breakdown_tol::Float64
        Arnoldi breakdown tolerance.
    reorthogonalize::Bool
        Use a second Gram-Schmidt pass for stability.
"""
function make_config(;
    n_steps::Int,
    dt::Real,
    coeff_threshold::Real = 1e-12,
    normalize_weights::Bool = false,
    propagation::Symbol = :arnoldi_matrixfree,
    krylov_m::Int = 40,
    krylov_breakdown_tol::Real = 1e-13,
    reorthogonalize::Bool = true,
)
    return (
        n_steps = n_steps,
        dt = Float64(dt),
        coeff_threshold = Float64(coeff_threshold),
        normalize_weights = normalize_weights,
        propagation = propagation,
        krylov_m = krylov_m,
        krylov_breakdown_tol = Float64(krylov_breakdown_tol),
        reorthogonalize = reorthogonalize,
    )
end

# ============================================================
# Metrics
# ============================================================

function normalized(w::AbstractVector)
    s = sum(w)
    s <= 0 && return copy(w)
    return w ./ s
end

"""
Entropy of the Majorana-weight distribution.
This is not the entropy over individual Pauli strings.
"""
function weight_entropy(w::AbstractVector; eps::Float64 = 1e-15)
    p = normalized(w)
    p = p[p .> eps]
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
# Pauli coefficient-vector utilities
# ============================================================

@inline function pauli_index_from_bits(z::Int128, x::Int128, ::Val{N}) where {N}
    return Int(z + (x << N) + 1)
end

@inline function pauli_index(p::PauliBasis{N}) where {N}
    return pauli_index_from_bits(p.z, p.x, Val(N))
end

@inline function pauli_weight(p::PauliBasis)
    return count_ones(p.z | p.x)
end

function coeff_vector(ps::PauliSum{N}; T = ComplexF64) where {N}
    D = Int(4)^N
    v = zeros(T, D)

    for (p, c) in ps
        v[pauli_index(p)] += convert(T, c)
    end

    return v
end

function singleton_paulisum(p::PauliBasis{N}; T = ComplexF64) where {N}
    out = PauliSum(N, T)
    out[p] = one(T)
    return out
end

# ============================================================
# Precomputed basis tables as NamedTuples, not structs
# ============================================================

function pauli_tables(::Val{N}) where {N}
    D = Int(4)^N

    zs = Vector{Int128}(undef, D)
    xs = Vector{Int128}(undef, D)
    nys = Vector{Int}(undef, D)
    mws = Vector{Int}(undef, D)

    for p in PauliBasis{N}
        j = pauli_index(p)

        zs[j] = p.z
        xs[j] = p.x
        nys[j] = count_ones(p.z & p.x)
        mws[j] = QuantumChemQC.majorana_weight(p)
    end

    return (
        N = N,
        D = D,
        z = zs,
        x = xs,
        ny = nys,
        majorana_weight = mws,
    )
end

# ============================================================
# Majorana-weight distribution from coefficient vector
# ============================================================

"""
    weight_profile(v, tables; normalize=true, coeff_threshold=0.0)

Builds

    p_k(t) = sum_{P: majorana_weight(P) = k} |c_P(t)|^2

from the Pauli coefficient vector v.

Important:
    coeff_threshold is diagnostic only. It does not prune the propagated
    vector, because pruning would perturb the Krylov evolution.
"""
function weight_profile(
    v::AbstractVector,
    tables;
    normalize::Bool = true,
    coeff_threshold::Float64 = 0.0,
)
    @assert length(v) == tables.D

    hist = zeros(Float64, 2 * tables.N + 1)
    total = 0.0

    @inbounds for j in eachindex(v)
        c = v[j]

        abs(c) > coeff_threshold || continue

        a2 = abs2(c)
        k = tables.majorana_weight[j]

        hist[k + 1] += a2
        total += a2
    end

    if normalize && total > 0
        hist ./= total
    end

    return hist, total
end

# Compatibility method matching the earlier call style.
function weight_profile(
    v::AbstractVector,
    ::Val{N};
    normalize::Bool = true,
    coeff_threshold::Float64 = 0.0,
) where {N}
    return weight_profile(
        v,
        pauli_tables(Val(N));
        normalize = normalize,
        coeff_threshold = coeff_threshold,
    )
end

# ============================================================
# Hamiltonian term preprocessing
# ============================================================

"""
    hamiltonian_terms(H)

Preprocess H into z/x/coefficient arrays. Identity terms are skipped because
they commute with every operator and do not contribute to L(O) = i[H,O].
"""
function hamiltonian_terms(H::PauliSum)
    hz = Int128[]
    hx = Int128[]
    hny = Int[]
    hc = ComplexF64[]

    for (p, c_raw) in H
        c = ComplexF64(c_raw)

        iszero(c) && continue
        p.z == 0 && p.x == 0 && continue

        push!(hz, p.z)
        push!(hx, p.x)
        push!(hny, count_ones(p.z & p.x))
        push!(hc, c)
    end

    return (
        z = hz,
        x = hx,
        ny = hny,
        c = hc,
        nterms = length(hc),
        opnorm_est = 2.0 * sum(abs, hc),
    )
end

# ============================================================
# Matrix-free Liouvillian action as a closure, not a type
# ============================================================

"""
    make_liouvillian_action(H, tables)

Returns `(mulL!, hterms)`, where `mulL!(y, x)` applies y = L*x with
L(O) = i[H,O] in the Pauli coefficient-vector basis.

For one Hamiltonian term h_a P_a and basis string P_b,

    i[h_a P_a, P_b] = -2*s*h_a*P_ab

when P_a and P_b anticommute.
"""
function make_liouvillian_action(H::PauliSum{N}, tables = pauli_tables(Val(N))) where {N}
    hterms = hamiltonian_terms(H)

    function mulL!(y::AbstractVector, x::AbstractVector)
        @assert length(x) == tables.D
        @assert length(y) == tables.D

        fill!(y, 0.0 + 0.0im)

        @inbounds for a in eachindex(hterms.c)
            za = hterms.z[a]
            xa = hterms.x[a]
            na = hterms.ny[a]
            ca = hterms.c[a]

            for col in 1:tables.D
                amp = x[col]
                iszero(amp) && continue

                zb = tables.z[col]
                xb = tables.x[col]

                m1 = count_ones(xa & zb)
                m2 = count_ones(za & xb)

                # Pauli strings commute iff this parity is even.
                isodd(m1 - m2) || continue

                zprod = xor(za, zb)
                xprod = xor(xa, xb)

                nb = tables.ny[col]
                nab = count_ones(zprod & xprod)

                # Pauli product phase:
                # P_a P_b = i^k P_prod
                k = mod(nab - na - nb + 2 * m1, 4)

                # For anticommuting strings, k should be 1 or 3.
                # k = 1 -> +i, k = 3 -> -i.
                sign = 2 - k

                row = pauli_index_from_bits(zprod, xprod, Val(N))

                # L = i[H, .]
                # [P_a, P_b] = 2im*sign*P_prod
                # i*h_a*[P_a, P_b] = -2*sign*h_a*P_prod
                y[row] += (-2.0 * sign) * ca * amp
            end
        end

        return y
    end

    return mulL!, hterms
end

# ============================================================
# Sparse Liouvillian builder for validation / dense reference
# ============================================================

function liouvillian_sparse_fast(
    H::PauliSum{N},
    tables = pauli_tables(Val(N));
    atol::Float64 = 0.0,
) where {N}

    hterms = hamiltonian_terms(H)

    rows = Int[]
    cols = Int[]
    vals = ComplexF64[]

    @inbounds for a in eachindex(hterms.c)
        za = hterms.z[a]
        xa = hterms.x[a]
        na = hterms.ny[a]
        ca = hterms.c[a]

        for col in 1:tables.D
            zb = tables.z[col]
            xb = tables.x[col]

            m1 = count_ones(xa & zb)
            m2 = count_ones(za & xb)

            isodd(m1 - m2) || continue

            zprod = xor(za, zb)
            xprod = xor(xa, xb)

            nb = tables.ny[col]
            nab = count_ones(zprod & xprod)

            k = mod(nab - na - nb + 2 * m1, 4)
            sign = 2 - k

            row = pauli_index_from_bits(zprod, xprod, Val(N))
            val = (-2.0 * sign) * ca

            abs(val) > atol || continue

            push!(rows, row)
            push!(cols, col)
            push!(vals, val)
        end
    end

    return sparse(rows, cols, vals, tables.D, tables.D)
end

# Backward-compatible name from earlier versions.
function liouvillian_matrix(H::PauliSum{N}; atol::Float64 = 0.0) where {N}
    tables = pauli_tables(Val(N))
    return liouvillian_sparse_fast(H, tables; atol = atol)
end

# ============================================================
# Non-Hermitian Arnoldi exp(dt * L) * v
# ============================================================

"""
    arnoldi_exp_step(mulA!, v, cfg)

Approximate exp(cfg.dt * A) * v using an m-dimensional Arnoldi Krylov basis.
This is matrix-exponential-times-vector propagation without forming exp(A).

This removes the dependency on ExponentialUtilities and works for non-Hermitian
Liouvillians.
"""
function arnoldi_exp_step(mulA!, v::AbstractVector, cfg)
    n = length(v)
    beta = norm(v)

    if beta == 0
        return zeros(ComplexF64, n)
    end

    mmax = min(cfg.krylov_m, n)

    V = zeros(ComplexF64, n, mmax + 1)
    Hsmall = zeros(ComplexF64, mmax + 1, mmax)
    w = zeros(ComplexF64, n)

    @views V[:, 1] .= v ./ beta

    m_actual = 0
    happy_breakdown = false

    for j in 1:mmax
        @views mulA!(w, V[:, j])

        # Modified Gram-Schmidt.
        for i in 1:j
            @views hij = dot(V[:, i], w)
            Hsmall[i, j] += hij
            @views w .-= hij .* V[:, i]
        end

        # Optional second pass improves stability for non-Hermitian problems.
        if cfg.reorthogonalize
            for i in 1:j
                @views hcorr = dot(V[:, i], w)
                Hsmall[i, j] += hcorr
                @views w .-= hcorr .* V[:, i]
            end
        end

        hnext = norm(w)
        Hsmall[j + 1, j] = hnext
        m_actual = j

        if hnext <= cfg.krylov_breakdown_tol
            happy_breakdown = true
            break
        end

        if j < mmax
            @views V[:, j + 1] .= w ./ hnext
        end
    end

    k = m_actual
    Hk = @views Hsmall[1:k, 1:k]

    e1 = zeros(ComplexF64, k)
    e1[1] = beta

    ysmall = exp(cfg.dt * Matrix(Hk)) * e1

    y = zeros(ComplexF64, n)
    @views mul!(y, V[:, 1:k], ysmall)

    return y
end

function arnoldi_exp_step_sparse(A::SparseMatrixCSC, v::AbstractVector, cfg)
    work = zeros(ComplexF64, length(v))
    mulA! = (y, x) -> mul!(y, A, x)
    return arnoldi_exp_step(mulA!, v, cfg)
end

# ============================================================
# Validation helpers
# ============================================================

"""
    validate_liouvillian_action(H; nsamples=10, atol=1e-10)

Checks the matrix-free Liouvillian action against the reference

    coeff_vector(i * commutator(H, P))

for selected Pauli basis vectors.
"""
function validate_liouvillian_action(
    H::PauliSum{N};
    nsamples::Int = 10,
    atol::Float64 = 1e-10,
) where {N}

    tables = pauli_tables(Val(N))
    mulL!, hterms = make_liouvillian_action(H, tables)

    D = tables.D
    ncheck = min(nsamples, D)

    maxerr = 0.0

    for col in round.(Int, range(1, D; length = ncheck))
        p = PauliBasis{N}(tables.z[col], tables.x[col])

        v = zeros(ComplexF64, D)
        y = zeros(ComplexF64, D)

        v[col] = 1.0 + 0.0im
        mulL!(y, v)

        ps = singleton_paulisum(p)
        ref_ps = 1im * commutator(H, ps)
        ref = coeff_vector(ref_ps)

        err = norm(y - ref)
        maxerr = max(maxerr, err)

        if err > atol
            @warn "Liouvillian validation failed" col err atol
        end
    end

    @printf("Max Liouvillian action validation error = %.6e\n", maxerr)

    return maxerr
end

"""
    validate_arnoldi_against_dense(O0, H, cfg)

Small-system check comparing one Arnoldi step against dense exp(dt * L) * v.
"""
function validate_arnoldi_against_dense(O0::PauliSum{N}, H::PauliSum{N}, cfg) where {N}
    tables = pauli_tables(Val(N))
    v0 = coeff_vector(O0; T = ComplexF64)

    L = liouvillian_sparse_fast(H, tables)
    dense_ref = exp(cfg.dt * Matrix(L)) * v0

    mulL!, hterms = make_liouvillian_action(H, tables)
    arn = arnoldi_exp_step(mulL!, v0, cfg)

    err_abs = norm(arn - dense_ref)
    err_rel = err_abs / max(norm(dense_ref), eps(Float64))

    @printf("Arnoldi vs dense one-step abs error = %.6e\n", err_abs)
    @printf("Arnoldi vs dense one-step rel error = %.6e\n", err_rel)

    return (
        abs_error = err_abs,
        rel_error = err_rel,
    )
end

# ============================================================
# Propagation and Metrics
# ============================================================

"""
Exact-in-time Liouville propagation, with a matrix-exponential-times-vector
Krylov approximation at each saved time step.

Available propagation modes:
    :arnoldi_matrixfree
        Does not build L as a matrix. Recommended for larger N.

    :arnoldi_sparse
        Builds sparse L, then applies Arnoldi expv. Useful for comparison.

    :dense_exp
        Old dense exp(dt * L) method. Only for small N and validation.
"""
function weight_dynamics(
    O0::PauliSum{N,TO},
    H::PauliSum{N,TH},
    cfg,
) where {N,TO,TH}

    snapshots = Vector{Vector{Float64}}()
    raw_norms = Float64[]

    tables = pauli_tables(Val(N))
    D = tables.D

    @printf("Liouville dimension: %d\n", D)
    @printf("Propagation mode: %s\n", string(cfg.propagation))

    vt = coeff_vector(O0; T = ComplexF64)

    w0, n0 = weight_profile(
        vt,
        tables;
        normalize = cfg.normalize_weights,
        coeff_threshold = 0.0,
    )

    push!(snapshots, w0)
    push!(raw_norms, n0)

    if cfg.propagation === :arnoldi_matrixfree
        println("Building matrix-free Liouvillian action")
        mulL!, hterms = make_liouvillian_action(H, tables)
        @printf("Hamiltonian terms used in Liouvillian: %d\n", hterms.nterms)
        @printf("Estimated ||L|| bound: %.6e\n", hterms.opnorm_est)
        @printf("Arnoldi dimension m = %d\n", cfg.krylov_m)

        for step in 1:cfg.n_steps
            vt = arnoldi_exp_step(mulL!, vt, cfg)

            w, norm2 = weight_profile(
                vt,
                tables;
                normalize = cfg.normalize_weights,
                coeff_threshold = cfg.coeff_threshold,
            )

            push!(snapshots, w)
            push!(raw_norms, norm2)
        end

        tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

        return (
            tgrid = tgrid,
            snapshots = snapshots,
            raw_norms = raw_norms,
            metrics = time_series_metrics(snapshots),
            final_coeffs = vt,
            propagation = cfg.propagation,
        )

    elseif cfg.propagation === :arnoldi_sparse
        println("Building sparse Liouvillian matrix")
        L = liouvillian_sparse_fast(H, tables)
        @printf("Sparse L nnz: %d\n", nnz(L))
        @printf("Arnoldi dimension m = %d\n", cfg.krylov_m)

        for step in 1:cfg.n_steps
            vt = arnoldi_exp_step_sparse(L, vt, cfg)

            w, norm2 = weight_profile(
                vt,
                tables;
                normalize = cfg.normalize_weights,
                coeff_threshold = cfg.coeff_threshold,
            )

            push!(snapshots, w)
            push!(raw_norms, norm2)
        end

        tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

        return (
            tgrid = tgrid,
            snapshots = snapshots,
            raw_norms = raw_norms,
            metrics = time_series_metrics(snapshots),
            final_coeffs = vt,
            propagation = cfg.propagation,
            liouvillian = L,
        )

    elseif cfg.propagation === :dense_exp
        println("Building sparse Liouvillian, then dense timestep propagator")
        L = liouvillian_sparse_fast(H, tables)

        @printf("Sparse L nnz: %d\n", nnz(L))
        println("Exponentiating dense one-step propagator")
        Udt = exp(cfg.dt * Matrix(L))

        tmp = similar(vt)

        for step in 1:cfg.n_steps
            mul!(tmp, Udt, vt)
            vt, tmp = tmp, vt

            w, norm2 = weight_profile(
                vt,
                tables;
                normalize = cfg.normalize_weights,
                coeff_threshold = cfg.coeff_threshold,
            )

            push!(snapshots, w)
            push!(raw_norms, norm2)
        end

        tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

        return (
            tgrid = tgrid,
            snapshots = snapshots,
            raw_norms = raw_norms,
            metrics = time_series_metrics(snapshots),
            final_coeffs = vt,
            propagation = cfg.propagation,
            liouvillian = L,
            propagator = Udt,
        )

    else
        error("Unknown cfg.propagation = $(cfg.propagation). Use :arnoldi_matrixfree, :arnoldi_sparse, or :dense_exp.")
    end
end

# ============================================================
# Running functions
# ============================================================

"""
Run a single point in the PES.

This function propagates a Z_i operator, where the i-th qubit is z_qubit.
"""
function run_single_distance(
    R::Float64,
    data_path::String,
    cfg;
    z_qubit::Int = 1,
    block::Bool = false,
    NOI::Bool = false,
    UNRESTRICTED::Bool = false,
    run_validation::Bool = false,
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
            data_path;
            NOI = NOI,
            block = block,
        )
    else
        println("Building RESTRICTED Hamiltonian")
        H = QuantumChemQC.molecular_hamiltonian(
            Norbs,
            data_path;
            NOI = NOI,
            block = block,
        )
    end

    O = PauliSum(Pauli(Nqubits, Z = [z_qubit]))

    if run_validation
        println("Validating matrix-free Liouvillian action")
        validate_liouvillian_action(H; nsamples = min(16, Int(4)^Nqubits))
        println("Validating one Arnoldi step against dense exp reference")
        validate_arnoldi_against_dense(O, H, cfg)
    end

    out = weight_dynamics(O, H, cfg)

    @printf("Initial norm^2 = %.12e\n", out.raw_norms[1])
    @printf("Final norm^2   = %.12e\n", out.raw_norms[end])
    @printf("Relative norm^2 = %.12e\n", out.raw_norms[end] / out.raw_norms[1])
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
function run_distance_scan(
    cases,
    cfg;
    z_qubit::Int = 1,
    block::Bool = false,
    NOI::Bool = false,
    UNRESTRICTED::Bool = false,
    run_validation::Bool = false,
)

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
            run_validation = run_validation,
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

    return plot(p_entropy)
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

BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H2-PP/tensors"

CASES = [
    (0.50, joinpath(BASE, "RHF_H2_R_0p5_tensors.npz")),
    (0.7414, joinpath(BASE, "RHF_H2_R_0p7414_tensors.npz")),
    (1.00, joinpath(BASE, "RHF_H2_R_1p0_tensors.npz")),
    (1.50, joinpath(BASE, "RHF_H2_R_1p5_tensors.npz")),
    (2.00, joinpath(BASE, "RHF_H2_R_2p0_tensors.npz")),
    (2.50, joinpath(BASE, "RHF_H2_R_2p5_tensors.npz")),
    (3.00, joinpath(BASE, "RHF_H2_R_3p0_tensors.npz")),
]

# ============================================================
# Main run
# ============================================================

cfg = make_config(
    n_steps = 500,
    dt = 0.1,
    coeff_threshold = 1e-12,
    normalize_weights = false,

    # Recommended larger-system path:
    propagation = :arnoldi_matrixfree,

    # Increase this to reduce Krylov truncation error.
    # For H2, 30-40 should usually be enough; for larger systems test norm drift.
    krylov_m = 40,
    krylov_breakdown_tol = 1e-13,
    reorthogonalize = true,
)

scan = run_distance_scan(
    CASES,
    cfg;
    z_qubit = 2,
    block = false,
    NOI = false,
    UNRESTRICTED = false,

    # Set true for small debugging runs. It builds dense references, so keep false
    # for larger systems.
    run_validation = false,
)

p_summary = plot_distance_scan_summary(scan)
display(p_summary)
savefig(p_summary, "majorana_distance_scan_summary.png")

p_heat = plot_selected_heatmaps(scan; selected_Rs = [0.50, 0.7414, 1.50, 3.00])
display(p_heat)
savefig(p_heat, "majorana_weight_heatmaps_selected_R.png")

results_scan_table(scan)