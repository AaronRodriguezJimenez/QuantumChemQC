#
# Majorana-weight dynamics using sparse-dictionary Liouville Arnoldi propagation.
#
# No module and no new top-level structs are defined, so this file is safe to
# re-include in a VS Code/REPL session.
#
# Main idea:
#   Keep O(t) as the PauliSum / Dict representation used by PauliOperators.jl.
#   Do not allocate a dense coefficient vector of length 4^N.
#   Do not scan the full Pauli basis at each Liouvillian application.
#
# Propagation:
#   O(t + dt) approx exp(dt * L) O(t),  L(O) = i[H, O]
# using a non-Hermitian Arnoldi Krylov method in sparse Pauli-dictionary space.
#
# This has no Trotter error. The controllable errors are Krylov truncation and,
# if enabled, coefficient pruning.
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

Important fields:
    n_steps, dt
    coeff_threshold
        Used only when recording the Majorana-weight profile.
    state_prune_threshold
        If > 0, coefficients smaller than this are deleted from O(t) after each
        time step and inside sparse Arnoldi vector algebra. This improves speed
        but introduces truncation error. Keep 0.0 for the least intrusive run.
    krylov_m
        Arnoldi subspace dimension.
    reorthogonalize
        Second Gram-Schmidt pass for stability.
"""
function make_config(;
    n_steps::Int,
    dt::Real,
    coeff_threshold::Real = 1e-12,
    normalize_weights::Bool = false,
    propagation::Symbol = :arnoldi_sparse_dict,
    krylov_m::Int = 30,
    krylov_breakdown_tol::Real = 1e-13,
    reorthogonalize::Bool = true,
    state_prune_threshold::Real = 0.0,
    arnoldi_drop_tol::Real = 0.0,
    progress_every::Int = 50,
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
        state_prune_threshold = Float64(state_prune_threshold),
        arnoldi_drop_tol = Float64(arnoldi_drop_tol),
        progress_every = progress_every,
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
# Pauli indexing and sparse PauliSum utilities
# ============================================================

@inline function pauli_linear_key(p::PauliBasis{N}) where {N}
    return p.z + (p.x << N) + Int128(1)
end

@inline function pauli_index(p::PauliBasis{N}) where {N}
    return Int(pauli_linear_key(p))
end

@inline function pauli_index_from_bits(z::Int128, x::Int128, ::Val{N}) where {N}
    return Int(z + (x << N) + Int128(1))
end

function complexify(ps::PauliSum{N}) where {N}
    out = PauliSum(N, ComplexF64)
    for (p, c) in ps
        out[p] = ComplexF64(c)
    end
    return out
end

function sparse_norm2(ps::PauliSum)
    s = 0.0
    for (_, c) in ps
        s += abs2(c)
    end
    return s
end

function sparse_norm(ps::PauliSum)
    return sqrt(sparse_norm2(ps))
end

"""
    sparse_dot(a, b)

Hilbert-Schmidt coefficient-space dot product:
    sum_P conj(a_P) * b_P
Iterates over the smaller dictionary.
"""
function sparse_dot(a::PauliSum, b::PauliSum)
    if length(a) <= length(b)
        s = 0.0 + 0.0im
        for (p, ca) in a
            s += conj(ca) * get(b, p, 0.0 + 0.0im)
        end
        return s
    else
        s = 0.0 + 0.0im
        for (p, cb) in b
            s += conj(get(a, p, 0.0 + 0.0im)) * cb
        end
        return s
    end
end

function sparse_scale_copy(ps::PauliSum{N}, alpha::Number; drop_tol::Float64 = 0.0) where {N}
    out = PauliSum(N, ComplexF64)
    alpha_c = ComplexF64(alpha)

    for (p, c) in ps
        v = alpha_c * ComplexF64(c)
        if drop_tol <= 0.0 || abs(v) > drop_tol
            out[p] = v
        end
    end

    return out
end

"""
    sparse_axpy!(y, alpha, x; drop_tol=0.0)

Sparse dictionary y <- y + alpha*x.
If drop_tol > 0, entries with |coefficient| <= drop_tol are deleted.
"""
function sparse_axpy!(y::PauliSum{N}, alpha::Number, x::PauliSum{N}; drop_tol::Float64 = 0.0) where {N}
    alpha_c = ComplexF64(alpha)

    for (p, cx_raw) in x
        cx = ComplexF64(cx_raw)
        newval = get(y, p, 0.0 + 0.0im) + alpha_c * cx

        if drop_tol > 0.0 && abs(newval) <= drop_tol
            if haskey(y, p)
                delete!(y, p)
            end
        else
            y[p] = newval
        end
    end

    return y
end

function cleanup_paulisum!(ps::PauliSum{N}; thresh::Float64 = 0.0) where {N}
    thresh <= 0.0 && return ps

    to_delete = Vector{PauliBasis{N}}()
    for (p, c) in ps
        if abs(c) <= thresh
            push!(to_delete, p)
        end
    end
    for p in to_delete
        delete!(ps, p)
    end

    return ps
end

function ps_difference_norm(a::PauliSum{N}, b::PauliSum{N}) where {N}
    tmp = complexify(a)
    sparse_axpy!(tmp, -1.0, complexify(b))
    return sparse_norm(tmp)
end

# Optional dense coefficient-vector utility for small validation only.
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
# Majorana-weight profile from active dictionary support
# ============================================================

@inline function cached_majorana_weight!(cache::Dict{Int128,Int}, p::PauliBasis)
    key = pauli_linear_key(p)
    return get!(cache, key) do
        QuantumChemQC.majorana_weight(p)
    end
end

"""
    weight_profile(Ot, Val(N), weight_cache; normalize=true, coeff_threshold=0.0)

Builds
    p_k(t) = sum_{P: majorana_weight(P)=k} |c_P(t)|^2
by iterating only over the Pauli strings present in the dictionary for O(t).

This is the fast path requested: it never scans all 4^N Pauli basis elements.
"""
function weight_profile(
    Ot::PauliSum{N},
    ::Val{N},
    weight_cache::Dict{Int128,Int};
    normalize::Bool = true,
    coeff_threshold::Float64 = 0.0,
) where {N}

    hist = zeros(Float64, 2 * N + 1)
    total = 0.0

    for (p, c) in Ot
        abs(c) > coeff_threshold || continue

        k = cached_majorana_weight!(weight_cache, p)
        a2 = abs2(c)

        hist[k + 1] += a2
        total += a2
    end

    if normalize && total > 0
        hist ./= total
    end

    return hist, total
end

# Compatibility helper if no cache is supplied.
function weight_profile(
    Ot::PauliSum{N},
    ::Val{N};
    normalize::Bool = true,
    coeff_threshold::Float64 = 0.0,
) where {N}
    return weight_profile(
        Ot,
        Val(N),
        Dict{Int128,Int}();
        normalize = normalize,
        coeff_threshold = coeff_threshold,
    )
end

# ============================================================
# Hamiltonian preprocessing
# ============================================================

"""
    hamiltonian_terms(H)

Preprocess H into z/x/coefficient arrays. Identity terms are skipped because
identity commutes with every operator and does not contribute to L(O)=i[H,O].
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
# Sparse dictionary Liouvillian action
# ============================================================

"""
    liouvillian_action_sparse(O, hterms; drop_tol=0.0)

Apply L(O)=i[H,O] in Pauli dictionary space.

For one Hamiltonian term h_a P_a and one operator term c_b P_b,
if P_a and P_b anticommute,

    i[h_a P_a, c_b P_b] = -2*s*h_a*c_b*P_ab,

where s = +1 or -1 is obtained from the z/x Pauli product phase.
"""
function liouvillian_action_sparse(
    O::PauliSum{N},
    hterms;
    drop_tol::Float64 = 0.0,
) where {N}

    out = PauliSum(N, ComplexF64)

    for (pb, amp_raw) in O
        amp = ComplexF64(amp_raw)
        iszero(amp) && continue
        if drop_tol > 0.0 && abs(amp) <= drop_tol
            continue
        end

        zb = pb.z
        xb = pb.x
        nb = count_ones(zb & xb)

        @inbounds for a in eachindex(hterms.c)
            za = hterms.z[a]
            xa = hterms.x[a]
            na = hterms.ny[a]
            ca = hterms.c[a]

            m1 = count_ones(xa & zb)
            m2 = count_ones(za & xb)

            # Pauli strings commute iff this parity is even.
            isodd(m1 - m2) || continue

            zprod = xor(za, zb)
            xprod = xor(xa, xb)
            nab = count_ones(zprod & xprod)

            # Pauli product phase:
            # P_a P_b = i^k P_prod
            k = mod(nab - na - nb + 2 * m1, 4)

            # For anticommuting strings, k should be 1 or 3.
            # k = 1 -> +i, k = 3 -> -i.
            sign = 2 - k

            pprod = PauliBasis{N}(zprod, xprod)
            val = (-2.0 * sign) * ca * amp

            newval = get(out, pprod, 0.0 + 0.0im) + val
            if drop_tol > 0.0 && abs(newval) <= drop_tol
                if haskey(out, pprod)
                    delete!(out, pprod)
                end
            else
                out[pprod] = newval
            end
        end
    end

    return out
end

# ============================================================
# Sparse Arnoldi exp(dt * L) * O
# ============================================================

"""
    arnoldi_exp_step_sparse_operator(O, hterms, cfg)

Approximate exp(cfg.dt * L)O using Arnoldi in sparse Pauli-dictionary space.
This avoids dense coefficient vectors and avoids scanning the full 4^N basis.
"""
function arnoldi_exp_step_sparse_operator(
    O::PauliSum{N},
    hterms,
    cfg,
) where {N}

    beta = sparse_norm(O)
    if beta == 0.0
        return PauliSum(N, ComplexF64)
    end

    mmax = cfg.krylov_m
    V = Vector{typeof(complexify(O))}(undef, mmax + 1)
    Hsmall = zeros(ComplexF64, mmax + 1, mmax)

    V[1] = sparse_scale_copy(O, 1.0 / beta; drop_tol = cfg.arnoldi_drop_tol)

    m_actual = 0

    for j in 1:mmax
        w = liouvillian_action_sparse(
            V[j],
            hterms;
            drop_tol = cfg.arnoldi_drop_tol,
        )

        # Modified Gram-Schmidt.
        for i in 1:j
            hij = sparse_dot(V[i], w)
            Hsmall[i, j] += hij
            sparse_axpy!(w, -hij, V[i]; drop_tol = cfg.arnoldi_drop_tol)
        end

        # Optional second pass for non-Hermitian stability.
        if cfg.reorthogonalize
            for i in 1:j
                hcorr = sparse_dot(V[i], w)
                Hsmall[i, j] += hcorr
                sparse_axpy!(w, -hcorr, V[i]; drop_tol = cfg.arnoldi_drop_tol)
            end
        end

        hnext = sparse_norm(w)
        Hsmall[j + 1, j] = hnext
        m_actual = j

        if hnext <= cfg.krylov_breakdown_tol
            break
        end

        if j < mmax
            V[j + 1] = sparse_scale_copy(w, 1.0 / hnext; drop_tol = cfg.arnoldi_drop_tol)
        end
    end

    k = m_actual
    Hk = @views Hsmall[1:k, 1:k]

    e1 = zeros(ComplexF64, k)
    e1[1] = beta

    ysmall = exp(cfg.dt * Matrix(Hk)) * e1

    out = PauliSum(N, ComplexF64)
    for i in 1:k
        sparse_axpy!(out, ysmall[i], V[i]; drop_tol = cfg.arnoldi_drop_tol)
    end

    cleanup_paulisum!(out; thresh = cfg.state_prune_threshold)
    return out
end

# ============================================================
# Optional dense/sparse-matrix reference builders for small validation only
# ============================================================

function pauli_tables(::Val{N}) where {N}
    D = Int(4)^N

    zs = Vector{Int128}(undef, D)
    xs = Vector{Int128}(undef, D)
    nys = Vector{Int}(undef, D)

    for p in PauliBasis{N}
        j = pauli_index(p)
        zs[j] = p.z
        xs[j] = p.x
        nys[j] = count_ones(p.z & p.x)
    end

    return (
        N = N,
        D = D,
        z = zs,
        x = xs,
        ny = nys,
    )
end

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
            nb = tables.ny[col]

            m1 = count_ones(xa & zb)
            m2 = count_ones(za & xb)

            isodd(m1 - m2) || continue

            zprod = xor(za, zb)
            xprod = xor(xa, xb)
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

function liouvillian_matrix(H::PauliSum{N}; atol::Float64 = 0.0) where {N}
    tables = pauli_tables(Val(N))
    return liouvillian_sparse_fast(H, tables; atol = atol)
end

function paulisum_from_coeffs(v::AbstractVector, ::Val{N}; atol::Float64 = 1e-14) where {N}
    out = PauliSum(N, ComplexF64)
    for p in PauliBasis{N}
        c = v[pauli_index(p)]
        if abs(c) > atol
            out[p] = ComplexF64(c)
        end
    end
    return out
end

# ============================================================
# Validation helpers
# ============================================================

function validate_liouvillian_action_sparse(
    H::PauliSum{N};
    nsamples::Int = 10,
    atol::Float64 = 1e-10,
) where {N}

    hterms = hamiltonian_terms(H)
    D = Int(4)^N
    ncheck = min(nsamples, D)
    maxerr = 0.0

    for col in round.(Int, range(1, D; length = ncheck))
        # Small validation only: constructs Pauli basis from dense ordering.
        k0 = Int128(col - 1)
        mask = (Int128(1) << N) - 1
        z = k0 & mask
        x = k0 >> N
        p = PauliBasis{N}(z, x)

        ps = singleton_paulisum(p)
        got = liouvillian_action_sparse(ps, hterms)
        ref = 1im * commutator(H, ps)

        err = ps_difference_norm(got, ref)
        maxerr = max(maxerr, err)

        if err > atol
            @warn "Sparse Liouvillian validation failed" col err atol
        end
    end

    @printf("Max sparse Liouvillian action validation error = %.6e\n", maxerr)
    return maxerr
end

function validate_arnoldi_against_dense(O0::PauliSum{N}, H::PauliSum{N}, cfg) where {N}
    hterms = hamiltonian_terms(H)

    L = liouvillian_sparse_fast(H)
    v0 = coeff_vector(O0; T = ComplexF64)
    dense_vec = exp(cfg.dt * Matrix(L)) * v0
    dense_ps = paulisum_from_coeffs(dense_vec, Val(N); atol = 0.0)

    arn_ps = arnoldi_exp_step_sparse_operator(complexify(O0), hterms, cfg)

    err_abs = ps_difference_norm(arn_ps, dense_ps)
    err_rel = err_abs / max(sparse_norm(dense_ps), eps(Float64))

    @printf("Arnoldi sparse-dict vs dense one-step abs error = %.6e\n", err_abs)
    @printf("Arnoldi sparse-dict vs dense one-step rel error = %.6e\n", err_rel)

    return (
        abs_error = err_abs,
        rel_error = err_rel,
    )
end

# ============================================================
# Propagation and metrics
# ============================================================

"""
    weight_dynamics(O0, H, cfg)

Sparse-dictionary Liouville propagation.

The state O(t) remains a PauliSum/Dict at all times. Majorana-weight profiles
are recorded by iterating over the active dictionary keys only.
"""
function weight_dynamics(
    O0::PauliSum{N,TO},
    H::PauliSum{N,TH},
    cfg,
) where {N,TO,TH}

    snapshots = Vector{Vector{Float64}}()
    raw_norms = Float64[]
    nnz_terms = Int[]

    @printf("Sparse-dict Liouville space, Nqubits = %d\n", N)
    @printf("Initial active Pauli terms: %d\n", length(O0))
    @printf("Propagation mode: %s\n", string(cfg.propagation))

    if cfg.propagation !== :arnoldi_sparse_dict
        error("This optimized script uses cfg.propagation = :arnoldi_sparse_dict")
    end

    hterms = hamiltonian_terms(H)
    @printf("Hamiltonian terms used in Liouvillian: %d\n", hterms.nterms)
    @printf("Estimated ||L|| bound: %.6e\n", hterms.opnorm_est)
    @printf("Arnoldi dimension m = %d\n", cfg.krylov_m)
    @printf("State prune threshold = %.3e\n", cfg.state_prune_threshold)
    @printf("Arnoldi drop tolerance = %.3e\n", cfg.arnoldi_drop_tol)

    Ot = complexify(O0)
    weight_cache = Dict{Int128,Int}()

    w0, n0 = weight_profile(
        Ot,
        Val(N),
        weight_cache;
        normalize = cfg.normalize_weights,
        coeff_threshold = 0.0,
    )
    push!(snapshots, w0)
    push!(raw_norms, n0)
    push!(nnz_terms, length(Ot))

    for step in 1:cfg.n_steps
        Ot = arnoldi_exp_step_sparse_operator(Ot, hterms, cfg)

        if cfg.state_prune_threshold > 0.0
            cleanup_paulisum!(Ot; thresh = cfg.state_prune_threshold)
        end

        w, norm2 = weight_profile(
            Ot,
            Val(N),
            weight_cache;
            normalize = cfg.normalize_weights,
            coeff_threshold = cfg.coeff_threshold,
        )

        push!(snapshots, w)
        push!(raw_norms, norm2)
        push!(nnz_terms, length(Ot))

        if cfg.progress_every > 0 && (step % cfg.progress_every == 0 || step == cfg.n_steps)
            @printf(
                "step %5d / %5d   active_terms = %9d   norm^2 = %.8e\n",
                step,
                cfg.n_steps,
                length(Ot),
                norm2,
            )
        end
    end

    tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

    return (
        tgrid = tgrid,
        snapshots = snapshots,
        raw_norms = raw_norms,
        metrics = time_series_metrics(snapshots),
        nnz_terms = nnz_terms,
        final_operator = Ot,
        weight_cache = weight_cache,
        propagation = cfg.propagation,
    )
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
        println("Validating sparse dictionary Liouvillian action")
        validate_liouvillian_action_sparse(H; nsamples = min(16, Int(4)^Nqubits))
        println("Validating one sparse Arnoldi step against dense exp reference")
        validate_arnoldi_against_dense(O, H, cfg)
    end

    out = weight_dynamics(O, H, cfg)

    @printf("Initial norm^2 = %.12e\n", out.raw_norms[1])
    @printf("Final norm^2   = %.12e\n", out.raw_norms[end])
    @printf("Relative norm^2 = %.12e\n", out.raw_norms[end] / out.raw_norms[1])
    @printf("Initial entropy = %.8f\n", out.metrics.entropy[1])
    @printf("Final entropy   = %.8f\n", out.metrics.entropy[end])
    @printf("Final active Pauli terms = %d\n", out.nnz_terms[end])

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
# Plots and tables
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

function plot_active_terms(scan)
    Rs = sorted_Rs(scan)

    p = plot(
        xlabel = "time",
        ylabel = "active Pauli terms",
        title = "Sparse dictionary support size",
        lw = 2,
        yscale = :log10,
    )

    for R in Rs
        out = scan[R].out
        label = @sprintf("R=%.4g", R)
        plot!(p, out.tgrid, out.nnz_terms; label = label, lw = 2)
    end

    return p
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
    println("R        final_S       max_S        norm_ret       final_terms")
    println("-"^105)

    for R in Rs
        out = scan[R].out
        final_S = out.metrics.entropy[end]
        max_S = maximum(out.metrics.entropy)
        norm_ret = out.raw_norms[end] / out.raw_norms[1]
        final_terms = out.nnz_terms[end]

        @printf(
            "%-8.4f %-13.6f %-13.6f %-14.8f %-12d\n",
            R,
            final_S,
            max_S,
            norm_ret,
            final_terms,
        )
    end
end

# ============================================================
# DATA INPUT AND PARAMETERS
# ============================================================

const BASE = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/MEAO-tensors" # MEAO natural orbitals
const CASES = [
    (0.50, joinpath(BASE, "h4_MEAO_LIN4_R_0p5_tensors.npz"))
    (0.7414, joinpath(BASE, "h4_MEAO_LIN4_R_0p7414_tensors.npz"))
    (1.00, joinpath(BASE, "h4_MEAO_LIN4_R_1_tensors.npz"))
    (1.50, joinpath(BASE, "h4_MEAO_LIN4_R_1p5_tensors.npz"))
    (2.00, joinpath(BASE, "h4_MEAO_LIN4_R_2_tensors.npz"))
    (2.50, joinpath(BASE, "h4_MEAO_LIN4_R_2p5_tensors.npz"))
    (3.00, joinpath(BASE, "h4_MEAO_LIN4_R_3_tensors.npz"))
]

# ============================================================
# Main run
# ============================================================

cfg = make_config(
    n_steps = 500,
    dt = 0.1,
    coeff_threshold = 1e-12,
    normalize_weights = false,

    propagation = :arnoldi_sparse_dict,

    # Start with 20-30 for speed, then compare norm drift / dense validation on
    # small systems. Larger m improves accuracy but costs more.
    krylov_m = 25,
    krylov_breakdown_tol = 1e-13,
    reorthogonalize = true,

    # Keep these at 0.0 for the least intrusive run. Increase to 1e-14 or 1e-12
    # only if support growth becomes too large and you accept truncation error.
    state_prune_threshold = 0.0,
    arnoldi_drop_tol = 0.0,

    progress_every = 50,
)

scan = run_distance_scan(
    CASES,
    cfg;
    z_qubit = 4,
    block = false,
    NOI = false,
    UNRESTRICTED = false,

    # Set true only for small debugging runs. It builds dense references.
    run_validation = false,
)

p_summary = plot_distance_scan_summary(scan)
display(p_summary)
savefig(p_summary, "majorana_distance_scan_summary.png")

p_terms = plot_active_terms(scan)
display(p_terms)
savefig(p_terms, "majorana_active_terms.png")

p_heat = plot_selected_heatmaps(scan; selected_Rs = [0.50, 0.7414, 1.50, 3.00])
display(p_heat)
savefig(p_heat, "majorana_weight_heatmaps_selected_R.png")

results_scan_table(scan)
