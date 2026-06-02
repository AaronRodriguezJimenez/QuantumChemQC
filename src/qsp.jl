#
# DMD and related QSP utilities.
#

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
