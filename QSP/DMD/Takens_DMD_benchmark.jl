# Takens + DMD comparison benchmark
#
# Goals:
#   - Compare DMD+Takens behavior on structured vs random signals.
#   - Use the same analysis pipeline across all cases.
#   - Check whether PP operator-weight dynamics behave more like structured
#     dynamics or like random / shuffled controls.
#
# Assumes the PP functions/types already exist in the session:
#   - evolution_op
#   - heisenberg_1D
#   - fit_dmd
#   - weight_stats
#   - plot_weight_heatmap
#   - plot_dmd_modes_abs
#
# This script is intentionally independent from the production PP code.
using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Statistics
using Printf
using Random
using Plots
using NPZ
#"""
#The workflow is:
#
#- Build a Heisenberg Hamiltonian as a PauliSum.
#- Start from an initial Pauli operator O.
#- Evolve that operator under the Hamiltonian using repeated conjugations by Pauli rotations.
#- After each Trotter layer, record the operator’s weight profile: 
#  how much coefficient norm lives on Pauli strings of weight 0, 1, …, N.
#- Assemble those weight profiles into a snapshot matrix.
#- Run DMD on that matrix to extract dominant modes, eigenvalues, and a low-rank linear map.
#- Plot the resulting weight dynamics and DMD modes.
# ------------------------------------------------------------
#
#"""

# ------------------------------------------------------------
neighbor(site::Int, N::Int) = site == N ? 1 : site + 1

# ------------------------------------------------------------
# Weight diagnostics
# ------------------------------------------------------------
function weight_profile(W::PauliSum{N, T}) where {N, T}
    #hist = zeros(Float64, N + 1)
    hist = zeros(Float64, 2*N + 1)
    total = 0.0

    for (P, c) in W
        #k = weight(P)
        k = QuantumChemQC.majorana_weight(P)
        a2 = abs2(c)
        hist[k + 1] += a2
        total += a2
    end
    
    # Normalization
#    if total > 0
#        hist ./= total
#    end

    return hist
end

function weight_stats(w::AbstractVector)
    ks = 0:length(w)-1
    μ = sum(ks .* w)
    σ2 = sum(((ks .- μ).^2) .* w)
    return μ, σ2
end

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

function dmd_bootstrap(snapshots::Vector{<:AbstractVector}; B::Int=100, r::Union{Nothing,Int}=nothing, seed::Int=1)
    rng = MersenneTwister(seed)
    m = length(snapshots)
    evals = Vector{Vector{ComplexF64}}(undef, B)

    for b in 1:B
        lo = rand(rng, 1:m-3)
        hi = rand(rng, lo+2:m)
        Xb = hcat(snapshots[lo:hi]...)
        evals[b] = QuantumChemQC.fit_dmd(Xb; r=r).evals
    end

    return evals
end

# ------------------------------------------------------------
# Main evolution loop
# ------------------------------------------------------------

function evolution_op(ket, o::PauliSum{N, T}, H::PauliSum{N, T}, n_intervals, dt;
    thresh::Float64=1e-3) where {N, T}

    O0 = deepcopy(o)
    Ot = deepcopy(o)

    corr_real = Float64[]
    corr_imag = Float64[]
    snapshots = Vector{Vector{Float64}}()

    push!(snapshots, weight_profile(Ot))

    c0 = expectation_value(O0 * Ot, ket)
    push!(corr_real, real(c0))
    push!(corr_imag, imag(c0))

    generators, angles = QuantumChemQC.gens_from_H(H)
    nt = length(angles)

    println("Total Rotations: ", nt * n_intervals)

    for _ in 1:n_intervals
        accumulated_error = 0.0 + 0.0im

        for j in 1:nt
            θ = 2 * dt * angles[j]
            QuantumChemQC.evolve!(Ot, generators[j], θ)

            QuantumChemQC.coeff_clip!(Ot; thresh=1e-12)
            before = expectation_value(O0 * Ot, ket)

            QuantumChemQC.coeff_clip!(Ot; thresh=thresh)
            after = expectation_value(O0 * Ot, ket)

            accumulated_error += after - before
        end

        push!(snapshots, weight_profile(Ot))

        c = expectation_value(O0 * Ot, ket) + accumulated_error
        push!(corr_real, real(c))
        push!(corr_imag, imag(c))
    end

    tgrid = collect(range(0.0, stop=n_intervals * dt, length=length(corr_real)))
    return corr_real, corr_imag, tgrid, snapshots
end

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

function plot_weight_heatmap(w_snapshots; tgrid=nothing)
    W = hcat(w_snapshots...)
    nweights, nt = size(W)

    if tgrid === nothing
        tgrid = 0:nt-1
    end

    heatmap(
        tgrid,
        0:nweights-1,
        W,
        xlabel = "time",
        ylabel = "Weight",
        title = "L2 norm Weight Dynamics",
        #legend = false
    )
end

function plot_weight_stack(w_snapshots; tgrid=nothing)
    W = hcat(w_snapshots...)
    nweights, nt = size(W)

    if tgrid === nothing
        tgrid = 0:nt-1
    end

    plt = plot(
        xlabel = "time",
        ylabel = "weight fraction",
        title = "Pauli Weight Distribution",
        legend = :right
    )

    for k in 1:nweights
        plot!(plt, tgrid, W[k, :], label = "k=$(k-1)", lw=2)
    end

    return plt
end

function plot_dmd_modes_abs(res::DMDResult; nmodes=4)
    r = min(nmodes, size(res.modes, 2))
    ks = 0:(size(res.modes, 1)-1)

    p = plot(
        xlabel = "Pauli weight",
        ylabel = "|mode amplitude|",
        title = "Magnitude of DMD modes",
        legend = :right
    )

    for j in 1:r
        plot!(p, ks, abs.(res.modes[:, j]), label = "mode $j", lw=2)
    end

    return p
end

# ------------------------------------------------------------
# Basic utilities
# ------------------------------------------------------------

mean_weight(w::AbstractVector) = sum((0:length(w)-1) .* w)

function weight_entropy(w::AbstractVector; ϵ::Float64=1e-15)
    p = clamp.(w, ϵ, 1.0)
    return -sum(p .* log.(p))
end

function dominant_rank(s::AbstractVector{<:Real}; energy::Float64=0.95)
    tot = sum(abs2, s)
    tot <= 0 && return 0
    acc = 0.0
    for (i, σ) in enumerate(s)
        acc += abs2(σ)
        if acc / tot >= energy
            return i
        end
    end
    return length(s)
end

function l2_relative(a::AbstractVector, b::AbstractVector)
    return norm(a .- b) / max(norm(b), eps())
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

    return (mean_weight = μ, variance = σ2, entropy = ent)
end

# ------------------------------------------------------------
# Signal generators
# ------------------------------------------------------------

function synthetic_signal(m::Int; ω1::Float64=0.12, ω2::Float64=0.31)
    t = collect(0:m-1)
    x1 = sin.(ω1 .* t) .+ 0.3 .* cos.(ω2 .* t)
    x2 = cos.(ω1 .* t) .- 0.2 .* sin.(ω2 .* t)
    x3 = 0.5 .* sin.(0.5 .* ω1 .* t .+ 0.1)
    return [Float64[x1[i], x2[i], x3[i]] for i in 1:m]
end

function white_noise_signal(m::Int, d::Int; seed::Int=1)
    rng = MersenneTwister(seed)
    return [randn(rng, d) for _ in 1:m]
end

"""
 Correlated random signal (AR process)

Pure white noise is almost too easy. 
A more realistic null model is correlated stochastic dynamics
AR(1) formula: x_t = α * x_{t-1} + noise
"""
function ar1_signal(m::Int, d::Int; α::Float64=0.9, seed::Int=1)
    rng = MersenneTwister(seed)
    x = zeros(Float64, d)
    out = Vector{Vector{Float64}}(undef, m)
    for t in 1:m
        x = α .* x .+ 0.1 .* randn(rng, d)
        out[t] = copy(x)
    end
    return out
end

"""
 Take PP snapshots and shuffle time order
"""
function shuffled_snapshots(snapshots::Vector{<:AbstractVector}; seed::Int=1)
    rng = MersenneTwister(seed)
    idx = randperm(rng, length(snapshots))
    return snapshots[idx]
end

# ------------------------------------------------------------
# Benchmark container
# ------------------------------------------------------------

struct BenchmarkResult
    name::String
    q_values::Vector{Int}
    fit_error::Vector{Float64}
    rank95::Vector{Int}
    rank99::Vector{Int}
    sv_first::Vector{Float64}
    sv_second::Vector{Float64}
end

function benchmark_signal(name::String, snapshots::Vector{<:AbstractVector}; dt = 0.1, q_values=1:10, r::Union{Nothing,Int}=2)
    qv = collect(q_values)
    fit_error = Float64[]
    rank95 = Int[]
    rank99 = Int[]
    sv_first = Float64[]
    sv_second = Float64[]
    @printf("\nBenchmarking %s...\n", name)
    for q in qv
        if length(snapshots) < q + 1
            @printf("Skipping q=%d for %s (not enough snapshots)\n", q, name)
            push!(fit_error, NaN)
            push!(rank95, 0)
            push!(rank99, 0)
            push!(sv_first, NaN)
            push!(sv_second, NaN)
            continue
        end

        out = QuantumChemQC.embed_fit_summary(snapshots; q=q, r=r)
        QuantumChemQC.print_dmd_summary(out.res; dt=dt, topk=5) # Summary for dominant topk modes
        
        push!(fit_error, out.rel_fit)
        push!(rank95, out.rank95)
        push!(rank99, out.rank99)
        push!(sv_first, isempty(out.res.singular_values) ? NaN : out.res.singular_values[1])
        push!(sv_second, length(out.res.singular_values) >= 2 ? out.res.singular_values[2] : NaN)
    end

    return BenchmarkResult(name, qv, fit_error, rank95, rank99, sv_first, sv_second)
end

# ------------------------------------------------------------
# PP-specific wrappers
# ------------------------------------------------------------

function run_pp_snapshots(ket, o, H, n_intervals::Int, dt::Real; threshold::Float64=1e-10)
    rRES, iRES, tgrid, w_snapshots = evolution_op(ket, o, H, n_intervals, dt; thresh=threshold)
    return (
        rRES = rRES,
        iRES = iRES,
        tgrid = tgrid,
        w_snapshots = w_snapshots,
        metrics = time_series_metrics(w_snapshots),
    )
end

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

function plot_benchmark_metric(results::Vector{BenchmarkResult}; metric::Symbol=:fit_error)
    p = plot(
        xlabel = "embedding q",
        ylabel = string(metric),
        title = "Takens + DMD benchmark",
        legend = :best,
        xscale = :linear,
    )

    for r in results
        vals = getproperty(r, metric)
        plot!(p, r.q_values, vals, marker=:circle, lw=2, label=r.name)
    end

    return p
end

function plot_benchmark_ranks(results::Vector{BenchmarkResult})
    p1 = plot(
        xlabel = "embedding q",
        ylabel = "rank for 95% energy",
        title = "Effective rank (95%)",
        legend = :best,
    )
    p2 = plot(
        xlabel = "embedding q",
        ylabel = "rank for 99% energy",
        title = "Effective rank (99%)",
        legend = :best,
    )

    for r in results
        plot!(p1, r.q_values, r.rank95, marker=:circle, lw=2, label=r.name)
        plot!(p2, r.q_values, r.rank99, marker=:circle, lw=2, label=r.name)
    end

    return p1, p2
end

function plot_singular_values_vs_q(results::Vector{BenchmarkResult})
    p = plot(
        xlabel = "embedding q",
        ylabel = "leading singular value",
        title = "Leading singular values vs embedding dimension",
        legend = :best,
    )

    for r in results
        plot!(p, r.q_values, r.sv_first, marker=:circle, lw=2, label="$(r.name): σ₁")
        plot!(p, r.q_values, r.sv_second, marker=:square, lw=2, label="$(r.name): σ₂")
    end

    return p
end

# ------------------------------------------------------------
# Main benchmark
# ------------------------------------------------------------

function main_benchmark()
    println("========================================")
    println("Takens + DMD comparison benchmark")
    println("========================================")

    q_values = 1:12
    r = 4

    # Structured synthetic control
    structured = synthetic_signal(250)

    # Null controls
    white = white_noise_signal(250, 3; seed=2)
    ar1 = ar1_signal(250, 3; α=0.9, seed=3)

    # Shuffled version of the structured signal
    shuffled = shuffled_snapshots(structured; seed=4)

    # PP data
    # Get Molecular Hamiltonian
    #data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P1-RHF_integrals.npz" #planar, 2 orbitals
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/C2H4_rotation/P1-RHF_integrals.npz" #twisted, 2 orbitals
    #data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/C2H4_rotation/0.0_tensors.npz" #14 orbitals,
    data = npzread(data_path)
    H0 = data["hc"][1]
    println("H0: ", H0)
    H1 = data["h1e"]
    H2 = data["h2e"]
    println("H0 shape: ", size(H0))
    println(typeof(H0))
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))
    Norbs = size(H1,1)  # number of spatial orbitals
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=false)
    
    # - - - Define Reference ket
    #ket, _  = QuantumChemQC.string_to_ket("1010") #Block
    ket, _  = QuantumChemQC.string_to_ket("1100") #Interleaved
    println("Initial state Expectation Value:")
    e1 = expectation_value(H,ket)
    @printf("E(0) = %.6f\n", e1)

    # - - - Define Excitation Operator
    O = QuantumChemQC.homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
    println("HOMO-LUMO excitation operator:")
    display(O)
    
    time_step = 0.5
    n_intervals = 250
    pp = run_pp_snapshots(ket, O, H, n_intervals, time_step; threshold=1e-10)
    shuffled_pp = shuffled = shuffled_snapshots(pp.w_snapshots; seed=4)

    # Benchmark all cases
    cases = BenchmarkResult[]
    push!(cases, benchmark_signal("synthetic", structured; dt=1.0, q_values=q_values, r=r))
    push!(cases, benchmark_signal("white noise", white; dt=1.0, q_values=q_values, r=r))
    #push!(cases, benchmark_signal("AR(1)", ar1; dt=1.0, q_values=q_values, r=r))
    #push!(cases, benchmark_signal("shuffled synthetic", shuffled; dt=1.0, q_values=q_values, r=r))
    push!(cases, benchmark_signal("PP weight snapshots", pp.w_snapshots; dt=time_step, q_values=1:min(12, length(pp.w_snapshots)-1), r=r))
    push!(cases, benchmark_signal("PP shuffled snapshots",shuffled_pp; dt=time_step, q_values=1:min(12, length(pp.w_snapshots)-1), r=r))

    println("\n--- Summary ---")
    for c in cases
        println("Case: ", c.name)
        println("  q values    : ", c.q_values)
        println("  fit error   : ", c.fit_error)
        println("  rank95      : ", c.rank95)
        println("  rank99      : ", c.rank99)
    end

    # Plots
    #p_fit = plot_benchmark_metric(cases; metric=:fit_error)
    #p_rank95, p_rank99 = plot_benchmark_ranks(cases)
    #p_sv = plot_singular_values_vs_q(cases)

    #display(p_fit)
    #display(p_rank95)
    #display(p_rank99)
    #display(p_sv)

    # PP-specific diagnostics
    display(plot_weight_heatmap(pp.w_snapshots; tgrid=pp.tgrid))

    return (; cases, pp)
end

function signals()
    println("========================================")
    println("PLOT SIGNALS")
    println("========================================")
    
    m = 250
    labels = ["x1" "x2" "x3"] # Horizontal matrix for series labeling

    # 1. Generate Signals (Vector of 3-element Vectors)
    structured_raw = synthetic_signal(m)
    white_raw      = white_noise_signal(m, 3; seed=2)
    ar1_raw        = ar1_signal(m, 3; α=0.9, seed=3)
    shuffled_raw   = shuffled_snapshots(structured_raw; seed=4)
    # PP data
    # Get Molecular Hamiltonian C2H4 rotation
    #data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P1-RHF_integrals.npz" #planar, 2 orbitals
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/C2H4_rotation/P1-RHF_integrals.npz" #twisted, 2 orbitals
    #data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/C2H4_rotation/0.0_tensors.npz" #14 orbitals
    data = npzread(data_path)
    H0 = data["hc"][1]
    println("H0: ", H0)
    H1 = data["h1e"]
    H2 = data["h2e"]
    println("H0 shape: ", size(H0))
    println(typeof(H0))
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))
    Norbs = size(H1,1)  # number of spatial orbitals
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=false)
    
    # - - - Define Reference ket
    #ket, _  = QuantumChemQC.string_to_ket("1010") #Block
    ket, _  = QuantumChemQC.string_to_ket("1100") #Interleaved
    println("Initial state Expectation Value:")
    e1 = expectation_value(H,ket)
    @printf("E(0) = %.6f\n", e1)

    # - - - Define Excitation Operator
    O = QuantumChemQC.homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
    println("HOMO-LUMO excitation operator:")
    display(O)

    pp = run_pp_snapshots(ket, O, H, 250, 0.5; threshold=1e-10)

    # 2. Conversion helper 
    # Transforms Vector{Vector{Float64}} into a Matrix{Float64} of size (m, 3)
    prepare(s) = reduce(hcat, s)'

    # 3. Plotting
    # We plot all 3 components on the same subplot for each signal type
    p1 = plot(prepare(structured_raw), title="Structured Signal", label=labels)
    p2 = plot(prepare(white_raw),      title="White Noise",       label=false)
    p3 = plot(prepare(ar1_raw),        title="AR(1) Process",     label=false)
    p4 = plot(prepare(shuffled_raw),   title="Shuffled Snapshots", label=false)
    p5 = plot(prepare(pp.w_snapshots), title="PP Weight Snapshots", label=false)

    # Combine into a 4-row stack
    combined_plot = plot(p1, p2, p3, p4, p5,
        layout=(5, 1), 
        size=(900, 1300), 
        link=:x,           # Links the x-axis for easier scrolling/comparison
        margin=5Plots.mm,
        ylabel="Value")

    display(combined_plot)
    return #combined_plot

end

# run 
results = main_benchmark()
signals()

#=
data_path = "/Users/admin/PycharmProjects/pyQCTools/DBF//N2_tensors_rhf_stable/1.1_tensors.npz" #Boys canonical localized
data = npzread(data_path)
Norbs = size(data["h1e"],1)  # number of spatial orbitals
H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=false)
QuantumChemQC.coeff_clip!(H, thresh=1e-6)
println("Hamiltonian terms:")
for (P, c) in H
        k = majorana_weight(P)
        a2 = abs2(c)
        println("P: ", P, " c: ", c, " |c|^2: ", a2, " Mweight: ", k)
end


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

println("O terms:")
for (P, c) in O
        k = majorana_weight(P)
        a2 = abs2(c)
        println("P: ", P, " |c|^2: ", a2, " Mweight: ", k)
end
=#