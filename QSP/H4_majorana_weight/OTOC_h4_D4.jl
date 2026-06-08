# ============================================================
# OTOC implementation using Pauli propagation
# ============================================================
using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Statistics
using Printf
using Random
using Plots
using NPZ

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

# ============================================================
# OTOC dynamics using the same PP loop as heisenberg_weight_dynamics
# ============================================================

function local_Z_op(Nqubits::Int, i::Int)
    return PauliSum(Pauli(Nqubits, Z = [i]))
end


function otoc_dynamics_pp(
    ket,
    H::PauliSum{N,T},
    cfg::PPConfig;
    j::Int = 1,
    probe_sites = collect(1:N),
    track_fidelity::Bool = true,
) where {N,T}

    Nqubits = N

    # W_j = exp(-i*pi*Z_j/2) = -i Z_j.
    # The global phase cancels in the OTOC, so propagate Z_j.
    Wt = -1im * local_Z_op(Nqubits, j)

    # Local OTOC probes V_i = Z_i
    Vs = Dict(i => local_Z_op(Nqubits, i) for i in probe_sites)

    generators, coeffs = QuantumChemQC.gens_from_H(H)

    tgrid = collect(0:cfg.dt:(cfg.n_steps * cfg.dt))

    F = zeros(ComplexF64, length(probe_sites), length(tgrid))
    fidelity = zeros(Float64, length(tgrid))
    raw_norms = zeros(Float64, length(tgrid))

    function record!(idx::Int, Wop)
        # Numerical norm check for W(t)
        _, norm2 = weight_profile(Wop; normalize = false)
        raw_norms[idx] = norm2

        if track_fidelity
            # Paper-style fidelity:
            # |<psi0|W_j(t)|psi0>|^2
            amp = expectation_value(Wop, ket)
            fidelity[idx] = abs2(amp)
        end

        for (a, i) in enumerate(probe_sites)
            V = Vs[i]

            # For alpha = pi, W is proportional to Z_j.
            # Global phase cancels, so this computes:
            # F_ij(t) = < W_j(t) Z_i W_j(t) Z_i >
            Ootoc = Wop * V * Wop * V
            F[a, idx] = expectation_value(Ootoc, ket)
        end
    end

    # t = 0
    record!(1, Wt)

    @printf("Total Pauli rotations: %d\n", length(coeffs) * cfg.n_steps)

    for step in 1:cfg.n_steps
        accumulated_error = 0.0 + 0.0im

        for a in eachindex(coeffs)
            θ = 2 * cfg.dt * coeffs[a]

            # Same PP update as your original code
            evolve!(Wt, generators[a], θ)

            QuantumChemQC.coeff_clip!(Wt; thresh = cfg.clip_rotation)
        end

        # Same step-level clipping as your original code
        QuantumChemQC.coeff_clip!(Wt; thresh = cfg.clip_step)

        record!(step + 1, Wt)
    end

    return (
        tgrid = tgrid,
        probe_sites = probe_sites,
        F = F,
        absF = abs.(F),
        realF = real.(F),
        imagF = imag.(F),
        fidelity = fidelity,
        raw_norms = raw_norms,
    )
end

#data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors/h4_RHF_H2H2_R_0p7414_tensors.npz"
data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/H4-PP/tensors/h4_RHF_H2H2_R_6p0_tensors.npz"

data = npzread(data_path)
H1 = data["h1e"]
Norbs = size(H1, 1)
Nqubits = 2 * Norbs

H = QuantumChemQC.molecular_hamiltonian(
    Norbs,
    data_path,
    NOI = false,
    block = false,
)

ket, _ = QuantumChemQC.string_to_ket("11110000")

cfg_otoc = PPConfig(
    50,       # n_steps
    0.1,      # dt
    1e-12,    # clip_rotation
    1e-12,    # clip_step
    false,
)

otoc = otoc_dynamics_pp(
    ket,
    H,
    cfg_otoc;
    j = 1,
    probe_sites = collect(1:Nqubits),
)

function plot_otoc_heatmap_pp(otoc; title = "OTOC |Fᵢⱼ(t)|")
    heatmap(
        otoc.tgrid,
        otoc.probe_sites,
        otoc.absF;
        xlabel = "time",
        ylabel = "probe qubit i",
        title = title,
        clims = (0, 1),
    )
end


function plot_otoc_fidelity_pp(otoc; title = "Fidelity |<ψ₀|Wⱼ(t)|ψ₀>|²")
    plot(
        otoc.tgrid,
        otoc.fidelity;
        xlabel = "time",
        ylabel = "fidelity",
        title = title,
        label = false,
        lw = 2,
        ylim = (0, 1.05),
    )
end


p1 = plot_otoc_heatmap_pp(otoc; title = "OTOC heatmap, R = 0.7414")
p2 = plot_otoc_fidelity_pp(otoc; title = "Fidelity, R = 0.7414")

display(p1)
display(p2)

function check_otoc_initial_value(otoc; atol = 1e-10)
    maxerr = maximum(abs.(otoc.F[:, 1] .- 1.0))
    @printf("max |F_i(0) - 1| = %.4e\n", maxerr)

    if maxerr > atol
        @warn "Initial OTOC is not close to 1. Check qubit indexing or operator definitions."
    end
end


function check_otoc_bounds(otoc; tol = 1e-8)
    max_abs = maximum(otoc.absF)
    @printf("max |F_i(t)| = %.8f\n", max_abs)

    if max_abs > 1 + tol
        @warn "OTOC magnitude exceeds 1. Check clipping, propagation convention, or norm drift."
    end
end


check_otoc_initial_value(otoc)
check_otoc_bounds(otoc)