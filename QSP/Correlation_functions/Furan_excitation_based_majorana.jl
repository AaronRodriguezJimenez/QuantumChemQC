#!/usr/bin/env julia

"""
General two-point correlation-function calculator based on PauliOperators.jl
and SparsePauliVector propagation.

The calculator evaluates

    C(t) = <psi | A * B(t) | psi>

where B(t) is propagated in the Heisenberg picture under H.  By default
A = B, reproducing the autocorrelation convention used in the source code:

    C(t) = <psi | O(0) * O(t) | psi>.

Main public entry point:

    run_correlation_calculation(H, O, ket; ...)

Required packages:
    PauliOperators
    Plots

The Hamiltonian/operator may be built with PauliSum for convenience.  They
are converted to SparsePauliVector before the propagation hot loop.
"""

using PauliOperators
using Plots
using Printf
using Dates

# ============================================================================
# 1. CORRELATION-SPECIFIC TRUNCATION CORRECTION
# ============================================================================

"""
Track the accumulated change in

    <psi | A * B | psi>

due only to truncation events applied to the propagated operator B.

PauliOperators defines a correction as

    accumulated += observable_after - observable_before

so the corrected observable is

    observable_corrected = observable_raw - accumulated.
"""
mutable struct CorrelationCorrection{S,K} <: PauliOperators.CorrectionAccumulator
    left_state::S
    ket::K
    accumulated::ComplexF64
end

CorrelationCorrection(left_state, ket) =
    CorrelationCorrection(left_state, ket, 0.0 + 0.0im)

function PauliOperators._measure(O, corr::CorrelationCorrection)
    return correlation_value(corr.left_state, O, corr.ket)
end

function PauliOperators._accumulate!(corr::CorrelationCorrection, before, after)
    corr.accumulated += after - before
    return nothing
end

# ============================================================================
# 2. HELPERS
# ============================================================================

"""Convert an operator to SparsePauliVector, preserving an existing SPV."""
function to_sparse_pauli_vector(
    O;
    capacity_factor::Real = 2.0,
    coefficient_type = nothing,
)
    capacity_factor >= 1 || throw(ArgumentError("capacity_factor must be >= 1"))

    if O isa SparsePauliVector
        if coefficient_type === nothing || typeof(O).parameters[3] == coefficient_type
            return copy(O)
        else
            # Rebuild through PauliSum only when the user explicitly requests a
            # different engine coefficient type.
            return SparsePauliVector(
                PauliSum(O);
                T = coefficient_type,
                capacity_factor = Float64(capacity_factor),
            )
        end
    elseif coefficient_type === nothing
        return SparsePauliVector(O; capacity_factor = Float64(capacity_factor))
    else
        return SparsePauliVector(
            O;
            T = coefficient_type,
            capacity_factor = Float64(capacity_factor),
        )
    end
end

"""
Validate the Hermitian-coefficient assumption used by `trotterize`.

Pauli basis elements are Hermitian, so a Hermitian Hamiltonian should have
numerically real coefficients in this representation.
"""
function validate_hamiltonian_coefficients(
    H::SparsePauliVector;
    imag_tol::Real = 1e-10,
)
    imag_tol >= 0 || throw(ArgumentError("hamiltonian_imag_tol must be >= 0"))

    for (p, c) in H
        if abs(imag(c)) > imag_tol * max(1.0, abs(c))
            throw(ArgumentError(
                "Hamiltonian contains a significantly complex coefficient $(c) " *
                "for Pauli term $(p). trotterize uses real Hamiltonian " *
                "coefficients; check Hermiticity or increase hamiltonian_imag_tol.",
            ))
        end
    end

    return true
end

"""
Construct the deterministic pruning strategy used during propagation.

The three controls are independent:

- `evol_thresh > 0` enables coefficient truncation.
- `max_weight !== nothing` enables ordinary Pauli-weight truncation.
- `max_majorana_weight !== nothing` enables Majorana-weight truncation.

Consequently, this factory supports all eight enabled/disabled combinations.
When more than one filter is active, `CompositeTruncation` applies their
intersection.  For these deterministic pure-drop filters, a term survives only
if it satisfies every enabled cutoff.
"""
function make_truncation_strategy(
    evol_thresh::Real,
    max_weight::Union{Nothing,Integer},
    max_majorana_weight::Union{Nothing,Integer},
)
    evol_thresh >= 0 || throw(ArgumentError("evol_thresh must be >= 0"))
    max_weight === nothing || max_weight >= 0 ||
        throw(ArgumentError("max_weight must be >= 0 or nothing"))
    max_majorana_weight === nothing || max_majorana_weight >= 0 ||
        throw(ArgumentError("max_majorana_weight must be >= 0 or nothing"))

    if max_majorana_weight !== nothing &&
       !isdefined(PauliOperators, :MajoranaWeightTruncation)
        throw(ArgumentError(
            "The installed PauliOperators version does not provide " *
            "MajoranaWeightTruncation. Update PauliOperators or set " *
            "max_majorana_weight=nothing.",
        ))
    end

    # Build a tuple so CompositeTruncation retains concrete strategy types.
    # This construction happens once, outside the propagation hot loop.
    strategies = ()

    if evol_thresh > 0
        strategies = (strategies..., CoeffTruncation(Float64(evol_thresh)))
    end

    if max_weight !== nothing
        strategies = (strategies..., WeightTruncation(Int(max_weight)))
    end

    if max_majorana_weight !== nothing
        strategies = (
            strategies...,
            PauliOperators.MajoranaWeightTruncation(Int(max_majorana_weight)),
        )
    end

    if isempty(strategies)
        return NoTruncation()
    elseif length(strategies) == 1
        return first(strategies)
    else
        return CompositeTruncation(strategies...)
    end
end

"""Human-readable list of the enabled deterministic pruning components."""
function truncation_description(
    evol_thresh::Real,
    max_weight::Union{Nothing,Integer},
    max_majorana_weight::Union{Nothing,Integer},
)
    components = String[]
    evol_thresh > 0 && push!(components, "coefficient")
    max_weight !== nothing && push!(components, "Pauli weight")
    max_majorana_weight !== nothing && push!(components, "Majorana weight")
    return isempty(components) ? "none" : join(components, " + ")
end

"""
Build the fixed bra-side state |L> = A^dagger|psi> used in
C(t) = <L|B(t)|psi>.
"""
function build_left_state(left_operator::SparsePauliVector, ket)
    # Materialize the adjoint with the SPV scalar-multiplication API, then act
    # on the basis ket. This is done once per propagation run.
    left_dagger = left_operator' * 1
    return left_dagger * ket
end

"""
Evaluate <psi|A*B|psi> as <A^dagger psi | B psi>.

The left state is precomputed once.  The measurement then walks the SPV terms
directly, avoiding construction of either A*B or a temporary B|psi> KetSum.
"""
function correlation_value(left_state, evolved_operator::SparsePauliVector, ket)
    value = 0.0 + 0.0im

    for (p, c) in evolved_operator
        phase, target_ket = p * ket
        left_amplitude = get(left_state, target_ket, 0.0 + 0.0im)
        value += conj(left_amplitude) * phase * c
    end

    return ComplexF64(value)
end

"""Return a file-system-safe run name."""
function safe_run_name(name::AbstractString)
    cleaned = replace(strip(name), r"[^A-Za-z0-9_.-]+" => "_")
    return isempty(cleaned) ? "correlation" : cleaned
end

"""Best-effort package version string."""
function package_version_string(mod)
    try
        return string(Base.pkgversion(mod))
    catch
        return "unknown"
    end
end

# ============================================================================
# 3. SPARSE PAULI PROPAGATION
# ============================================================================

"""
Propagate B(t) interval-by-interval with SparsePauliVector.

`generators` and `angles` describe one time interval.  The same ordered
rotation sequence is repeated for every interval.

When `track_pruning_correction=true`, the returned corrected correlation is

    C_corrected(t) = C_raw(t) - sum(delta C_truncation)

with the sum accumulated from t=0 up to the requested time.
"""
function propagate_correlation(
    ket,
    initial_operator::SparsePauliVector,
    left_operator::SparsePauliVector,
    generators,
    angles,
    n_intervals::Integer;
    truncation = NoTruncation(),
    window::Integer = 1,
    track_pruning_correction::Bool = true,
)
    n_intervals >= 0 || throw(ArgumentError("n_intervals must be >= 0"))
    window >= 1 || throw(ArgumentError("window must be >= 1"))
    length(generators) == length(angles) ||
        throw(DimensionMismatch("generators and angles must have equal length"))

    Wt = copy(initial_operator)
    left_state = build_left_state(left_operator, ket)

    raw_Ct = Vector{ComplexF64}(undef, n_intervals + 1)
    corrected_Ct = Vector{ComplexF64}(undef, n_intervals + 1)
    accumulated_correction = Vector{ComplexF64}(undef, n_intervals + 1)
    term_counts = Vector{Int}(undef, n_intervals + 1)

    c0 = correlation_value(left_state, Wt, ket)
    raw_Ct[1] = c0
    corrected_Ct[1] = c0
    accumulated_correction[1] = 0.0 + 0.0im
    term_counts[1] = length(Wt)

    corr = CorrelationCorrection(left_state, ket)
    use_correction = track_pruning_correction && !(truncation isa NoTruncation)

    for k in 1:n_intervals
        if use_correction
            evolve!(
                Wt,
                generators,
                angles;
                window = Int(window),
                truncation = truncation,
                correction = corr,
            )
        else
            evolve!(
                Wt,
                generators,
                angles;
                window = Int(window),
                truncation = truncation,
            )
        end

        raw = correlation_value(left_state, Wt, ket)
        delta = use_correction ? corr.accumulated : (0.0 + 0.0im)

        idx = k + 1
        raw_Ct[idx] = raw
        accumulated_correction[idx] = delta
        corrected_Ct[idx] = raw - delta
        term_counts[idx] = length(Wt)
    end

    return (
        raw_Ct = raw_Ct,
        Ct = corrected_Ct,
        accumulated_correction = accumulated_correction,
        term_counts = term_counts,
        final_operator = Wt,
    )
end

# ============================================================================
# 4. OUTPUT WRITERS
# ============================================================================

function write_numerical_data(path, times, result)
    open(path, "w") do io
        println(io, "# PauliOperators.jl SparsePauliVector correlation data")
        println(io, "# C(t) = <psi | A B(t) | psi>")
        println(io, "# Corrected C(t) = raw C(t) - accumulated truncation change")
        println(io, "#")
        println(
            io,
            "# time  Re_C  Im_C  Abs_C  Phase_C  Re_C_raw  Im_C_raw  " *
            "Re_accumulated_correction  Im_accumulated_correction  n_terms",
        )

        for i in eachindex(times)
            c = result.Ct[i]
            raw = result.raw_Ct[i]
            corr = result.accumulated_correction[i]

            @printf(
                io,
                "%.4f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %d\n",
                times[i],
                real(c),
                imag(c),
                abs(c),
                angle(c),
                real(raw),
                imag(raw),
                real(corr),
                imag(corr),
                result.term_counts[i],
            )
        end
    end

    return path
end

function make_correlation_plot(
    path,
    times,
    result;
    title::AbstractString = "Correlation function",
    ylabel::AbstractString = "C(t)",
    show_raw::Bool = false,
)
    p = plot(
        times,
        real.(result.Ct);
        lw = 2,
        label = "Re C(t)",
        xlabel = "Time",
        ylabel = ylabel,
        title = title,
    )

    plot!(
        p,
        times,
        imag.(result.Ct);
        lw = 2,
        label = "Im C(t)",
    )

    if show_raw
        plot!(
            p,
            times,
            real.(result.raw_Ct);
            lw = 1.5,
            linestyle = :dash,
            label = "Re C_raw(t)",
        )
        plot!(
            p,
            times,
            imag.(result.raw_Ct);
            lw = 1.5,
            linestyle = :dash,
            label = "Im C_raw(t)",
        )
    end

    savefig(p, path)
    return path
end

function write_summary(
    path;
    run_name,
    H_sparse,
    O_sparse,
    left_sparse,
    ket,
    tmax,
    n_intervals,
    dt,
    trotter_order,
    n_trotter,
    generators,
    evol_thresh,
    max_weight,
    max_majorana_weight,
    truncation,
    window,
    track_pruning_correction,
    capacity_factor,
    coefficient_type,
    hamiltonian_imag_tol,
    stats,
    result,
    spv_ok,
    data_path,
    plot_path,
    metadata,
)
    allocated_mib = stats.bytes / 2.0^20

    open(path, "w") do io
        println(io, "============================================================")
        println(io, "PAULI CORRELATION CALCULATION SUMMARY")
        println(io, "============================================================")
        println(io, "Run name:                    ", run_name)
        println(io, "Timestamp:                   ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "Julia version:               ", VERSION)
        println(io, "PauliOperators version:      ", package_version_string(PauliOperators))
        println(io, "Plots version:               ", package_version_string(Plots))
        println(io)

        println(io, "CORRELATION DEFINITION")
        println(io, "  C(t) = <psi | A B(t) | psi>")
        println(io, "  B(t) is evolved in the Heisenberg picture.")
        println(io, "  The default call uses A = B = O.")
        println(io, "  Reference ket type:        ", typeof(ket))
        println(io)

        println(io, "SPARSE REPRESENTATION")
        println(io, "  Hamiltonian representation: ", typeof(H_sparse))
        println(io, "  Initial B representation:   ", typeof(O_sparse))
        println(io, "  Left A representation:      ", typeof(left_sparse))
        println(io, "  Hamiltonian terms:           ", length(H_sparse))
        println(io, "  Initial B terms:             ", length(O_sparse))
        println(io, "  Left A terms:                ", length(left_sparse))
        println(io, "  SPV capacity factor:         ", capacity_factor)
        println(io, "  Requested coefficient type: ", coefficient_type === nothing ? "preserve/input default" : string(coefficient_type))
        println(io, "  Hamiltonian imag tolerance:  ", hamiltonian_imag_tol)
        println(io, "  Final SPV integrity check:   ", spv_ok)
        println(io)

        println(io, "TIME EVOLUTION")
        @printf(io, "  Total time:                  %.16g\n", tmax)
        println(io, "  Number of output intervals:  ", n_intervals)
        @printf(io, "  dt per output interval:      %.16g\n", dt)
        println(io, "  Trotter order:               ", trotter_order)
        println(io, "  Trotter substeps/interval:   ", n_trotter)
        println(io, "  Rotations per interval:      ", length(generators))
        println(io, "  Total rotations:             ", length(generators) * n_intervals)
        println(io)

        println(io, "PAULI PROPAGATION / PRUNING")
        println(io, "  Active components:           ", truncation_description(evol_thresh, max_weight, max_majorana_weight))
        println(io, "  Strategy type:               ", typeof(truncation))
        println(io, "  Coefficient threshold:       ", evol_thresh > 0 ? @sprintf("%.6e", evol_thresh) : "disabled")
        println(io, "  Maximum Pauli weight:        ", max_weight === nothing ? "disabled" : string(max_weight))
        println(io, "  Maximum Majorana weight:     ", max_majorana_weight === nothing ? "disabled" : string(max_majorana_weight))
        println(io, "  Majorana convention:         Jordan-Wigner ordering of the Pauli string")
        println(io, "  Sparse evolution window:     ", window)
        println(io, "  Pruning correction tracked:  ", track_pruning_correction)
        println(io, "  Correction convention:       C_corrected = C_raw - sum(C_after - C_before)")
        println(io, "  Correction scope:            removes the direct measured change at pruning events;")
        println(io, "                               it is not a rigorous bound on later propagated error")
        println(io)

        println(io, "COMPUTATIONAL COST")
        @printf(io, "  Elapsed time (s):             %.6f\n", stats.time)
        @printf(io, "  Allocated memory (MiB):       %.3f\n", allocated_mib)
        @printf(io, "  GC time (s):                  %.6f\n", stats.gctime)
        println(io, "  Initial propagated terms:     ", result.term_counts[1])
        println(io, "  Final propagated terms:       ", result.term_counts[end])
        println(io, "  Max sampled terms:            ", maximum(result.term_counts))
        println(io)

        println(io, "FINAL CORRELATION")
        @printf(io, "  Re C(tmax):                   %.16e\n", real(result.Ct[end]))
        @printf(io, "  Im C(tmax):                   %.16e\n", imag(result.Ct[end]))
        @printf(io, "  |C(tmax)|:                    %.16e\n", abs(result.Ct[end]))
        @printf(io, "  Re raw C(tmax):               %.16e\n", real(result.raw_Ct[end]))
        @printf(io, "  Im raw C(tmax):               %.16e\n", imag(result.raw_Ct[end]))
        @printf(io, "  Re accumulated correction:    %.16e\n", real(result.accumulated_correction[end]))
        @printf(io, "  Im accumulated correction:    %.16e\n", imag(result.accumulated_correction[end]))
        println(io)

        if !isempty(metadata)
            println(io, "USER METADATA")
            for key in sort!(collect(keys(metadata)); by = string)
                println(io, "  ", key, ": ", metadata[key])
            end
            println(io)
        end

        println(io, "OUTPUT FILES")
        println(io, "  Numerical data:              ", abspath(data_path))
        println(io, "  Correlation plot:            ", abspath(plot_path))
        println(io, "  Summary:                     ", abspath(path))
        println(io, "============================================================")
    end

    return path
end

# ============================================================================
# 5. MAIN CALCULATOR API
# ============================================================================

"""
    run_correlation_calculation(H, O, ket; kwargs...)

Compute a two-point correlation function using SparsePauliVector propagation.

Required inputs
---------------
- `H`: Hamiltonian as a PauliSum or SparsePauliVector.
- `O`: operator B that is propagated in time.
- `ket`: reference state.

Important keyword inputs
------------------------
- `left_operator=O`: static operator A in C(t)=<psi|A B(t)|psi>.
- `tmax`: final time.
- `n_intervals`: number of output/evolution intervals.
- `evol_thresh=1e-3`: coefficient pruning threshold; set `0.0` to disable it.
- `max_weight=nothing`: maximum ordinary Pauli weight; `nothing` disables it.
- `max_majorana_weight=nothing`: maximum Majorana weight under the Jordan-Wigner
  encoding; `nothing` disables it. It may be combined independently with the
  coefficient and ordinary Pauli-weight filters.
- `trotter_order=1`: first- or second-order Trotterization.
- `n_trotter=1`: Trotter substeps inside each output interval.
- `window=1`: SPV merge/truncation cadence. Keep 1 for strict per-rotation pruning.
- `track_pruning_correction=true`: accumulate the directly measured change in the
  correlation caused by pruning and subtract it from the sampled raw value. This
  is a first-order bookkeeping correction, not a rigorous later-time error bound.
- `coefficient_type=nothing`: preserve input coefficient types. Set to `Float64`
  for Hermitian/real H, A, and B when valid; correlation measurement remains
  complex-safe because it is evaluated through state overlaps rather than A*B.
- `capacity_factor=2.0`: SPV preallocation headroom.
- `hamiltonian_imag_tol=1e-10`: tolerance used to verify that Hamiltonian
  coefficients are numerically real before Trotter decomposition.
- `output_dir="correlation_output"`: directory for outputs.
- `run_name="correlation"`: prefix for output filenames.
- `show_raw_on_plot=false`: overlay uncorrected correlation.
- `warmup=true`: compile a short run before timing.
- `metadata=Dict()`: arbitrary information written to the summary file.

Outputs
-------
Creates:
1. `<run_name>_correlation.txt` numerical data.
2. `<run_name>_correlation.pdf` plot.
3. `<run_name>_summary.txt` calculation information.

Returns a named tuple containing the arrays, final operator, timing statistics,
and output paths.
"""
function run_correlation_calculation(
    H,
    O,
    ket;
    left_operator = O,
    tmax::Real,
    n_intervals::Integer,
    evol_thresh::Real = 1e-3,
    max_weight::Union{Nothing,Integer} = nothing,
    max_majorana_weight::Union{Nothing,Integer} = nothing,
    trotter_order::Integer = 1,
    n_trotter::Integer = 1,
    window::Integer = 1,
    track_pruning_correction::Bool = true,
    coefficient_type = nothing,
    capacity_factor::Real = 2.0,
    hamiltonian_imag_tol::Real = 1e-10,
    output_dir::AbstractString = "correlation_output",
    run_name::AbstractString = "correlation",
    show_raw_on_plot::Bool = false,
    warmup::Bool = true,
    metadata = Dict{String,Any}(),
)
    tmax >= 0 || throw(ArgumentError("tmax must be >= 0"))
    n_intervals > 0 || throw(ArgumentError("n_intervals must be > 0"))
    trotter_order in (1, 2) || throw(ArgumentError("trotter_order must be 1 or 2"))
    n_trotter >= 1 || throw(ArgumentError("n_trotter must be >= 1"))
    window >= 1 || throw(ArgumentError("window must be >= 1"))

    dt = Float64(tmax) / Int(n_intervals)
    label = safe_run_name(run_name)
    mkpath(output_dir)

    # Convert both the Hamiltonian and operators to flat sorted storage.
    # For a PauliSum Hamiltonian this also makes the Trotter generator order
    # deterministic with respect to the SPV packed-key ordering.
    H_sparse = to_sparse_pauli_vector(
        H;
        capacity_factor = capacity_factor,
        coefficient_type = coefficient_type,
    )
    O_sparse = to_sparse_pauli_vector(
        O;
        capacity_factor = capacity_factor,
        coefficient_type = coefficient_type,
    )
    left_sparse = to_sparse_pauli_vector(
        left_operator;
        capacity_factor = capacity_factor,
        coefficient_type = coefficient_type,
    )

    validate_hamiltonian_coefficients(
        H_sparse;
        imag_tol = hamiltonian_imag_tol,
    )

    # One ordered rotation sequence for one output interval.
    generators, angles = trotterize(
        H_sparse,
        dt;
        n_trotter = Int(n_trotter),
        order = Int(trotter_order),
    )

    truncation = make_truncation_strategy(
        evol_thresh,
        max_weight,
        max_majorana_weight,
    )

    println()
    println("============================================================")
    println("SPARSE PAULI CORRELATION CALCULATION")
    println("============================================================")
    println("Run name:                 ", run_name)
    println("Hamiltonian terms:        ", length(H_sparse))
    println("Initial operator terms:   ", length(O_sparse))
    println("Left operator terms:      ", length(left_sparse))
    @printf("Total time:               %.8g\n", tmax)
    println("Intervals:                ", n_intervals)
    @printf("dt:                       %.8g\n", dt)
    println("Trotter order:            ", trotter_order)
    println("Trotter substeps:         ", n_trotter)
    println("Rotations per interval:   ", length(generators))
    println("Total rotations:          ", length(generators) * n_intervals)
    println("Pruning components:       ", truncation_description(evol_thresh, max_weight, max_majorana_weight))
    println("Truncation strategy type: ", typeof(truncation))
    println("Coefficient threshold:    ", evol_thresh > 0 ? @sprintf("%.6e", evol_thresh) : "disabled")
    println("Maximum Pauli weight:     ", max_weight === nothing ? "disabled" : max_weight)
    println("Maximum Majorana weight:  ", max_majorana_weight === nothing ? "disabled" : max_majorana_weight)
    println("SPV window:               ", window)
    println("Pruning correction:       ", track_pruning_correction)

    if warmup
        # Compile the actual SparsePauliVector evolution/correlation path without
        # contaminating the reported timing. One interval is sufficient.
        propagate_correlation(
            ket,
            O_sparse,
            left_sparse,
            generators,
            angles,
            1;
            truncation = truncation,
            window = window,
            track_pruning_correction = track_pruning_correction,
        )
    end

    GC.gc()
    stats = @timed propagate_correlation(
        ket,
        O_sparse,
        left_sparse,
        generators,
        angles,
        n_intervals;
        truncation = truncation,
        window = window,
        track_pruning_correction = track_pruning_correction,
    )

    result = stats.value
    times = collect(range(0.0; stop = Float64(tmax), length = n_intervals + 1))

    spv_ok = try
        PauliOperators.check_spv(result.final_operator)
    catch
        "check_spv unavailable"
    end

    data_path = joinpath(output_dir, "$(label)_correlation.txt")
    plot_path = joinpath(output_dir, "$(label)_correlation.pdf")
    summary_path = joinpath(output_dir, "$(label)_summary.txt")

    write_numerical_data(data_path, times, result)
    make_correlation_plot(
        plot_path,
        times,
        result;
        title = "Correlation function: $(run_name)",
        ylabel = "< A(0) B(t) >",
        show_raw = show_raw_on_plot,
    )

    write_summary(
        summary_path;
        run_name = run_name,
        H_sparse = H_sparse,
        O_sparse = O_sparse,
        left_sparse = left_sparse,
        ket = ket,
        tmax = Float64(tmax),
        n_intervals = Int(n_intervals),
        dt = dt,
        trotter_order = Int(trotter_order),
        n_trotter = Int(n_trotter),
        generators = generators,
        evol_thresh = Float64(evol_thresh),
        max_weight = max_weight,
        max_majorana_weight = max_majorana_weight,
        truncation = truncation,
        window = Int(window),
        track_pruning_correction = track_pruning_correction,
        capacity_factor = Float64(capacity_factor),
        coefficient_type = coefficient_type,
        hamiltonian_imag_tol = Float64(hamiltonian_imag_tol),
        stats = stats,
        result = result,
        spv_ok = spv_ok,
        data_path = data_path,
        plot_path = plot_path,
        metadata = metadata,
    )

    println()
    println("Calculation complete.")
    @printf("Elapsed time:             %.6f s\n", stats.time)
    @printf("Allocated memory:         %.3f MiB\n", stats.bytes / 2.0^20)
    println("Final propagated terms:   ", result.term_counts[end])
    @printf("Final Re C(t):            %.12e\n", real(result.Ct[end]))
    @printf("Final Im C(t):            %.12e\n", imag(result.Ct[end]))
    println("Numerical data:           ", abspath(data_path))
    println("Plot:                     ", abspath(plot_path))
    println("Summary:                  ", abspath(summary_path))
    println("============================================================")

    return (
        times = times,
        ReCt = real.(result.Ct),
        ImCt = imag.(result.Ct),
        Ct = result.Ct,
        raw_Ct = result.raw_Ct,
        accumulated_correction = result.accumulated_correction,
        term_counts = result.term_counts,
        final_operator = result.final_operator,
        stats = stats,
        truncation = truncation,
        paths = (
            data = data_path,
            plot = plot_path,
            summary = summary_path,
        ),
    )
end

# ============================================================================
# 6. EXAMPLE: H2 MOLECULE
# ============================================================================
#
# Keep this calculator file general. In a small driver script or in the REPL:
#
using QuantumChemQC
using NPZ
using PauliOperators
#include("correlation_calculator.jl")

data_path =  "/home/aaron/Vertical_excitation/furan_tensors_sto3g/Furan-STO3g_integrals.npz"
data = npzread(data_path)
Norbs = size(data["h1e"], 1)

H = QuantumChemQC.molecular_hamiltonian(
       Norbs,
       data_path;
       NOI = false,
       block = false,
    )

# Excitation operator
function excitation_operator_hermitian(N::Int, i::Int, j::Int)
    adag_i = jordan_wigner(i, N)
    a_j    = jordan_wigner(j, N)'

    O = adag_i * a_j

    return O + O'
end

O = excitation_operator_hermitian(58, 37, 35)

display(O)

hf_string = "1"^36 * "0"^22
println("HF string length: ", length(hf_string))
ket, _ = QuantumChemQC.string_to_ket(hf_string)

result = run_correlation_calculation(
       H,
       O,
       ket;
       tmax = 20.0,
       n_intervals = 100,
       evol_thresh = 0.0001,
       max_weight = nothing,
       max_majorana_weight = 4,
       trotter_order = 1,
       n_trotter = 1,
       window = 1,
       track_pruning_correction = true,
       coefficient_type = Float64,  # fast path when H/O coefficients are real
       output_dir = joinpath(@__DIR__, "correlation_output"),
       run_name = "furan_RHF",
       show_raw_on_plot = false,
       metadata = Dict(
           "Hamiltonian file" => data_path,
           "Reference state" => "1100",
           "Operator" => "T 18->19",
       ),
   )

# Truncation controls are independent. Examples:
#
#   evol_thresh = 0.0,    max_weight = nothing, max_majorana_weight = nothing # none
#   evol_thresh = 1e-4,   max_weight = nothing, max_majorana_weight = nothing # coefficient only
#   evol_thresh = 0.0,    max_weight = 4,       max_majorana_weight = nothing # Pauli weight only
#   evol_thresh = 0.0,    max_weight = nothing, max_majorana_weight = 4       # Majorana weight only
#   evol_thresh = 1e-4,   max_weight = 4,       max_majorana_weight = 4       # all three
#
# Any other enabled/disabled combination is accepted as well.
#
# For the exact autocorrelation convention of the supplied code, no
# `left_operator` keyword is needed because A defaults to O.
#
# For a general two-operator correlator <A B(t)>, pass:
#
#   left_operator = A
# ============================================================================
