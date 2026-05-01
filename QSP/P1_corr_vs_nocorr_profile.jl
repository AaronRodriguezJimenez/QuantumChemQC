using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots
using BenchmarkTools
using Profile
using ProfileView   # optional but recommended
using StatsPlots    # for plotting


function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P1-RHF_integrals.npz"
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
    
    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=true)

    return H
end

O = Pauli(4, X=[1,3])
O = PauliSum(O)
#O += Pauli(4, X=[2,3])
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("0000")
#display(H)

# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
n_intervals = 200
t = 10.0
dt = 0.1#t/n_intervals

evol_thresh = 0.0001 

# * * * PROFILE WITH AND WITHOUT CORRECTIONS * * *
println("========== WARM-UP ==========")

# Warm-up (compilation)
QuantumChemQC.QSP_evolution_op(ket, O, H, n_intervals, dt;
    thresh=evol_thresh, checkfile="P1_checkfile")

QuantumChemQC.QSP_evolution_op_no_corr(ket, O, H, n_intervals, dt;
    thresh=evol_thresh, checkfile="P1_checkfile")

println("Warm-up complete.\n")

# ============================================================
# ⏱️ BENCHMARKING (TIME + MEMORY)
# ============================================================

println("========== BENCHMARKING ==========")

bench_corr = @benchmark QuantumChemQC.QSP_evolution_op(
    $ket, $O, $H, $n_intervals, $dt;
    thresh=$evol_thresh,
    checkfile="P1_checkfile"
)

bench_nocorr = @benchmark QuantumChemQC.QSP_evolution_op_no_corr(
    $ket, $O, $H, $n_intervals, $dt;
    thresh=$evol_thresh,
    checkfile="P1_checkfile"
)

println("\n--- With corrections ---")
display(median(bench_corr))

println("\n--- Without corrections ---")
display(median(bench_nocorr))

# ============================================================
# 📊 SIMPLE PLOT (runtime comparison)
# ============================================================

times = [
    median(bench_corr).time / 1e6,     # ms
    median(bench_nocorr).time / 1e6
]

labels = ["With corrections", "No corrections"]

bar(labels, times, ylabel="Time (ms)", title="Runtime Comparison")

# ============================================================
# 🔍 PROFILING (WHERE TIME IS SPENT)
# ============================================================

println("\n========== PROFILING ==========")

# --- With corrections ---
Profile.clear()

@profile begin
    QuantumChemQC.QSP_evolution_op(
        ket, O, H, n_intervals, dt;
        thresh=evol_thresh,
        checkfile="P1_checkfile"
    )
end

println("\n--- Profile: With corrections ---")
Profile.print()

# Optional flame graph:
# ProfileView.view()


# --- Without corrections ---
Profile.clear()

@profile begin
    QuantumChemQC.QSP_evolution_op_no_corr(
        ket, O, H, n_intervals, dt;
        thresh=evol_thresh,
        checkfile="P1_checkfile"
    )
end

println("\n--- Profile: Without corrections ---")
Profile.print()

# Optional:
# ProfileView.view()