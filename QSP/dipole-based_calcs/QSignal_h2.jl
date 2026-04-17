using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots

"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/dipole_moment/h2-RHF_dip_mo.npz"
    data = npzread(data_path)
    d_mo = data["dip_op"]
    println("d_mo shape: ", size(d_mo))
    println(typeof(d_mo))

    N = size(d_mo)[2]  # number of spatial orbitals
    println("N:", N)
    @time D = QuantumChemQC.R_dipole_moment_op(N, data_path, block=false)
    coeff_clip!(D, thresh=1.0e-6)
    return D
end

function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/dipole_moment/h2-RHF_tensors.npz"
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
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=false)

    return H
end

D = run_dip()
H = run_H()
display(H)
# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
n_intervals = 200
t = 50.0 #Total Time Evolution
dt = t/n_intervals
ket, _ = QuantumChemQC.string_to_ket("1100")
evol_thresh = 0.00001

println(typeof(D))
println(typeof(H))
rRES, iRES, tgrid = QuantumChemQC.QSP_evolution_op(ket, D, H, n_intervals, dt,
                                                        thresh=evol_thresh)

# Number of snapshots actually returned
nsnap = length(rRES)
println("* * * * Number of snapshots collected: $nsnap")

# Print C(t) results
plt = plot(tgrid, rRES, lw=2, #seriestype=:scatter,
          label="Re(C(t), th=$evol_thresh")
plt = plot!(tgrid, iRES, lw=2, #seriestype=:scatter,
          label="Im(C(t), th=$evol_thresh")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
savefig(plt, "QSP_H2.pdf")

println("- - - Sanity Check: |C(t)|^2 - - - ")
@printf("idx    dt     Re(C(t))    Im(C(t))   |C(t)|^2   |C(t)|\n")
for (i,interval) in enumerate(tgrid)
    normop2 = rRES[i]^2 + iRES[i]^2
    normop = sqrt(normop2)
    @printf("%.2s    %.4f    %.6f   %.6f   %.6f  %.6f\n", i, interval, rRES[i], iRES[i], normop2, normop)
end

println("Dipole operator:")
display(D)
println("Initial eigenstate |0> ")
display(ket)
println("- - - <0|H|0> - - - -")
function compute_ref_expval(H, k)
    ref_expval = 0
    for (p, c) in H
        ref_expval += c * expectation_value(p, ket)
    end
    return ref_expval
end

ref_expval = compute_ref_expval(H, ket)
println(ref_expval)

#- - - - - PHASE CORRECTION (ONLY WHEN REF STATE IS ALSO AN EIGENSTATE) - - - - -
# SIGNAL PROCESSING multiply by exp(-iE_0t) to correct signal
Ek = 0.0#-8 #Eigenvalue associated to ref state in case of correction.
signal = rRES .+ 1im * iRES;
phase = exp.(1im * Ek .* tgrid); #-1 is the eigenvalue associated with the eigenvector (|0>)

#corrected signal F(t) = exp(-iE_0t)*C(t)
F = phase .* signal

# Print C(t) results
plt2 = plot(tgrid, real(F), lw=2, #seriestype=:scatter,
          label="Re(F(t), th=$evol_thresh")
plt2 = plot!(tgrid, imag(F), lw=2, #seriestype=:scatter,
          label="Im(F(t), th=$evol_thresh")

xlabel!(plt2, "Time"); ylabel!(plt2, "exp(-iE_t) * < O(0)O(t) >")
savefig(plt2, "QSP_H2_corrct.pdf")

println("- - - Compare signals  - - - ")
@printf("dt   Re(C(t))    Im(C(t))   Re(F(t))    Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    @printf("%.4f     %.6f    %.6f    %.6f   %.6f\n", interval, rRES[i], iRES[i], real(F[i]), imag(F[i]))
end

plt3 = plot(tgrid, rRES, lw=2, #seriestype=:scatter,
          label="Re(C(t), th=$evol_thresh")
plt3 = plot!(tgrid, iRES, lw=2, #seriestype=:scatter,
          label="Im(C(t), th=$evol_thresh")
plt3 = plot!(tgrid, real(F), lw=2, #seriestype=:scatter,
          label="Re(F(t), th=$evol_thresh")
plt3 = plot!(tgrid, imag(F), lw=2, #seriestype=:scatter,
          label="Im(F(t), th=$evol_thresh")

xlabel!(plt3, "Time"); ylabel!(plt3, "exp(-iE_kt) * < O(0)O(t) >")
title!(plt3, "Signal Comparison dt=$dt")
savefig(plt3, "QSP_H2_comparison.pdf")

# Lowest FCI singlet states:
# [-1.13727017 -0.53247901 -0.16990139  0.47983612]
# S0-S1 = -0.60479116
# S0-S2 = -0.96736878
# S0-S3 = -1.61710629

## Compare signal with "theoretical one"
# We assume that the dipole moment will connect ground and excited states so
# according to thi hypothesis, the correct signal must be
#f(x) = (1/2)*cos(-0.60479116x) + (1/4)*cos(-0.96736878x) + (1/4)*cos(-1.61710629x)
f(x) = 0.38*cos(0.6048x) + 0.05*cos(1.0123x) + 0.10*cos(0.9673x) -0.04*cos(0.6497x)
#
# Evaluate on the grid
y = f.(tgrid)

plt4 = plot(tgrid, rRES, lw=2, seriestype=:scatter,
          label="Re(C(t)")
plt4 = plot!(tgrid, y, lw=2, xlabel="t", ylabel="C(t)",  
             title="Signal vs theoretical f(x)", label="f(t)",
             legend=true)

savefig(plt4, "QSP_H2_signal_vs_theory.pdf")

# Write an output text file
open("output_H2_signal.txt", "w") do io
    println(io, "dt   Re(C(t))    Im(C(t))   Re(F(t))    Im(F(t))")
    for (i, interval) in enumerate(tgrid)
        line = @sprintf("%8.4f  %10.6f  %10.6f  %10.6f  %10.6f",
                        interval, rRES[i], iRES[i], real(F[i]), imag(F[i]))
        println(io, line)
    end
end