using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots

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
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=true)

    return H
end

O = Pauli(4, X=[1,3])
O = PauliSum(O)
O += Pauli(4, X=[2,3])
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("0000")
#display(H)

# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
n_intervals = 200
t = 50.0
dt = t/n_intervals

evol_thresh = 0.00001 
rRES, iRES, tgrid = QuantumChemQC.QSP_evolution_op(ket, O, H, n_intervals, dt,
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
savefig(plt, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_H2_deterministic_$evol_thresh.pdf")

println("- - - Sanity Check: |C(t)|^2 - - - ")
@printf("idx    dt     Re(C(t))    Im(C(t))\n")
for (i,interval) in enumerate(tgrid)
    normop2 = rRES[i]^2 + iRES[i]^2
    normop = sqrt(normop2)
    @printf("%.4f    %.6f    %.6f\n", interval, rRES[i], iRES[i])
end

println("Initial eigenstate |0> ")
display(ket)
println("- - - <0|H|0> - - - -")
function compute_ref_expval(H, ket)
    ref_expval = 0
    for (p, c) in H
        ref_expval += c * expectation_value(p, ket)
    end
    return ref_expval
end

ref_expval = compute_ref_expval(H, ket)
println("Reference energy: ", ref_expval)

#- - - - - PHASE CORRECTION (ONLY WHEN REF STATE IS ALSO AN EIGENSTATE) - - - - -
# SIGNAL PROCESSING multiply by exp(-iE_0t) to correct signal
Ek = ref_expval
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
savefig(plt2, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_H2_deterministic_F_$evol_thresh.pdf")

println("- - - Corrected signal  - - - ")
@printf("dt   Re(F(t))    Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    @printf("%.4f   %.6f   %.6f\n", interval, real(F[i]), imag(F[i]))
end
