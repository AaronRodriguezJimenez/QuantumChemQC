using PauliOperators
using Printf
using QuantumChemQC
using NPZ

"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/dipole_moment/c2h4-RHF_dip_mo.npz"
    data = npzread(data_path)
    d_mo = data["dip_op"]
    println("d_mo shape: ", size(d_mo))
    println(typeof(d_mo))

    N = size(d_mo)[2]  # number of spatial orbitals
    println("N:", N)
    @time D = QuantumChemQC.R_dipole_moment_op(N, data_path, block=false)
#    display(D)
    return D
end

function run_H()
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/dipole_moment/c2h4-RHF_tensors.npz"
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
n_intervals = 100
t = 2.5
dt = t/n_intervals
ket, _ = QuantumChemQC.string_to_ket("111111111111000000000000")
evol_thresh = 1e-1

println(typeof(D))
println(typeof(H))
rCtvals, iCtvals, tgrid = QuantumChemQC.QSP_evolution_op(H, D, n_intervals, 
                                                        dt, ket;
                                                        thresh=evol_thresh)

# Number of snapshots actually returned
nsnap = length(rRES)
println("* * * * Number of snapshots collected: $nsnap")

# Print C(t) results
plt = plot(tgrid, rRES, lw=2, seriestype=:scatter,
          label="Re(C(t), th=$threshold")
plt = plot!(tgrid, iRES, lw=2, seriestype=:scatter,
          label="Im(C(t), th=$threshold")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
title!(plt, "N=$N, J=$Jx, Jz=$Jz,dt=$dt")
savefig(plt, "QSP_C2H4.pdf")

println("- - - Sanity Check: |C(t)|^2 - - - ")
@printf("idx    dt     Re(C(t))    Im(C(t))   |C(t)|^2   |C(t)|\n")
for (i,interval) in enumerate(tgrid)
    normop2 = rRES[i]^2 + iRES[i]^2
    normop = sqrt(normop2)
    @printf("%.2s    %.4f    %.6f   %.6f   %.6f  %.6f\n", i, interval, rRES[i], iRES[i], normop2, normop)
end

println("Initial operator O:")
display(o)
println("Initial eigenstate |0> ")
display(ket)
display(Vector(ket))
println("- - - H|0> - - - -")
Hpsi0 = Hmat*Vector(ket) 
display(Hpsi0) 
println("- - - <0|H|0> - - - -")
function compute_ref_expval(H, k)
    ref_expval = 0
    for (p, c) in H
        ref_expval += c * expectation_value(p, ket)
    end
    return ref_expval
end

ref_expval = compute_ref_expval(H, k)
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
          label="Re(F(t), th=$threshold")
plt2 = plot!(tgrid, imag(F), lw=2, #seriestype=:scatter,
          label="Im(F(t), th=$threshold")

xlabel!(plt2, "Time"); ylabel!(plt2, "exp(-iE_t) * < O(0)O(t) >")
title!(plt2, "N=$N, J=$Jx, Jz=$Jz,dt=$dt")
savefig(plt2, "QSP_C2H4_corrct.pdf")

println("- - - Compare signals  - - - ")
@printf("dt   Re(C(t))    Im(C(t))   Re(F(t))    Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    normop2 = rRES[i]^2 + iRES[i]^2
    normop = sqrt(normop2)
    @printf("%.4f     %.6f    %.6f    %.6f   %.6f\n", interval, rRES[i], iRES[i], real(F[i]), imag(F[i]))
end

plt3 = plot(tgrid, rRES, lw=2, #seriestype=:scatter,
          label="Re(C(t), th=$threshold")
plt3 = plot!(tgrid, iRES, lw=2, #seriestype=:scatter,
          label="Im(C(t), th=$threshold")
plt3 = plot!(tgrid, real(F), lw=2, #seriestype=:scatter,
          label="Re(F(t), th=$threshold")
plt3 = plot!(tgrid, imag(F), lw=2, #seriestype=:scatter,
          label="Im(F(t), th=$threshold")

xlabel!(plt3, "Time"); ylabel!(plt3, "exp(-iE_kt) * < O(0)O(t) >")
title!(plt3, "Signal Comparison N=$N, J=$Jx, Jz=$Jz,dt=$dt")
savefig(plt3, "QSP_C2H4_comparison.pdf")