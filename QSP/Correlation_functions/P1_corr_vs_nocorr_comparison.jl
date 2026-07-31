using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots

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
rRESc, iRESc, tgridc = QuantumChemQC.QSP_evolution_op(ket, O, H, n_intervals, dt,
                                                        thresh=evol_thresh,
                                                        checkfile="P1_checkfile")

rRES, iRES, tgrid = QuantumChemQC.QSP_evolution_op_no_corr(ket, O, H, n_intervals, dt,
                                                           thresh=evol_thresh, 
                                                           checkfile="P1_checkfile")

# Number of snapshots actually returned
nsnap = length(rRES)
println("* * * * Number of snapshots collected: $nsnap")

# Print C(t) results
plt = plot(tgrid, rRESc, lw=2, seriestype=:scatter,color=:red,
          label="corr Re(C(t), th=$evol_thresh")
#plt = plot!(tgrid, iRESc, lw=2, seriestype=:scatter, color=:red,
#          label="corrIm(C(t), th=$evol_thresh")

plt = plot!(tgrid, rRES, lw=2, #seriestype=:scatter,
          label="Im(C(t), th=$evol_thresh", color=:black)
#plt = plot!(tgrid, iRES, lw=2, #seriestype=:scatter,
#          label="Im(C(t), th=$evol_thresh", color=:green)          

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
savefig(plt, "/Users/admin/VSCProjects/QuantumChemQC/QSP/P1_PP_comparison_$evol_thresh.pdf")

@printf("dt     Re(C(t))    Im(C(t))\n")
for (i,interval) in enumerate(tgrid)
    @printf("%.4f    %.6f    %.6f\n", interval, rRESc[i], iRESc[i])
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
signalc = rRESc .+ 1im * iRESc;
signal = rRES .+ 1im * iRES;


#corrected signal F(t) = exp(-iE_0t)*C(t)
phasec = exp.(1im * Ek .* tgridc);
Fc = phasec .* signalc
phase = exp.(1im * Ek .* tgrid);
F = phase .* signal

# Print C(t) results
plt2 = plot(tgrid, real(Fc), lw=2, seriestype=:scatter, color=:red,
          label="corr Re(F(t), th=$evol_thresh")
#plt2 = plot!(tgrid, imag(Fc), lw=2, seriestype=:scatter,
#          label="corr Im(F(t), th=$evol_thresh")

plt2 = plot!(tgrid, real(F), lw=2, #seriestype=:scatter,
          label="Re(F(t), th=$evol_thresh", color=:black)
#plt2 = plot!(tgrid, imag(F), lw=2, #seriestype=:scatter,
 #         label="Im(F(t), th=$evol_thresh")
          
xlabel!(plt2, "Time"); ylabel!(plt2, "exp(-iE_t) * < O(0)O(t) >")
savefig(plt2, "/Users/admin/VSCProjects/QuantumChemQC/QSP/P1_PP_F_comparison_$evol_thresh.pdf")

println("- - - Corrected signal  - - - ")
@printf("dt   Re(F(t))    Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    @printf("%.4f   %.6f   %.6f\n", interval, real(F[i]), imag(F[i]))
end
