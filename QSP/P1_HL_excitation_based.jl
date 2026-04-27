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
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=false)

    return H
end

# - - - Get Molecular Hamiltonian
H = run_H()
# - - - Define Reference ket
ket, _ = QuantumChemQC.string_to_ket("1100")
# - - - Define Excitation Operator
O = QuantumChemQC.homo_lumo_excitation_op(ket, spin_conserving=true, block=false)
display(O)
display(O*ket)
# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
n_intervals = 400
t = 100.0
dt = t/n_intervals

evol_thresh = 0.00001 
rRESc, iRESc, tgridc = QuantumChemQC.QSP_evolution_op(ket, O, H, n_intervals, dt,
                                                        thresh=evol_thresh)

rRES, iRES, tgrid = QuantumChemQC.QSP_evolution_op_no_corr(ket, O, H, n_intervals, dt,
                                                           thresh=evol_thresh)

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
savefig(plt, "/Users/admin/P1_PP_$evol_thresh.pdf")

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

println("- - - Signal  - - - ")
@printf("dt   Re(C(t))    Im(C(t))\n")
for (i,interval) in enumerate(tgrid)
    @printf("%.4f   %.6f   %.6f\n", interval, rRES[i], iRES[i])
end
