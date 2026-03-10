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
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P3-RHF_dip_mo.npz"
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
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/Ethylene_and_polyenes/P3-RHF_integrals.npz"
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
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=false, block=false)

    return H
end

D = run_dip()
H = run_H()
Hmat = Matrix(H)
Omat = Matrix(D)

n_intervals = 100
t = 50.0 #Total Time Evolution
dt = t/n_intervals
ket, _ = QuantumChemQC.string_to_ket("111111000000")
V0 = Vector(ket)

# Compute O matrix elements:
# O_{a0} = <v_a | O | v_0>
#Oa0 = V' * (Omat * V0)
#println("O_{a0} = <v_a | O | v_0>")
#display(Oa0)

# Time evolution unitaries:
function Ut(Hmat, dt)
    return exp(-im * dt * Hmat )
end

# Time grid
# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
k = 100
t = 50.00
dt = t/k
nsamp = Int(k) + 1
tgrid = collect(range(0.0, stop=k*dt, length=nsamp))

# Step 2: Compute C(t)
W = deepcopy(Omat)
Ct = Vector{ComplexF64}([])
for time in tgrid
    U = Ut(Hmat, time)
    Udag = U'
    UdagWU = Udag * W * U
    WWt = W * UdagWU
    res = V0' * WWt * V0
   # res = res * (1/norm(res))
    push!(Ct, res)
end

#Step 3: Plot C(t)
println("- - - C(t) - - - ")
#display(Ct)

rRES = real(Ct)
iRES = imag(Ct)
plt = plot(tgrid, rRES, lw=2,# seriestype=:scatter,
          label="Re(C(t)")
plt = plot!(tgrid, iRES, lw=2,# seriestype=:scatter,
          label="Im(C(t)")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
title!(plt, "P3(6e,6o) ,dt=$dt")
savefig(plt, "QSP_P3_exact.pdf")

println("- - - Sanity Check: saving signal - - - ")
open("/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_P3_exact_signals.txt", "w") do io
    @printf(io,"dt     Re(C(t))    Im(C(t))   |C(t)|^2   |C(t)|\n")
    for (i, interval) in enumerate(tgrid)
        normop2 = rRES[i]^2 + iRES[i]^2
        normop = sqrt(normop2)
        @printf(io, "%.4f     %.6f    %.6f    %.6f   %.6f\n",
                interval, rRES[i], iRES[i], normop2, normop)
    end
end