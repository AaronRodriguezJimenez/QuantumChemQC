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
    println(H0)
    println("H1 shape: ", size(H1))
    println(typeof(H1))
    println("H2 shape: ", size(H2))
    println(typeof(H2))

    Norbs = size(H1,1)  # number of spatial orbitals
    
    #H  = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)
    @time H = QuantumChemQC.molecular_hamiltonian(Norbs, data_path, NOI=true, block=true)

    return H
end

# Evolved operator O:
O = Pauli(4, X=[1,3])
O = PauliSum(O)
O += Pauli(4, X=[2,3])

H = run_H()
Hmat = Matrix(H)
Omat = Matrix(O)
ket, _ = QuantumChemQC.string_to_ket("0000") #highest spin state in Block encoding (the triplet state)
V0 = Vector(ket)

#1) When selecting the highest possible spin state, this is an eigenstate of the Hamiltonian,
#   thus we need to apply a phase correction to the signal, multiplying by exp(-iE_0t) where E_0 
#   is the eigenvalue associated to the reference state (Ek in the code).
Ek = V0' * Hmat * V0
println("Reference state energy: ", Ek)

# To make sure that ODMD can be used successfully, we need an optimal time step (dt) that allows 
# to capture the dynamics of the system without losing information. This is ~2pi/|H|max, where |H|max is the largest eigenvalue of the Hamiltonian.
# Based on Klymko's paper, if we know that E_min <= <H> <= E_max, then we can set dt = 2pi/(E_max - E_min) to capture the dynamics effectively.
# In our case, we can Estimate an "optimal" based on the energy of Hf and that of the high spin state El
Hf_ket, _ = QuantumChemQC.string_to_ket("1010") #Singlet S0 in Block encoding
eHF = Vector(Hf_ket)' * Hmat * Vector(Hf_ket)
println("HF energy: ", eHF)
optimal_dt = 2*pi/abs(Ek - eHF)
println("Optimal dt: 2π/|Ek - eHF| = ", optimal_dt)

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
dt = 0.25#t/k
#dt = optimal_dt
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

rRES = real(Ct)
iRES = imag(Ct)

#- - - - - PHASE CORRECTION (ONLY WHEN REF STATE IS ALSO AN EIGENSTATE) - - - - -
# SIGNAL PROCESSING multiply by exp(-iE_0t) to correct signal
signal = rRES .+ 1im * iRES;
phase = exp.(1im * Ek .* tgrid); #-1 is the eigenvalue associated with the eigenvector (|0>)
#corrected signal F(t) = exp(-iE_0t)*C(t)
F = phase .* signal
rF = real(F)
iF = imag(F)

# Print C(t) results
plt2 = plot(tgrid, real(F), lw=2, label="Re(F(t)")
plt2 = plot!(tgrid, imag(F), lw=2,label="Im(F(t)")

xlabel!(plt2, "Time"); ylabel!(plt2, "exp(-iE_0t) * < O(0)O(t) >")
title!(plt2, "Corrected signal F(t) for C2H4(2e, 2o) ,dt=$dt")
savefig(plt2, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_P1_F_exact.pdf")


#Step 3: Plot Full C(t)
println("- - - C(t) - - - ")
#display(Ct)
plt = plot(tgrid, rRES, lw=2,# seriestype=:scatter,
          label="Re(C(t)")
plt = plot!(tgrid, iRES, lw=2,# seriestype=:scatter,
          label="Im(C(t)")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
title!(plt, "P1(2e, 2o) ,dt=$dt")
savefig(plt, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_P1_exact.pdf")

println("- - - Sanity Check: |C(t)|^2 - - - ")
@printf("dt     Re(C(t))    Im(C(t))\n")
for (i,interval) in enumerate(tgrid)
    #normop2 = rRES[i]^2 + iRES[i]^2
    #normop = sqrt(normop2)
    @printf("%.4f    %.6f   %.6f\n", interval, rRES[i], iRES[i])
end

println("- - - Sanity Check: |C(t)|^2 - - - ")
@printf("dt     Re(F(t))  Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    #normop2 = rRES[i]^2 + iRES[i]^2
    #normop = sqrt(normop2)
    @printf("%.4f    %.6f   %.6f\n", interval, rF[i], iF[i])
end