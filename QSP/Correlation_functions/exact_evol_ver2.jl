using PauliOperators
using Printf
using QuantumChemQC
using NPZ
using Plots
using LinearAlgebra

function run_H(n)
    # n is the number of H atoms in the chain
    # Get Molecular Hamiltonian
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/small_molecules/HChains/h$n-RHF_tensors.npz"
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

n = 6 #4,6,8,10 
H = run_H(n)
#In the sto3g basis set each atom contributes with 1 spatial orbital, thus for n H atoms we have n spatial orbitals and 2n spin-orbitals = n qubits.
O = Pauli(2*n, X=[1,2,3,7,8,9])
O = PauliSum(O)
#O += Pauli(2*n, X=[2,3])
Hmat = Matrix(H)
Omat = Matrix(O)
ket = Ket{2*n}(0)#All-zero state
V0 = Vector(ket)

#1) When selecting the all-zero state this is an eigenstate of the Hamiltonian,
#   thus we need to apply a phase correction to the signal, multiplying by exp(-iE_0t) where E_0 
#   is the eigenvalue associated to the reference state (Ek in the code).
Ek = V0' * Hmat * V0
println("Reference state energy: ", Ek)

# One-time diagonalization
t1 = time()
E = eigen(Hmat)
V = E.vectors
Λ = E.values
println("Diagonalization time: ", time() - t1, " seconds")

# Transform state and operator to eigenbasis
psi0 = V0
psi0_eig = V' * psi0
Oeig = V' * Omat * V

function Ct_from_eigenbasis(t, psi0_eig, Oeig, Λ)
    phase_forward  = exp.(-im * Λ * t)
    phase_backward = exp.(+im * Λ * t)

    tmp = phase_forward .* psi0_eig
    tmp = Oeig * tmp
    tmp = phase_backward .* tmp
    tmp = Oeig * tmp
    
    return dot(conj(psi0_eig), tmp)
end

# Time grid
# Define time evolution parameters
# Circuit divided in k layers
# Thus total time (t) is divided in dt = t/k time intervals
k = 100
t = 10.00
dt = t/k
#dt = optimal_dt
nsamp = Int(k) + 1
tgrid = collect(range(0.0, stop=k*dt, length=nsamp))

# Step 2: Compute C(t)
Ct = ComplexF64[]
for time in tgrid
    push!(Ct, Ct_from_eigenbasis(time, psi0_eig, Oeig, Λ))
end

rRES = real(Ct)
iRES = imag(Ct)

#- - - - - PHASE CORRECTION (ONLY WHEN REF STATE IS ALSO AN EIGENSTATE) - - - - -
# SIGNAL PROCESSING multiply by exp(-iE_0t) to correct signal
signal = rRES .+ 1im * iRES;
phase = exp.(1im * Ek .* tgrid); #-1 is the eigenvalue associated with the eigenvector (|0>)
#corrected signal F(t) = exp(-iE_0t)*C(t)
F_exact = phase .* signal
rF = real(F_exact)
iF = imag(F_exact)

# Plot and print C(t) results
plt2 = plot(tgrid, real(F_exact), lw=2, label="Re(F(t)")
plt2 = plot!(tgrid, imag(F_exact), lw=2,label="Im(F(t)")

xlabel!(plt2, "Time"); ylabel!(plt2, "exp(-iE_0t) * < O(0)O(t) >")
title!(plt2, "H$n ,dt=$dt")
savefig(plt2, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_H$n-F_exact.pdf")


#Step 3: Plot Full C(t)
println("- - - C(t) - - - ")
#display(Ct)
plt = plot(tgrid, rRES, lw=2,# seriestype=:scatter,
          label="Re(C(t)")
plt = plot!(tgrid, iRES, lw=2,# seriestype=:scatter,
          label="Im(C(t)")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
title!(plt, "H$n ,dt=$dt")
savefig(plt, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_H$n-exact.pdf")

println("- - - Full signal C(t) - - - ")
@printf("dt     Re(C(t))    Im(C(t))   Re(F(t)  Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    #normop2 = rRES[i]^2 + iRES[i]^2
    #normop = sqrt(normop2)
    @printf("%.4f    %.6f   %.6f   %.6f  %.6f\n", interval, rRES[i], iRES[i], rF[i], iF[i])
end

println("- - - Corrected exact signal F(t) - - - ")
@printf("dt     Re(F(t))  Im(F(t))\n")
for (i,interval) in enumerate(tgrid)
    #normop2 = rRES[i]^2 + iRES[i]^2
    #normop = sqrt(normop2)
    @printf("%.4f    %.6f   %.6f\n", interval, rF[i], iF[i])
end

