using QuantumChemQC
using NPZ
using PauliOperators
#include("correlation_calculator.jl")

#Define datapath for the npz file containing molecular tensors
#data_path =  "/home/aaron/Vertical_excitation/furan_tensors_sto3g/Furan-STO3g_integrals.npz"
data_path =  "/home/aaron/Vertical_excitation/acetaldehyde_tensors_sto3g/acetaldehyde-STO3g_integrals.npz"
data = npzread(data_path)
Norbs = size(data["h1e"], 1)

H = QuantumChemQC.molecular_hamiltonian(
       Norbs,
       data_path;
       NOI = false,
       block = false,
    )

#display(H)
println("Total orbitals :", Norbs)
println("Total H terms :", length(H))

#hf_string = "1"^36 * "0"^22 #Furan
hf_string = "1"^24 * "0"^14 #Acetaldehyde

println("HF string length: ", length(hf_string))
ket, _ = QuantumChemQC.string_to_ket(hf_string)

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