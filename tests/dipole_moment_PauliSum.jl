using PauliOperators
using Printf
using QuantumChemQC
using NPZ

"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/tests/dipole_moment/c2h4-RHF_dip_mo.npz"
    data = npzread(data_path)
    d_mo = data["dip_op"]
    println("d_mo shape: ", size(d_mo))
    println(typeof(d_mo))

    N = size(d_mo)[2]  # number of spatial orbitals
    println("N:", N)
    @time D = QuantumChemQC.R_dipole_moment_op(N, data_path, block=false)
    display(D)

    return D
end

D = run_dip()