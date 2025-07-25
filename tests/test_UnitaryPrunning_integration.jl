using PauliOperators
using UnitaryPruning
using QuantumChemQC
using Tullio
using Printf
using IterTools
using LinearAlgebra
using Printf

"""
  Function returnign the Hamiltonian for H2 at a given interatomic distance
    d -> The interatomic distance
"""
function H2_hamiltonian_diss(d)
    mol = """
    0 1
    H  0.0000   0.000  0.000
    H  0.0000   0.000  $d
    """

    bset, p = QuantumChemQC.molecule(mol, "sto-3g", spherical = false)
    scf_obj = QuantumChemQC.SCF(mol, bset, p)

    # Atomic integrals
    ao_hcore = scf_obj.T + scf_obj.V
    ao_eris = scf_obj.I

    # Transform Atomic integrals into Molecular integrals
    C = scf_obj.C # MO coeffs
    mo_hcore , mo_eris = QuantumChemQC.ao2mo_coefficients(C, ao_hcore, ao_eris)

    # Spinorbital tensors
    core =0 
    h0 = scf_obj.Enuc + core 
    h1, h2 = QuantumChemQC.get_spin_orbital_tensors(mo_hcore, mo_eris)
    N = size(h1, 1)
    qubit_paulis, qubit_coeffs = QuantumChemQC.qubit_hamiltonian(N, h0, h1, h2) 

    gen, coeff = QuantumChemQC.combine_terms(qubit_paulis, qubit_coeffs)
return gen, coeff;
end


#- - - Test UnitaryPrunning 
function run(; d = 0.74, w_type = "Majorana", max_weight=1)
    N = 4
    ket = Ket(N, 3) #HF STATE
    display(ket)
    o = Pauli(N, Z=[1,2])

    generators, parameters = H2_hamiltonian_diss(d)

    for (pauli, coeff) in zip(generators, parameters)
        println(@sprintf("%+.6f %+.6fim | %s", real(coeff), imag(coeff), string(pauli)))
        #display(pauli)
    end

    #Call to bfs bfs_evolution_test based on weight
    ei, nops = UnitaryPruning.bfs_evolution_weight(generators, parameters, PauliSum(o), ket, w_type, max_weight=max_weight)
    
    # Exact evolution (Heisenberg picture)
    U = UnitaryPruning.build_time_evolution_matrix(generators, parameters)
    o_mat = Matrix(o)
    m = diag(U'*o_mat*U)
    expval = m[1]

    abs_err = abs(real(expval) - real(ei) )
    totops = sum(nops)
    println("Exact :", real(expval), " Approx :", real(ei))
    println("Total Operators: $totops,  Absolute Error: $abs_err")
    return abs_err;

end 

run(d = 4.00,  w_type = "Pauli", max_weight=1)

