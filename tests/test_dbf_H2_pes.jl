using PyCall
using PauliOperators
using DBF
using Printf
using QuantumChemQC

"""
Compute PES of H2 molecule using DBF for ground state energy.
This script sweeps H-H bond distances (in Angstrom) and reports HF and DBF energies.
"""

# Helper bit utilities
get_bit(k::Union{Ket{N}, Bra{N}}, i::Integer) where N = Int((k.v >> (i-1)) & Int128(1))

function get_on_bits(v::Int128)
    out = Int[]
    w = v
    pos = 1
    while w != 0
        if (w & Int128(1)) != 0
            push!(out, pos)
        end
        w >>= 1
        pos += 1
    end
    return out
end

function occvec(k::Union{Ket{N}, Bra{N}}) where N
    [ get_bit(k, i) for i in 1:N ]
end

# Build Pauli Hamiltonian for H2 at bond length R (Angstrom)
function build_hamiltonian(R::Float64)
    # ensure pyqctools reload and imports are fresh
    importlib = pyimport("importlib")
    importlib.reload(pyimport("pyqctools"))
    pq = pyimport("pyqctools")
    Hfcns = pyimport("pyqctools.ham_fcns")

    # Build geometry string for PySCF (Angstrom)
    geom = "H 0 0 0; H 0 0 $(R)"
    println("Building molecule at R = $(R) Å; geometry: ", geom)

    mol = pyimport("pyscf.gto").Mole()
    mol.atom = geom
    mol.basis = "3-21g" # "sto-3g"
    mol.unit = "Angstrom"   # explicitly set units
    mol.build()

    # Run RHF
    mf = pyimport("pyscf.scf").RHF(mol).run()
    E_HF = mf.e_tot
    println("  Hartree-Fock Energy: ", E_HF)

    # get spin-orbital tensors
    H0, H1, H2 = Hfcns.get_tensors(mol, mf, localized=false)
    n = size(H1,1)  # number of spin-orbitals (qubits)
    N_total = n

    # Build Pauli Hamiltonian from QuantumChemQC helper
    H = QuantumChemQC.PauliSum_hamiltonian(n, H0, H1, H2)

    return H, N_total, E_HF
end

# occupation vector of length N (Int elements 0/1)
function occvec(k::Union{Ket{N}, Bra{N}}) where N
    [ get_bit(k, i) for i in 1:N ]
end

# DBF ground-state wrapper: returns energy (real) or nothing on failure
function dbf_gstate_energy(H, N)
    # For H2 neutral: 2 electrons total
    Nparticles = 2
    # initial particle ket (choose first reference)
    #ket, occ = DBF.particle_ket(N, Nparticles, 0.0; mode=:first)
    kidx = argmin([real(expectation_value(H,Ket{N}(ψi))) for ψi in 1:2^N])
    ket = Ket{N}(kidx)
    occ = occvec(ket)

    # Transform H so that the occupied bitstring becomes |000...>
    H_trans = deepcopy(H)
    for i in 1:N
        if occ[i] == 1
            H_trans = Pauli(N, X=[i]) * H_trans * Pauli(N, X=[i])
        end
    end

    # reference product state |0...0>
    ψ = Ket([0 for _ in 1:N])

    e_ref = real(expectation_value(H_trans, ψ))
    @printf("    Reference (expectation on |0..0>) = %12.8f\n", e_ref)

    # run DBF groundstate optimization inside try-catch so PES loop can continue if something fails
    try
        res = DBF.dbf_groundstate(H_trans, ψ;
                     max_iter=400, conv_thresh=1e-5,
                     evolve_coeff_thresh=1e-2,
                     grad_coeff_thresh=1e-10,
                     energy_lowering_thresh=1e-10)

        H_evolved = res["hamiltonian"]
        # compute energy after evolution (expectation with same ψ)
        e_dbf = real(expectation_value(H_evolved, ψ))
        return e_dbf
    catch err
        @printf("    DBF failed: %s\n", repr(err))
        return nothing
    end
end

# Compute PES over a vector of R distances (Å). Returns tuple of arrays (R, E_HF, E_DBF_or_nan)
function compute_pes(Rs::Vector{Float64})
    nR = length(Rs)
    E_HF_list = fill(NaN, nR)
    E_DBF_list = fill(NaN, nR)
    for (i,R) in enumerate(Rs)
        println("\n=== R = $(R) Å  (point $(i)/$(nR)) ===")
        H, N_total, E_HF = build_hamiltonian(R)
        E_HF_list[i] = E_HF

        # run DBF and get energy
        e_dbf = dbf_gstate_energy(H, N_total)
        if e_dbf === nothing
            E_DBF_list[i] = NaN
        else
            E_DBF_list[i] = e_dbf
        end

        @printf("  Summary: R=%5.3f  HF=%12.8f  DBF=%12.8f\n", R, E_HF_list[i], E_DBF_list[i])
    end
    return Rs, E_HF_list, E_DBF_list
end

# Main entry: sweep and display results
function run_pes()
    # Range of bond distances (Å)
    Rs = collect(0.5:0.1:3.5)

    Rs_out, E_HF, E_DBF = compute_pes(Rs)

    println("\nFinal PES table (R [Å], E_HF [Ha], E_DBF [Ha])")
    println("-------------------------------------------------------")
    @printf("%8s %18s %18s\n", "R (Å)", "E_HF (Ha)", "E_DBF (Ha)")
    for i in eachindex(Rs_out)
        if isnan(E_DBF[i])
            @printf("%8.3f %18.8f %18s\n", Rs_out[i], E_HF[i], "failed")
        else
            @printf("%8.3f %18.8f %18.8f\n", Rs_out[i], E_HF[i], E_DBF[i])
        end
    end
    println("\nDBF Energies only:")
    for i in eachindex(Rs_out)
            @printf("%2.8f\n", E_DBF[i])
    end

    return Rs_out, E_HF, E_DBF
end

# Run script
Rs, E_HF, E_DBF = run_pes()