using QuantumChemQC
using PauliOperators
using LinearAlgebra
using Printf
using Plots
using NPZ

"""
 Read and create a dipole moment operator using PauliOperators
"""
function run_dip()
    # Get precomputed active space spinorbitals tensors
    data_path = "/Users/admin/PycharmProjects/pyQCTools/QSP/benzene_and_acenes/benz-RHF_dip_mo.npz"
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
    data_path =  "/Users/admin/PycharmProjects/pyQCTools/QSP/benzene_and_acenes/benz-RHF_integrals.npz"
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

# return a PauliOperators Ket equivalent to a given bitstring
function string_to_ket(bits::String)
    b = collect(bits)
    v = parse.(Int128, b)
    N = length(v)
    out = 0
    count = 0

    for bit in v
        if bit%2 == 1
            out += 2^count
        end
        count +=1
    end
    ket = Ket{N}(out)
    return ket, out
end

"""
 Evolve function from DBF code
"""
function evolve!(O::PauliSum{N, T}, G::PauliBasis{N}, θ::Real) where {N,T}
    _cos = cos(θ)
    _sin = 1im*sin(θ)
    sin_branch = PauliSum(N)
    for (p,c) in O
        if PauliOperators.commute(p,G) == false
            # replace sum! with more efficient version
            # sum!(sin_branch, c*_sin*G*p)
            tmp = c*_sin*G*p
            curr = get(sin_branch, PauliBasis(tmp), 0.0) + PauliOperators.coeff(tmp)
            sin_branch[PauliBasis(tmp)] = curr 
            O[p] *= _cos
        end
    end
    sum!(O, sin_branch)
    return O 
end

# Helper: extract PauliBasis operators and coefficients from PauliSum
function extract_hamiltonian_coeffs_and_ops(H::PauliSum{N,T}) where {N,T}
    ops = PauliBasis{N}[]
    coeffs = Float64[]
    for (p, c) in H
        push!(ops, p)
        push!(coeffs, float(c))
    end
    return ops, coeffs
end

"""
Pauli-propagation time series with pruning controls.
Arguments:
- generators, angles  : from UnitaryPruning.heisenberg_1D(...)
- o                   : PauliSum (used as both V and W by default)
- ket                 : state for estimator (your current approach)
- thresh              : magnitude threshold (|coeff|)
- n_intervals : number of intervals,  n_intervals =tot_time/dt
"""
function evolution_op(ket, o::PauliSum{N,T}, H::PauliSum{N,T}, n_intervals, dt;
                      thresh::Float64=1e-3) where {N,T}

    #err = Vector{Float64}([])
    Wt = deepcopy(o)        # evolve W ≡ U*OU
    W  = deepcopy(o)        # initial operator 
    rCtvals = Vector{Float64}([])#(undef, nsamp) # vector for C(t) values storing
    iCtvals = Vector{Float64}([])#(undef, nsamp)

    # t = 0 
    WW = W*W
    expval0 = expectation_value(WW, ket)

    C0real = real(expval0)
    C0imag = imag(expval0)

    push!(rCtvals, C0real)
    push!(iCtvals, C0imag)

    # Extract Pauli strings and coefficients
    generators, angles = extract_hamiltonian_coeffs_and_ops(H)
    #@printf("Hamiltonian has %d terms \n", length(coeffs))

    nt = length(angles)                            
    println("Total Rotations:", nt * n_intervals)
    # Evolve W under Trotterization

    #accumulated_error = 0
    for ki in 1:n_intervals
            # Trotter terms U_1 U_2, ..., U_k
            
            WWt = W * Wt # OTOC-like product
            accumulated_error = 0
            e1 = expectation_value(WWt,ket)

            for j in 1:nt
                # Access to the evolution of the operator by H = Sum(theta_i * P_i)
                Pi  = generators[j]
                theta = 2*dt*angles[j]
                pb = PauliBasis(Pi)
                
                evolve!(Wt, pb, theta)

                # --- POST pruning ---
                coeff_clip!(Wt, thresh=thresh)
                WWt = W * Wt 
                e2 = expectation_value(WWt, ket)
                accumulated_error += e2-e1
            end           
            
            #WWt = W * Wt (at time interval)
            expval = expectation_value(WWt, ket) - accumulated_error # Contraction with reference ket
            Ctreal = real(expval) # Real part of C(t) = <O(0) * (U_i^ O U_i)>
            Ctimag = imag(expval)
            push!(rCtvals, Ctreal)
            push!(iCtvals, Ctimag)
    end

    tgrid = collect(range(0.0, stop=n_intervals*dt, length=length(rCtvals)))

    return rCtvals, iCtvals, tgrid
end


D = run_dip()
H = run_H()
ket, _ = QuantumChemQC.string_to_ket("111111000000")

n_intervals = 100
t = 50.0 #Total Time Evolution
dt = t/n_intervals
threshold = 1e-4 #pruning threshold based on coeff.


t1 = time()
rRES, iRES, tgrid = evolution_op(ket, D, H, n_intervals, dt; thresh=threshold)

                                  # Code block to measure
elapsed_time = time() - t1
println("Elapsed time: ", elapsed_time, " seconds")

# Number of snapshots actually returned
nsnap = length(rRES)
println("* * * * Number of snapshots collected: $nsnap")

# Print C(t) results
plt = plot(tgrid, rRES, lw=2, seriestype=:scatter,
          label="Re(C(t), th=$threshold")
plt = plot!(tgrid, iRES, lw=2, seriestype=:scatter,
          label="Im(C(t), th=$threshold")

xlabel!(plt, "Time"); ylabel!(plt, "< O(0)O(t) >")
title!(plt, "Benzene (6e,6o)")
savefig(plt, "/Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_benzene_PP_thesh_$threshold.pdf")

println("- - - <0|H|0> - - - -")
function compute_ref_expval(H, psi0)
    ref_expval = 0
    for (p, c) in H
        ref_expval += c * expectation_value(p, psi0)
    end
    return ref_expval
end

ref_expval = compute_ref_expval(H, ket)
println(ref_expval)

open("Users/admin/VSCProjects/QuantumChemQC/QSP/QSP_benzene_PP_thesh_$threshold.txt", "w") do io
    println(io, "- - - Output signals  - - - ")
    @printf(io, "dt   Re(C(t))    Im(C(t))\n")
    for (i, interval) in enumerate(tgrid)
        @printf(io, "%.4f     %.6f    %.6f\n",
                interval, rRES[i], iRES[i])
    end
end

