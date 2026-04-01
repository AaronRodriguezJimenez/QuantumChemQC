using LinearAlgebra
using DBF
using PauliOperators
using Statistics
using Plots
using Printf

#"""
# Simple test for checking accuracy at ground state on a simple model.
#
#- - - Hamiltonian - - -
#
function heisenberg_1D(N, Jx, Jy, Jz; x=0, y=0, z=0)
    H = PauliSum(N, Float64)
    for i in 0:N-1
        H += -Jx * Pauli(N, X=[i+1,(i+1)%(N)+1])
        H += -Jy * Pauli(N, Y=[i+1,(i+1)%(N)+1])
        H += -Jz * Pauli(N, Z=[i+1,(i+1)%(N)+1])
    end 
    for i in 1:N
        H += x * Pauli(N, X=[i])
        H += y * Pauli(N, Y=[i])
        H += z * Pauli(N, Z=[i])
    end 
    coeff_clip!(H)
    return H
end

# Define model parameters
Jx = 1.0
Jy = 1.0
Jz = 1.0
gx = 0.0
gy = 0.0
gz = 0.0
g = 0.00#-0.01
N = 4 #Total number of qubits
H = heisenberg_1D(N, Jx, Jy, Jz)
@printf("1D-Heisenberg Hamiltonian (J= %.2f, Jz= %.2f): \n", Jx, Jz)
display(H)

# 1) diagonalize H
E = eigvals(Hmat)
E = sort(real(E))   # sorted ascending

# 2) print first few eigenvalues and consecutive gaps
Kshow = min(20, length(E))
@printf("\nFirst %d eigenvalues:\n", Kshow)
for i in 1:Kshow
    @printf("  %2d  E = %12.8f\n", i-1, E[i])
end

@printf("\nConsecutive gaps (first %d):\n", Kshow-1)
for i in 1:Kshow-1
    gap = E[i+1] - E[i]
    @printf("  gap %2d->%2d = %12.8f\n", i-1, i, gap)
end


#- - - Helpers for getting the occupation from bitstring - - - -
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
    return ket, v, out
end

# Now use DBF to compute ground state energy
function dbf_gstate(H::PauliSum{N, T}, ket::Ket{N}) where {N,T}

    # = = = DBF parameters = = = 
    max_iter=2000
    conv_thresh=1e-10
    evolve_coeff_thresh=1e-6
    grad_coeff_thresh=1e-10
    energy_lowering_thresh=1e-10
    max_rots_per_grad=50

    println(occ)
    println(N)

    println("Initial state:")
    display(ket)
    #e1 = expectation_value(H,ket)
    #@printf(" Reference = %12.8f\n", e1)

    # Transform H to make |000> the most stable bitstring
    for i in 1:N
        if occ[i] == 1
            H = Pauli(N, X=[i]) * H * Pauli(N, X=[i])
        end
    end

    H0 = deepcopy(H)
    ψ = Ket([0 for i in 1:N])
    display(ψ)

    e1 = expectation_value(H,ψ)
    @printf(" Reference = %12.8f\n", e1)

    g = Vector{PauliBasis{N}}([]) 
    θ = Vector{Float64}([]) 

    println("\n ########################")
    res = DBF.dbf_groundstate(H0, ψ,
                    max_iter= max_iter,
                    conv_thresh= conv_thresh, 
                    evolve_coeff_thresh= evolve_coeff_thresh,
                    grad_coeff_thresh= grad_coeff_thresh,
                    energy_lowering_thresh= energy_lowering_thresh,
                    max_rots_per_grad = max_rots_per_grad,
                    n_body=1
                    )                        
    
    H = res["hamiltonian"]
    #gi = res["generators"]
    #θi = res["angles"]
    e2 = real(expectation_value(H,ψ))

    #- - - Get Calculation Details
    println("\n===============================")
    println("= = = DBF Parameters = = =")
    println("max_iter : ", max_iter) 
    println("conv_thresh : ",conv_thresh)
    println("evolve_coeff_thresh :", evolve_coeff_thresh)
    println("grad_coeff_thresh : ", grad_coeff_thresh)
    println("energy_lowering_thresh : ", energy_lowering_thresh)

    println("Number of qubits: ", N)
    println("Input Hamiltonian #terms: ", length(H))

    # Initial energy
    println("\n===============================")
    println("\n Initial energy:")
    println("Index of min energy bitstring: ", kidx)
    println("Initial state:")
    display(ket)
    println("Occuation vector: ", occ)
    @printf("<H> = %12.8f <U'HU> = %12.8f \n", e1, e2)

                    
    return
end

ket, occ, kidx  = string_to_ket("0000") #Leading CAS/sto3g configuration SINGLET
println("Initial state:")
display(ket)
e1 = expectation_value(H,ket)
@printf(" Reference = %12.8f\n", e1)
dbf_gstate(H, ket)