
"""
  G Fock Matrix Eris, P density matrix, I are the computed ERIs
"""
function two_e_fock_matrix!(G::Matrix{Float64}, P::Matrix{Float64}, I::AbstractArray{Float64,4})
    n = LinearAlgebra.checksquare(G)
    for i=1:n, j=i:n
        G[i,j] = 0.0
        for k=1:n, l=1:n
            G[i,j] += (2 * I[i,j,k,l] - I[i,k,j,l]) * P[k,l]
        end
        G[j,i] = G[i,j]
    end
    return nothing
end



function two_e_fock_matrix(P::Matrix{Float64}, I::AbstractArray{Float64, 4})
    G = similar(P)
    two_e_fock_matrix!(G,P,I)
    return G
end

function density_matrix!(P::Matrix{Float64}, C::Matrix{Float64}, N::Int)
    n = LinearAlgebra.checksquare(P)
    for i=1:n, j=1:n
        P[i,j] = 0.0
        for a=1:N
            P[i,j] += C[i,a] * C[j,a]
        end
    end
    return nothing
end

function density_matrix(C::Matrix{Float64}, N::Int)
    P = similar(C)
    density_matrix!(P, C, N)
    return P
end

function scf(
    S::Matrix{Float64},
    T::Matrix{Float64},
    V::Matrix{Float64},
    I::AbstractArray{Float64, 4},
    P0::Matrix{Float64},
    N::Int,
    e_tol::Float64=1e-7,
    d_tol::Float64=1e-7,
    max_cycle::Int=64
    )
    
    n = size(S,1)
    Hcore = T + V

    #Guess Fock matrix
    P = copy(P0)
    G = two_e_fock_matrix(P, I)
    F = Hcore + G

    #Guess electron energy
    Eel = sum(P .* (Hcore + F)) #Frobenius inner product
    println("Cycle 0: Eel = $(Eel)")
    #@info "Cycle 0: Eel = $(Eel)"

    diis = DIIS(n) # Call diis function defined in diis.jl

    for cycle = 1:max_cycle
        #Update density matrix
        e::Vector{Float64}, C::Matrix{Float64} = eigen(F,S)
        density_matrix!(P, C, N)

        #Update Fock Matrix
        two_e_fock_matrix!(G, P, I)
        F .= Hcore .+ G

        #Calculate electron energy
        Eold = Eel
        Eel = sum(P .* (Hcore + F))

        #Calculate DIIS error
        SPF = S * P * F
        residual = SPF' - SPF
        diis_err = norm(residual) / sqrt(length(residual))

        #Test convergence
        if abs(Eel - Eold) < e_tol && diis_err < d_tol
            println("Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err) Converged!")
            #@info "Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err) Converged!"
            return Eel, P, e, C
        else
            println("Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err)")
            #@info "Cycle $(cycle): Eel = $(Eel), EDelta = $(Eold - Eel), DIIS Error = $(diis_err)"
        end

        # Build DIIS Fock matrix
        update!(diis, F, residual)
    end

    @error "SCF failed after $(max_cycle) cycles"
end

struct SCF
    mol::String
    S::Matrix{Float64}
    T::Matrix{Float64}
    V::Matrix{Float64}
    I::Array{Float64, 4}
    P::Matrix{Float64}
    C::Matrix{Float64}
    e::Vector{Float64}
    Eel::Float64
    Enuc::Float64
end

# SCF  constructors
function SCF(mol_str, bset, p)
     # Compute integrals
    S, T, V, _, I = QuantumChemQC.compute_integrals(bset, p)

    # Initial density guess
    P0 = zeros(size(S))
    
    # Number of doubly occupied orbitals (closed-shell HF)
    N = div(p.nalpha + p.nbeta, 2)

    # Perform SCF
    Eel, P, e, C = scf(S, T, V, I, P0, N)

    # Nuclear repulsion energy
    Enuc = QuantumChemQC.nuclear_repulsion_energy(bset, p)

    # Return SCF object
    return SCF(mol_str, S, T, V, I, P, C, e, Eel, Enuc)
end


function SCF(mol_str::String; basis="sto-3g", spherical=true)
    bset, p = QuantumChemQC.molecule(mol_str, basis; spherical)
    
    # Compute integrals
    S, T, V, _, I = QuantumChemQC.compute_integrals(bset, p)

    # Initial density guess
    P0 = zeros(size(S))
    
    # Number of doubly occupied orbitals (closed-shell HF)
    N = div(p.nalpha + p.nbeta, 2)

    # Perform SCF
    Eel, P, e, C = scf(S, T, V, I, P0, N)

    # Nuclear repulsion energy
    Enuc = QuantumChemQC.nuclear_repulsion_energy(bset, p)

    # Return SCF object
    return SCF(mol_str, S, T, V, I, P, C, e, Eel, Enuc)
end

function total_energy(scf::SCF)
    return scf.Eel + scf.Enuc
end

function electron_energy(
    S::Matrix{Float64},
    T::Matrix{Float64},
    V::Matrix{Float64},
    I::AbstractArray{Float64,4},
    P0::Matrix{Float64},
    N::Int,
    )
    Eel, P, e, C = scf(S, T, V, I, P0, N)
    return Eel
end

"""
 Compute molecular integrals from atomic integrals
 Transform Atomic integrals into Molecular integrals
 two-body coefficients are comptued following the convention from OpenFermion:
 h[p,q,r,s] = (ps|qr)
"""
function ao2mo_coefficients(C::Matrix{Float64}, ao_hcore::Matrix{Float64}, ao_eris::AbstractArray{Float64,4})
    #C = scf_obj.C # MO coeffs
    mo_hcore = C' * ao_hcore * C 
    mo_eris = Array{Float64, 4}(undef, size(ao_eris)...)

    # Openfermion convention:
    @tullio mo_eris[p,q,r,s] := C[μ, p] * C[λ, s] * C[ν, q] * C[σ, r] * ao_eris[μ, λ, ν, σ]
    return mo_hcore, mo_eris
end