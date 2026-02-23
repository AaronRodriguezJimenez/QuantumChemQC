
"""
 # # # Operator evolution related functions
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

"""
 Converts the PauliSum Operator into generators and angles vectors
 H :: A PauliSum Hamiltonian
"""
function gens_from_H(H::PauliSum{N,T}) where {N,T}
    generators = Vector{PauliBasis{N}}([])
    angles = Vector{Float64}([])
    for (p,c) in H
        println(c)
        push!(generators, p)
        push!(angles, c)
    end
    return generators, angles
end


"""
Quantum signal from Pauli-propagation time series with pruning controls.
Arguments:
- generators, angles  : from UnitaryPruning.heisenberg_1D(...)
- o                   : PauliSum (used as both V and W by default)
- ket                 : state for estimator (your current approach)
- thresh              : magnitude threshold (|coeff|)
- n_intervals : number of intervals,  n_intervals =tot_time/dt
"""
function QSP_evolution_op(H::PauliSum{N,T}, o::PauliSum{N,T},
                      n_intervals, 
                      dt, 
                      ket;
                      thresh::Float64=1e-3) where {N, T}

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

    #generators, angles = gens_from_H(H)

    #angles = (2*dt/trott_steps) .* angles 
    #nt = length(angles)                            
    
    # Evolve W under Trotterization
    for ki in 1:n_intervals
     #   for rep in 1:trott_steps
            # Trotter terms U_1 U_2, ..., U_k
            for (p,c) in H #1:nt
                # Access to the evolution of the operator by H = Sum(theta_i * P_i)
                #Pi  = generators[j]
                #display(Pi)
                #theta = 2*dt*angles[j]
                #pb = PauliBasis(Pi)
                theta = real(c)
                evolve!(Wt, p, theta)
            end
            # --- POST pruning ---
            coeff_clip!(Wt, thresh=thresh)
            # -----------------------------
      #  end

        WWt = W * Wt # OTOC-like product
        expval = expectation_value(WWt, ket) # Contraction with reference ket
        Ctreal = real(expval) # Real part of C(t) = <O(0) * (U_i^ O U_i)>
        Ctimag = imag(expval)
        push!(rCtvals, Ctreal)
        push!(iCtvals, Ctimag)
    end

    tgrid = collect(range(0.0, stop=n_intervals*dt, length=length(rCtvals)))

    return rCtvals, iCtvals, tgrid
end
