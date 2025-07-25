struct DIIS
    trial_vector::Array{Float64,3}
    residual_vector::Array{Float64,3}
    L::Matrix{Float64}
    depth::Int
    cycle::MVector{1,Int} #Mutable statically sized vector altenative is Ref{Int}
    function DIIS(n::Integer, depth::Integer=10)
        trial_vector = zeros(n, n, depth)
        residual_vector =  zeros(n, n, depth)
        L = zeros(depth+1, depth+1)
        L[1, 2:end] .= 1
        L[2:end, 1] .= 1
        new(trial_vector, residual_vector, L, depth, [0])
    end
end

function update!(
    diis::DIIS,
    trial::AbstractMatrix{Float64},
    residual::AbstractMatrix{Float64},
    )
    ptr = mod(diis.cycle[], diis.depth) + 1
    diis.trial_vector[:, :, ptr] = trial
    diis.residual_vector[:, :, ptr] = residual

    n = diis.cycle[] < diis.depth ? ptr : diis.depth
    #
    #for i = 1:n
    #    diis.L[i+1, ptr+1] = diis.residual_vector[:, :, i] ⋅ residual
    #    diis.L[ptr+1, i+1] = diis.L[i+1, ptr+1]
    #end
    #
    #coeff = inv(diis.L[1:n+1, 1:n+1])[2:n+1, 1]
    #
    #REPLACED BY (JULY 2025):
    # inv() is trying to invert the (n+1)x(n+1) DIIS B matrix, which becomes
    # singular in some cases when the residual vectors used to build it are
    # linearly depeendent or poorly conditioned, instead, we do the folllowing
    # The \ operator avoids full matrix inversion and is numberically safer,
    # it uses LU with pivoting, then:
    #- - - - - - - - - - - -
    # Fill the DIIS matrix L
    diis.L[1, 1] = 0.0
    
    for i = 1:n
        diis.L[1, i+1] = -1.0
        diis.L[i+1, 1] = -1.0
    end

    Lsub = diis.L[1:n+1, 1:n+1]
    rhs = zeros(n+1)
    rhs[1] = -1.0

    # Check condition number or wrap in try block
    coeff_full = try
        Lsub \ rhs
    catch e
        @warn "Singular DIIS matrix. Skipping update."
        return
    end
    
    coeff = coeff_full[2:end]
    ##########

    m = LinearAlgebra.checksquare(trial)
    for i = 1:m, j = i:m
        trial[i, j] = 0.0
        for k = 1:n
            trial[i, j] += coeff[k] * diis.trial_vector[i, j, k]
        end
        trial[j, i] = trial[i, j]
    end
    diis.cycle[] += 1
end