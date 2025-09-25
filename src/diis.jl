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

#=
#  Original update! function using inv() 
#  We can make the DIIS update more robust by building the B matrix from
#  flattened residuals and solving via SVD + pseudoinverse with thresholding
#  and Tikhonov regularization fallback. See notes below.

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
=#

function update!(
    diis::DIIS,
    trial::AbstractMatrix{Float64},
    residual::AbstractMatrix{Float64},
)
    # store new trial & residual in circular buffer
    ptr = mod(diis.cycle[], diis.depth) + 1
    diis.trial_vector[:, :, ptr] = trial
    diis.residual_vector[:, :, ptr] = residual

    # number of stored vectors (n)
    n = diis.cycle[] < diis.depth ? ptr : diis.depth

    # Build matrix of flattened residuals as columns (m x n)
    # where m = number of independent matrix elements (upper triangle)
    # We will vectorize using full matrix flattened; if symmetric and storage is upper-triangle,
    # adapt accordingly. Here we flatten the full matrix (works always).
    m = length(vec(residual))
    R = zeros(m, n)
    for i = 1:n
        R[:, i] = vec(diis.residual_vector[:, :, i])
    end

    # Optional: prune nearly linearly dependent residuals using QR
    # Compute thin QR to detect small diagonal R entries => dependent columns
    Q, R_q = qr(R, Val(true))  # economy QR: R_q is upper triangular
    diagR = abs.(diag(R_q))
    tol_prune = maximum(size(R)) * eps(Float64) * maximum(diagR)
    keep_idx = findall(x -> x > tol_prune, diagR)
    if length(keep_idx) < n
        @debug "Pruned $(n - length(keep_idx)) nearly dependent DIIS residual(s)."
        R = R[:, keep_idx]
        n = size(R, 2)
        # Also reorder the trial & residual buffers to match keep_idx:
        temp_trial = zeros(size(trial,1), size(trial,2), n)
        temp_resid  = zeros(size(residual,1), size(residual,2), n)
        for (k, idx) in enumerate(keep_idx)
            temp_trial[:, :, k] = diis.trial_vector[:, :, idx]
            temp_resid[:, :, k]  = diis.residual_vector[:, :, idx]
        end
        diis.trial_vector[:, :, 1:n] = temp_trial
        diis.residual_vector[:, :, 1:n] = temp_resid
    end

    # Build the (n x n) overlap matrix B where B_ij = <r_i | r_j>
    B = zeros(n, n)
    for i = 1:n
        for j = i:n
            B[i, j] = dot(R[:, i], R[:, j])
            B[j, i] = B[i, j]
        end
    end

    # Build augmented (n+1)x(n+1) matrix Lsub and rhs
    Lsub = zeros(n+1, n+1)
    Lsub[2:end, 2:end] = B
    Lsub[1, 2:end] .= -1.0
    Lsub[2:end, 1] .= -1.0
    rhs = zeros(n+1)
    rhs[1] = -1.0

    # Solve robustly via SVD + pseudo-inverse with small-σ threshold + optional Tikhonov
    # If SVD indicates near-singular matrix, we regularize using small reg*I and retry.
    svd_tol_factor = maximum(size(Lsub)) * eps(Float64)
    coeff_full = nothing
    try
        S = svd(Lsub)  # returns U, S, V
        sigma = S.S
        tol = svd_tol_factor * maximum(sigma)
        # form pseudoinverse of Lsub with thresholding
        Sinv = Diagonal(map(s -> s > tol ? 1.0/s : 0.0, sigma))
        Lpinv = S.V * Sinv * S.U'
        coeff_full = Lpinv * rhs

        # If solution has NaNs (extremely ill-conditioned), fallback to Tikhonov
        if any(isnan, coeff_full) || any(isinf, coeff_full)
            throw(ErrorException("pseudo-inverse produced NaN/Inf"))
        end
    catch e
        # fallback: Tikhonov regularization and solve via \ (LU)
        @warn "DIIS matrix ill-conditioned (SVD/pinv failed). Applying Tikhonov regularization."
        # choose reg scale proportional to trace of Lsub
        reg = 1e-8 * max(1.0, trace(abs.(Lsub))/n)
        success = false
        for attempt = 1:5
            try
                Lreg = copy(Lsub)
                Lreg += reg * I # regularize diagonal
                coeff_full = Lreg \ rhs
                success = true
                break
            catch
                reg *= 10.0
            end
        end
        if !success
            @warn "Singular DIIS matrix after regularization attempts. Skipping DIIS update."
            return
        end
    end

    # extract coefficients for residual/trial combination
    coeff = coeff_full[2:end]

    # form the new trial (symmetric) from the stored trial vectors
    mtr = LinearAlgebra.checksquare(trial)
    for i = 1:mtr, j = i:mtr
        val = 0.0
        for k = 1:length(coeff)
            val += coeff[k] * diis.trial_vector[i, j, k]
        end
        trial[i, j] = val
        trial[j, i] = val
    end

    diis.cycle[] += 1
end

#=
Notes:
- DIIS struct holds the trial and residual vectors, the B matrix (L),
  the depth, and a mutable cycle counter.
- update! function updates the DIIS object with new trial and residual matrices,
  computes the new trial matrix using DIIS extrapolation.
- The \ operator is used to solve the linear system instead of direct inversion
  to improve numerical stability.
- The code includes error handling to manage cases where the DIIS matrix
    becomes singular.

HOWEVER THIS IMPLEMENTATION STILL CAN BE REFINED, THE FOLLOWING ARE NOTES FOR THE
CURRENT update! FUNCTION:

EXPLANATION: 
Flatten residuals and build B from dot products:
 this is the canonical DIIS formulation (B_ij = <r_i|r_j>). 
 Using flattened residuals is robust and consistent.

Prune near-dependent residuals: QR identifies columns with tiny R_ii. 
Those indicate linear dependence; drop them before building B.
 This prevents B from being rank deficient.

Solve via SVD/pseudoinverse first: SVD gives stable pseudo-inverse 
and allows us to threshold tiny singular values (effectively performing 
Tikhonov implicitly on the small subspace).

Tikhonov fallback: if pinv yields NaN/Inf or is otherwise unstable,
 we add a tiny diagonal reg*I and try solving with \. 
 If that fails, we increase reg progressively (robust).

If everything fails, skip update: better to keep the previous iterate than to blow up.

Keep the trial construction symmetric.

ADDITINAL STRATEGIES TO CONSIDER:
These can be used alongside the robust DIIS above:

Adaptive damping: apply mixing F_new = (1-α) F_old + α F_DIIS with α starting small (0.1)
 and increasing if iterations progress smoothly.

Level shifting: add a positive shift to virtual–virtual Fock blocks (e.g. 0.5–1.0 a.u.) 
for early iterations to avoid small-gap oscillations.

DIIS on Fock vs densities: try both. For some problems DIIS on the Fock matrix converges 
better, for others DIIS on density residuals is better.

Dynamic DIIS depth: reduce depth if numerical issues appear (e.g. limit to 5 when instability appears).

Restart DIIS after failure: if n falls below 2 or DIIS skipped multiple times, 
clear buffers and continue with damping for a few iterations before re-enabling DIIS.

Precondition residuals: scale residual vectors by 1/(ε_i + ε_j - 2*ɛF) 
approximate denominators (like common orbital energy preconditioning)
— helps when near-degenerate orbitals exist.

Fractional occupations / smearing: for metallic/near-metallic systems, 
use finite-temperature smearing to stabilize occupations.

Orthogonalize residuals: you already prune via QR; performing an explicit
 orthonormalization (Gram-Schmidt) and using orthonormal residuals to 
 build B helps conditioning.

=# 