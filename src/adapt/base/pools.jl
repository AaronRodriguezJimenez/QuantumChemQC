using PauliOperators


"""
  The following functions implement standard pools (GSD for fermions)
  Part of them has been taken or based from Kyle Sherbert's ADAPT.jl package:
  https://github.com/kmsherbertvt/ADAPT.jl/blob/main/src/base/pools.jl

  Qubit excitation pool is based on: https://www.nature.com/articles/s42005-021-00730-0

  # N is the number of qubits
    # i,j are the qubits indices as defined Yordanov et. al. 2021
    # Note that Yordanov's unitaries are defined as `exp(iθG)` rather than `exp(-iθG)`,
    # so variational parameters will be off by a sign.
    # Returns 
    # -PauliOperators.PauliSum{N} : The qubit excitation operator as a PauliSum
    # Note that all pauli terms in any single qubit excitation operator commute, 
    #so we can return a single PauliSum
"""
function qubitexcitation(n::Int, i::Int, k::Int)
    return 0.5 .* [Pauli(n, X=[i], Y=[k]),
                   -Pauli(n, X=[k], Y=[i])]
end

function qubitexcitation(n::Int, i::Int, j::Int, k::Int, l::Int)
        return (1/8) .* [Pauli(n; X=[i,k,l], Y=[j]),
                         Pauli(n; X=[j,k,l], Y=[i]),
                         Pauli(n; X=[l], Y=[i,j,k]),
                         Pauli(n; X=[k], Y=[i,j,l]),
                         -Pauli(n; X=[i,j,l], Y=[k]),
                         -Pauli(n; X=[i,j,k], Y=[l]),
                         -Pauli(n; X=[j], Y=[i,k,l]),
                         -Pauli(n; X=[i], Y=[j,k,l])]
    end

"""                
    qubitexcitationpool(n_system::Int)
          
    The number of singles excitations = (n 2), and the doubles = 3*(n 4).
            
   # Parameters
    - `n_system`: Number of qubits in the system

    # Returns
    - `pool`: the qubit-excitation-based pool as defined in Communications Physics 4, 1 (2021).
    - `target_and_source`: Dict mapping each pool operator to the target and source orbitals involved in the excitation. 
"""               
function qubitexcitationpool(n_system::Int)
    pool = Vector{Pauli{n_system}}()
    coeffs = Vector{Float64}()
    #target_and_source = Dict{PauliSum{n_system}, Vector{Vector{Int64}}}()
                            
    for i in 1:n_system
        for j in i+1:n_system
            # singles excitations
            singles = qubitexcitation(n_system, i, j)
            for op in singles
#                println("Operator in Single excitation")
#                display(op)
                push!(pool, op)
                push!(coeffs, 1.0) # Initial coefficient
            end
#            target_and_source[singles] = [[i,j]]

            # doubles excitations
            for k in j+1:n_system
                for l in k+1:n_system
                    target_pair = [i,j]; source_pair = [k,l]
                    doubles = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                    for op in doubles
#                        println("Operator in Double excitation")
#                        display(op)
                        push!(pool, op)
                        push!(coeffs, 1.0) # Initial coefficient
                    end
#                    arget_and_source[doubles] = [target_pair,source_pair]

                    target_pair = [i,k]; source_pair = [j,l]
                    doubles = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                    for op in doubles
#                        println("Operator in Double excitation")
#                        display(op)
                        push!(pool, op)
                        push!(coeffs, 1.0) # Initial coefficient
                    end
#                    target_and_source[doubles] = [target_pair,source_pair]
               
                    target_pair = [j,k]; source_pair = [i,l]     
                    doubles = qubitexcitation(n_system, target_pair[1], target_pair[2], source_pair[1], source_pair[2])
                    for op in doubles
#                        println("Operator in Double excitation")
#                        display(op)
                        push!(pool, op)
                        push!(coeffs, 1.0) # Initial coefficient
                    end
#                    target_and_source[doubles] = [target_pair,source_pair]

                    end
                end
            end
        end
        return pool,coeffs #, target_and_source                                            
    end



