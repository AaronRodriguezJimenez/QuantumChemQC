"""
    FermionOp(index::Int, creation::Bool)

A single fermionic operator, either creation (true) or annihilation (false).
"""
struct FermionOp
    mode::Int64       #fermionic mode 
    creation::Bool      #0-destruction 1-creation
 #   coefficient::ComplexF64  #Some associated coefficient
end

struct FermionOperator
    terms::Dict{Vector{FermionOp}, ComplexF64}
end

#Constructors
"""
 This function constructs a FermionOperator structure resembling a product of 
 fermionic operators represented as a vector with an associated amplitude
"""
function FermionOperator(term::Vector{Tuple{Int64, Bool}}, coeff::T) where T<:Number
    ops = [FermionOp(i, c) for (i, c) in term]
    return FermionOperator(Dict(ops => ComplexF64(coeff)))
end

#vaccuum case
FermionOperator(coeff::ComplexF64) = FermionOperator(Dict(Vector{FermionOp}() => coeff))

#Pretty printing
function Base.show(io::IO, op::FermionOp)
    print(io, op.creation ? "$(op.mode)^" : "$(op.mode)")
end

function Base.show(io::IO, op::FermionOperator)
    first = true
    for (ops, coeff) in op.terms
        if !first
            print(io, " + ")
        end
        print(io, coeff)
        print(io, " [")
        if !isempty(ops)
            join(io, ops, " , ")
        end
        first = false
    end
    print(io, "]")
end