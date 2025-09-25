"""
 This module Implement ADAPT-VQE with standard pools (GSD for fermions, spin-operators for spins).
"""
module ADAPT_PP

    include("core/types.jl")
    include("core/protocols.jl")

    export run_adapt!

    module base
        include("base/optimizer.jl")
        include("base/pools.jl") 

    end # module base

# re-export pools at this level for convenience
using .base: pools
export pools


end # module ADAPT_PP

"""
 Using example:

 
using QuantumChemQC

# Access through parent:
QuantumChemQC.ADAPT_PP.run_adapt!(...)

# Or directly (thanks to export ADAPT_PP):
ADAPT_PP.run_adapt!(...)
ADAPT_PP.pools(...)
"""