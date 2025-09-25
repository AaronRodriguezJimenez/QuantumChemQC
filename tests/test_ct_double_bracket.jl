"""
 Here we test the idea for the computation of the Correlation function
    C(t) = <psi|exp(-iHt)|psi> using double bracket notation.
    psi is a HF state in the qubit bitstring notation.
    H is the qubit Hamiltonian in PauliOperator format.
    t is a time parameter.

    It can be shown that:
    C(t) = <psi|exp(-iHt)|psi> = exp(-iϕt) + 2(<0|[A†, [exp(-iH't),A]]|0>)
    
    where ϕ is the eigenvalue from H|0> = ϕ|0>
"""

