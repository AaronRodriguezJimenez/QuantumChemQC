# QuantumChemQC.jl

`QuantumChemQC.jl` is a Julia package for performing basic **quantum chemistry calculations** with an eye toward applications in **quantum computing**.  

Main Purpose: QuantumChemQC converts quantum-chemistry data into PauliOperators objects and provides chemistry-oriented workflows such as molecular Hamiltonian construction and correlation-function simulations. 

---

## Features

- Compute **molecular integrals** using Gaussian basis sets
- Perform **Restricted Hartree–Fock (RHF)** calculations
- Access to orbital energies, MO coefficients, Fock matrix, etc.
- Integration with [`PauliOperators.jl`](https://github.com/nmayhall/PauliOperators.jl) for efficient representation of Pauli strings and associated methods.
- Perfoms quantum signal processing through ODMD, and DMD pipelines.
- Modular structure: easily extensible for post-HF methods or quantum mapping

---

## Current Limitations

- Only supports **closed-shell (RHF)** systems

---

## Installation

In your Julia environment:

Download and unzip or clone from github

```julia
using Pkg
Pkg.develop(path="/path/to/QuantumChemQC")