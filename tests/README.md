# QuantumChemQC.jl

`QuantumChemQC.jl` is a Julia package for performing basic **quantum chemistry calculations** with an eye toward applications in **quantum computing**.  
The current implementation is based on ideas and tools from [HartreeFock.jl](https://github.com/panxl/HartreeFock.jl.git) and [DoNOF.jl](https://github.com/felipelewyee/DoNOF.jl.git).

---

## Features

- Compute **molecular integrals** using Gaussian basis sets
- Perform **Restricted Hartree–Fock (RHF)** calculations
- Access to orbital energies, MO coefficients, Fock matrix, etc.
- Integration with packages like [`PauliOperators.jl`](https://github.com/QuantumBFS/PauliOperators.jl) for qubit Hamiltonians
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