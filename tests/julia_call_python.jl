using PyCall
"""
    This test script demonstrates calling Python code from Julia using PyCall.
    It imports a custom Python module 'pyqctools' which wraps some PySCF functionality,
    and performs a simple SCF energy calculation for H2 molecule.
    WE ARE USING A VIRTUAL ENVIRONMENT FOR PYTHON LOCATED AT:
"""

println("PyCall.python = ", PyCall.python)   # should be /Users/admin/VSCProjects/py4julia/bin/python
#If you edit Python code and want Julia to see updates without restarting
importlib = pyimport("importlib")
importlib.reload(pyimport("pyqctools"))

# import the module (lowercase)
pq = pyimport("pyqctools")
println("pyqctools imported:", pq)

# run a tiny test with PySCF to validate round-trip mutation
@pyimport pyscf.gto as gto
@pyimport pyscf.scf as scf

mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g")
println("initial basis: ", mol.basis)
pq.set_basis_and_build(mol, "6-31g")
println("new basis: ", mol.basis)

mf = scf.RHF(mol)
energy = pq.compute_scf_energy(mf)
println("SCF energy: ", energy)



