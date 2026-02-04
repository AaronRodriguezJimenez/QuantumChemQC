using NPZ

H = npzread("/Users/admin/VSCProjects/QuantumChemQC/tests/1.5_tensors.npz")

println("------------------------------------")
H0 = H["hc"][1]
println("H0 ", H0)
H1 = H["h1e"]
println(typeof(H0))
println("------------------------------------")
println(length(H1))
#H2 = H["h2e"]
#println(H2)
#println("------------------------------------")