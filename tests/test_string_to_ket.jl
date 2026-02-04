using PauliOperators

# return a PauliOperators Ket equivalent to a given bitstring
function string_to_ket(bits::String)
    b = collect(bits)
    v = parse.(Int128, b)
    N = length(v)
    out = 0
    count = 0

    for bit in v
        if bit%2 == 1
            out += 2^count
        end
        count +=1
    end
    ket = Ket{N}(out)
    return ket, out
end

A = "1011101000"
k, o = string_to_ket(A)

println("Original bitstring")
println(A)
println("- - - - - - - - - ")
println("Resulting Ket")
display(k)
println("With associated decimal number:", o)
