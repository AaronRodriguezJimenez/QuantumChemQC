# inspect_file.jl
using Printf

function inspect_file(path; nbytes=512)
    bytes = open(read, path)
    n = min(length(bytes), nbytes)
    b = bytes[1:n]

    # Hex string (format each byte individually)
    hex_tokens = Vector{String}(undef, n)
    for i in 1:n
        hex_tokens[i] = @sprintf("%02x", b[i])
    end
    println("Hex (first $n bytes):")
    println(join(hex_tokens, " "), "\n")

    # Printable preview: printable ASCII shows, others become '.'
    preview_chars = Char[]
    for x in b
        push!(preview_chars, (0x20 <= x <= 0x7e) ? Char(x) : '.')
    end
    println("Printable preview:\n", String(preview_chars), "\n")

    println("Contains NUL (0x00)? ", any(x -> x == 0x00, b))
    println("Total file size (bytes): ", length(bytes))
end

inspect_file("output_H2_signal.txt", nbytes=512)