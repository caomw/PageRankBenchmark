# an in-place equivalent to KronGraph500NoPerm in the MATLAB implementation
function generate_edges!(edges, rand_array, scale)
    # set R-MAT (2x2 Kronecker) coefficients
    a = 0.57
    b = 0.19
    c = 0.19
    d = 1.0 - (a + b + c)

    # calculate normalization coefficients
    a_plus_b = a + b
    c_norm = c / (1.0 - a_plus_b)
    a_norm = a / a_plus_b

    # loop over each scale
    @inbounds for j in 1:length(edges)
        rand!(rand_array)
        start_node = one(Int)
        end_node = one(Int)
        @inbounds for i in 1:scale
            k = 1 << (i - 1)
            start_bit = rand_array[i] > a_plus_b
            end_bit = rand_array[i + scale] > ifelse(start_bit, c_norm, a_norm)
            start_node += k * start_bit
            end_node += k * end_bit
        end
        edges[j] = (start_node, end_node)
    end

    return edges
end

function kernel0(path, nfiles, scale, edges_per_vertex)
    edges = Vector{Tuple{Int,Int}}(k0_edges_per_file(nfiles, scale, edges_per_vertex))
    rand_array = Vector{Float64}(2 * scale)
    k0path = kernel_path(path, 0)
    buffer = Vector{UInt8}(64)
    for file_index in my_file_range(nfiles)
        srand(file_index) # uniquely seed RNG specifically for this file
        generate_edges!(edges, rand_array, scale)
        open(joinpath(k0path, "$(file_index).tsv"), "w") do file
            write_edges!(file, edges, buffer)
        end
    end
end
