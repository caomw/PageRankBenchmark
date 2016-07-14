function kernel1(path, nfiles, scale, edges_per_vertex)
    numprocs = MPI.Comm_size(MPI.COMM_WORLD)
    files = my_file_range(nfiles, numprocs)
    k0path = kernel_path(path, 0)
    k1path = kernel_path(path, 1)

    # read in this process's edges
    nedges = k0_edges_per_file(nfiles, scale, edges_per_vertex) * length(files)
    edges = Vector{Tuple{Int,Int}}(nedges)
    offset = 0
    for file_index in files
        offset = open(joinpath(k0path, "$(file_index).tsv"), "r") do file
            read_edges!(file, edges, offset)
        end
    end

    # sort this process's edges
    sort_edges!(edges)

    # calculate the splitters by having every process sample its edges
    splitters = select_splitters(numprocs, edges)

    # redistribute edges based on the splitters
    balanced_edges = redistribute_edges(splitters, edges)

    # write out this process's edges
    buffer = Vector{UInt8}(64)
    nedges, nfiles = length(balanced_edges), length(files)
    file_counter = 1
    for file_index in files
        i = div((file_counter - 1) * nedges, nfiles) + 1
        j = div(file_counter * nedges, nfiles)
        open(joinpath(k1path, "$(file_index).tsv"), "w") do file
            write_edges!(file, view(balanced_edges, i:j), buffer)
        end
        file_counter += 1
    end
end

function select_splitters(numprocs, edges)
    local_sample = map(first, evenly_spaced_sample(numprocs - 1, edges))
    global_sample = MPI.Gather(local_sample, 0, MPI.COMM_WORLD)
    if am_i_root_process()
        splitters = evenly_spaced_sample(numprocs - 1, sort!(global_sample))
        prune_splitters!(splitters)
    else
        splitters = similar(local_sample)
    end
    MPI.Bcast!(splitters, 0, MPI.COMM_WORLD)
    return splitters
end


# assumes `edges` is already sorted, that `splitters` has monotonically increasing elements
# and that length(splitters) == (numprocs - 1)
function redistribute_edges(splitters, edges)
    procs = 0:length(splitters)
    upper_bounds = vcat(splitters, typemax(Int))
    i = 1
    for proc in procs
        # Retrieve the region of edges where `lower <= first(edge) < upper`.
        # This relies on `edges` already being sorted, and `splitters` being
        # monotonically increasing.
        upper_bound = upper_bounds[proc + 1]
        j = i + 1
        while (j < length(edges)) && (first(edges[j + 1]) < upper_bound)
            j += 1
        end
        # send the process its edges
        MPI.Isend(edges[i:j], proc, 0, MPI.COMM_WORLD)
        i = j + 1
    end

    # retrieve the local process's edges from the other processes
    received_edges = similar(edges, 0)
    for proc in procs
        # make sure to divide by two, since Get_count will count the number of nodes, not
        # edges, sent to this process (i.e. it counts how many integers there are total
        # without considering that the elements are actually 2-tuples)
        received_count = div(MPI.Get_count(MPI.Probe(proc, 0, MPI.COMM_WORLD), Int), 2)
        edges_from_proc = Vector{Tuple{Int,Int}}(received_count)
        MPI.Irecv!(edges_from_proc, proc, 0, MPI.COMM_WORLD)
        append!(received_edges, edges_from_proc)
    end

    # Instead of the above, it would be better to preallocate the full array and then use
    # views as our recieve buffers. The commented-out code below implements this, but MPI.jl
    # doesn't allow views for reception buffers :(
    #
    # received_edge_counts = [MPI.Get_count(MPI.Probe(proc, 0, MPI.COMM_WORLD), Int) for proc in procs]
    # received_edges = similar(edges, sum(received_edge_counts))
    # offset = 1
    # for proc in procs
    #     count = received_edge_counts[proc + 1]
    #     MPI.Irecv!(view(received_edges, offset:count), proc, 0, MPI.COMM_WORLD)
    #     offset += count
    # end

    # chunks from each process are already sorted w.r.t. themselves, but are not sorted
    # w.r.t. each other, so we still need to do one final sort here.
    sort_edges!(received_edges)

    return received_edges
end

# Sort edges in-place by starting node. Using the `by` keyword argument will unfortunately
# make this slower than it should be until JuliaLang/julia#16580 gets merged.
sort_edges!(edges) = sort!(edges, by = first)

# Sample `n` evenly spaced elements from `elements`.
evenly_spaced_sample(n, elements) = elements[(1:n) .* div(length(elements), n + 1)]

# ensure that splitters are all suitable upper bounds for a given region
# and fix any potential overlaps
function prune_splitters!(splitters)
    if first(splitters) == 1
        splitters[1] = 2
    end
    for i in 2:length(splitters)
        if splitters[i - 1] == splitters[i]
            splitters[i] += 1
        end
    end
    return splitters
end
