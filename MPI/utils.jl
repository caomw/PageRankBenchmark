#################
# MPI Utilities #
#################

am_i_root_process() = MPI.Comm_rank(MPI.COMM_WORLD) == 0

#################################
# derived experiment parameters #
#################################

max_node(scale) = 2^scale

function k0_edges_per_file(nfiles, scale, edges_per_vertex)
    return div(edges_per_vertex * max_node(scale), nfiles)
end

function my_file_range(nfiles, numprocs = MPI.Comm_size(MPI.COMM_WORLD),
                       rank = MPI.Comm_rank(MPI.COMM_WORLD))
    start_file = ceil(Int, 1 + (nfiles / numprocs) * rank)
    end_file = ceil(Int, (nfiles / numprocs) * (rank + 1))
    return start_file:end_file
end

kernel_path(path, i) = joinpath(path, string("kernel", i))

#########################
# IO + string wrangling #
#########################

const CHAR_TAB = UInt8('\t')
const CHAR_NEWLINE = UInt8('\n')
const CHAR_ZERO = UInt8('0')

# writing edges #
#---------------#

function write_edges!(io, edges, buffer = Vector{UInt8}(64))
    for (start_node, end_node) in edges
        write_node!(io, buffer, int2str!(buffer, start_node))
        write(io, CHAR_TAB)
        write_node!(io, buffer, int2str!(buffer, end_node))
        write(io, CHAR_NEWLINE)
    end
end

function write_node!(io, buffer, nchars)
    @inbounds for i in 1:nchars
        write(io, buffer[i])
    end
end

function int2str!(buffer, x)
    nchars = i = Base.ndigits0z(x)
    i > length(buffer) && resize!(buffer, i)
    while i > 0
        buffer[i] = CHAR_ZERO + ((x % 10) % UInt8)
        x = div(x, 10)
        i -= 1
    end
    return nchars
end

# reading edges #
#---------------#

function read_edges!(io, edges, offset)
    i = offset
    while !(eof(io))
        i += 1
        start_node = parse_node(io, CHAR_TAB)
        end_node = parse_node(io, CHAR_NEWLINE)
        edges[i] = (start_node, end_node)
    end
    return i
end

function parse_node(io, delim)
   node = 0
   char = read(io, UInt8)
   while char != delim
      node *= 10
      node += char - CHAR_ZERO
      char = read(io, UInt8)
   end
   return node
end
