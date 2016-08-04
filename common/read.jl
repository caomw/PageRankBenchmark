"""
    read_edges(filename, delim, linesep)

Reads in a list of edges given by tuples of vertices. The vertices are separated by `delim`
and the edges are seperated by `linesep`. The input file is MMap'd to give the best speed and
alloctions are avoided where possible.

# Asumptions:
- All vertices are positive numbers

# Arguments:
- `delim :: Char`: defaults to '\t'
- `linesep :: Char`: defaults to '\n'
"""
function read_edges(filename, delim = '\t', linesep = '\n')
   read_edges(Tuple{Int64, Int64}, filename, delim, linesep)
end

function read_edges(:: Type{Tuple{Int64, Int64}}, filename :: String, delim = '\t', linesep = '\n')
   data = Mmap.mmap(open(filename), Vector{UInt8}, (filesize(filename),))
   file = IOBuffer(data)
   ij = Vector{Tuple{Int64, Int64}}(0)
   read_edges!(ij, file, UInt8(delim), UInt8(linesep))
   close(file)
   finalize(data)
   data = nothing
   return ij
end

function read_edges(:: Type{Int64}, filename, delim ='\t', linesep = `\n`)
   file = IOBuffer(Mmap.mmap(open(filename), Vector{UInt8}, (filesize(filename),)))
   I = Vector{Int64}(0)
   J = Vector{Int64}(0)
   read_edges!(I, J, file, UInt8(delim), UInt8(linesep))
   close(file)
   return (I, J)
end

"""
    read_edges!(ij, file, delim, linesep)

Expert version of `read_edges`.

# Arguments:
- `ij :: Vector{Tuple{Int, Int}}` The reading edges are appended to this vector
- `file :: IO` The caller is responsible for opening and closing the underlying IO
- `delim :: UInt8`
- `linesep :: UInt8`
"""
function read_edges!(ij, file :: IO, delim :: UInt8, linesep :: UInt8)
   while !eof(file)
      i = parseuntil(file, delim)
      j = parseuntil(file, linesep)
      push!(ij, (i, j))
   end
end

function read_edges!(I, J, file :: IO, delim :: UInt8, linesep :: UInt8)
   while !eof(file)
      i = parseuntil(file, delim)
      j = parseuntil(file, linesep)
      push!(I, i)
      push!(J, i)
   end
end

# Exploit the domain knowledge we have
# Only positive numbers
# Taken from CSV.jl
function parseuntil(file, delim::UInt8)
   v = 0
   b = read(file, UInt8)
   while b != delim
      v *= 10
      v += b - ('0' % UInt8)
      eof(file) && break
      b = read(file, UInt8)
   end
   v
end

