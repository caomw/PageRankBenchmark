using BufferedStreams

"""
Optimized function that writes list of edges as tab seperated values
Tries to be perform as little allocation as possible
"""
function writetsv(filename, ij)
   file = BufferedOutputStream(open(filename, "w"))
   arr = Vector{UInt8}(64) # scratch space for dec!
   for (i, j) in ij
      nchars = dec!(arr, i)
      writeto(file, arr, nchars)
      write(file, '\t' % UInt8)
      nchars = dec!(arr, j)
      writeto(file, arr, nchars)
      write(file, '\n' % UInt8)
   end
   close(file)
end

function readtsv(filename)
   file = IOBuffer(Mmap.mmap(open(filename), Vector{UInt8}, (filesize(filename),)))
   ij = Vector{Tuple{Int64, Int64}}(0)
   while !eof(file)
      i = parseuntil(file, '\t' % UInt8)
      j = parseuntil(file, '\n' % UInt8)
      push!(ij, (i, j))
   end
   close(file)
   return ij
end

function writeto(io, arr, nchars)
   @inbounds for i in 1:nchars
      write(io, arr[i])
   end
end

# Assumes positive number
function dec!(arr::Vector{UInt8}, x::Int64)
   nchars = i = Base.ndigits0z(x)
   if i > length(arr) # only resize if array is to small
      resize!(arr, i)
   end
   while i > 0
      arr[i] = ('0'+rem(x,10)) % UInt8
      x ÷= 10
      i -= 1
   end
   nchars
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
      b = read(file, UInt8)
   end
   v
end

