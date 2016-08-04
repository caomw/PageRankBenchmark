using BufferedStreams

function write_edges(filename ::String, I :: Vector{Int64}, J :: Vector{Int64}, delim = '\t', linesep = '\n')
   write_edges(filename, zip(I, J), delim, linesep)
end

"""
    write_edges(filename, ij, delim, linesep)

Writes the edges given by `ij` to a file with `delim and `linesep`.
"""
function write_edges(filename :: String, ij, delim = '\t', linesep = '\n')
   file = BufferedOutputStream(open(filename, "w"))
   write_edges(file, ij, UInt8(delim), UInt8(linesep))
   close(file)
end

function write_edges(file :: IO, ij, delim :: UInt8, linesep :: UInt8)
   arr = Vector{UInt8}(64) # scratch space for dec!
   for (i, j) in ij
      nchars = dec!(arr, i)
      writeto(file, arr, nchars)
      write(file, delim)
      nchars = dec!(arr, j)
      writeto(file, arr, nchars)
      write(file, linesep)
   end
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
