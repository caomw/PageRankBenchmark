#
# Read/Write Edges from file
#

using BufferedStreams

function writeuv(fname, u, v)
    n = length(u)
    file = BufferedOutputStream(open(fname, "w"))
    for j in 1:n-1
         write(file, dec(u[j]))
         write(file, ' ')
         write(file, dec(v[j]))
         write(file, ' ')
    end
    write(file, dec(u[n]))
    write(file, ' ')
    write(file, dec(v[n]))

    close(file)
end

function StrFileRead(fname)
    file=IOBuffer(Mmap.mmap(open(fname), Vector{UInt8}, (filesize(fname),)))
    ut,vt = readtsv(file)
    close(file)
    return ut,vt
end

function readtsv(file)
    ij1 = Array{Int64}(0)
    ij2 = Array{Int64}(0)
    const SPACE = UInt8(' ')
    while !eof(file)
        push!(ij1, parseuntil(file,SPACE))
        push!(ij2, parseuntil(file,SPACE))
    end
    return ij1, ij2
end


function parseuntil(file, delim::UInt8)
    v = 0
    b = read(file, UInt8)
    const ZERO = UInt8('0')
    while b != delim
        v *= 10
        v += b - ZERO
        if eof(file)
            break
        end
        b = read(file, UInt8)
    end
    v
end
