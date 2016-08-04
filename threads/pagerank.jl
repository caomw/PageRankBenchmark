module PageRankThreads

include("KronGraph500NoPerm_jthread.jl")
include("io.jl")

#include("../common/common.jl")
#using .Common

function setup()
end

function teardown(state)
end

#
# Pipeline
#

function kernel0(dir, scl, EdgesPerVertex, state=nothing)
    Nfile=4
    files = collect(joinpath(dir, "$i.tsv") for i in 1:Nfile)

    n = 2^scl # Total number of vertices
    m = EdgesPerVertex * n # Total number of edges

    EdgesPerFile = EdgesPerVertex./Nfile

    for fname in files
        ut,vt=KronGraph500NoPerm(scl,EdgesPerFile);  # Generate data.
        writeuv(fname, ut, vt)
    end
    
    return dir, files, scl, EdgesPerVertex, Nfile
end

function kernel1(dir, files, scl, EdgesPerVertex, Nfile, state=nothing)
    i=1
    for fname in files
        ut,vt=StrFileRead(fname)
        if i==1
            u=ut
            v=vt
        else
            append!(u, ut)
            append!(v, vt)
        end
        i+=1 
    end
    sortIndex = sortperm(u)                      # Sort starting vertices.
    u = u[sortIndex]                                  # Get starting vertices.
    v = v[sortIndex]                                  # Get ending vertices.

    j=1
    c = size(u,1)/Nfile        # Compute first edge of file.
    for i in 1:Nfile
      jEdgeStart = round(Int, (j-1)*c+1)# Compute first edge of file.
      jEdgeEnd = round(Int, j*c)          # Compute last edge of file.
      uu = view(u,jEdgeStart:jEdgeEnd)                                 # Select start vertices.
      vv = view(v,jEdgeStart:jEdgeEnd)                                 # Select end vertices.
      fname = joinpath(dir, "chunk_"*string(i) * ".tsv")

      writeuv(fname, uu, vv)

      j = j + 1                                                   # Increment file counter.
    end
end

end
