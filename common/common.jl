module Common
    # I/O
    export write_edges, read_edges, uncache
    include("read.jl")
    include("write.jl")
    include("util.jl")
end
