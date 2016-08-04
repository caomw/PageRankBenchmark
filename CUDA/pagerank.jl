using CUDAdrv, CUDAnative

module PageRankCUDA

include("kronGraph500NoPerm.jl")
include("../common/common.jl")
using .Common


#
# Global state
#

type BenchmarkState
    dev::CuDevice
    ctx::CuContext
end

function setup()
    dev = CuDevice(0)
    ctx = CuContext(dev)

    return BenchmarkState(dev, ctx)
end

function teardown(state)
    destroy(state.ctx)
end


#
# Pipeline
#

function kernel0(dir, scl, EdgesPerVertex, state=nothing)
    file = joinpath(dir, "0.tsv")

    ij1, ij2 = kronGraph500NoPerm_shuffle(scl, EdgesPerVertex)
    write_edges(file, ij1, ij2)

    return file
end

end
