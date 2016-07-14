# execute by running `➜ mpiexec -n $NUM_PROCS julia mpi_pagerank.jl`

import MPI

MPI.Init()

include("utils.jl")
include("kernel0.jl")
include("kernel1.jl")

function main(scale, edges_per_vertex, nfiles)
    @assert isinteger(log(2, edges_per_vertex))
    @assert isinteger(log(2, nfiles))
    @assert edges_per_vertex >= nfiles

    is_root_process = am_i_root_process()

    # set up directories
    path = joinpath(dirname(@__FILE__), "data")
    if is_root_process
        isdir(path) && rm(path, recursive = true)
        mkdir(path)
        for kernel in 0:3
            mkdir(kernel_path(path, kernel))
        end
    end

    # kernel 0
    is_root_process && print("running kernel 0...")
    is_root_process && tic()

    MPI.Barrier(MPI.COMM_WORLD)
    kernel0(path, nfiles, scale, edges_per_vertex)
    MPI.Barrier(MPI.COMM_WORLD)

    is_root_process && println("done (took $(toq()) seconds)")

    # kernel 1
    is_root_process && print("running kernel 1...")
    is_root_process && tic()

    MPI.Barrier(MPI.COMM_WORLD)
    kernel1(path, nfiles, scale, edges_per_vertex)
    MPI.Barrier(MPI.COMM_WORLD)

    is_root_process && println("done (took $(toq()) seconds)")
end

main(4, 8, 4)

MPI.Finalize()
