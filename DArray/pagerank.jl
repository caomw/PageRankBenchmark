module Pagerank
using DistributedArrays

include("kronGraph500NoPerm.jl")
include("../io/io.jl")
using .PagerankIO
include("dsparse.jl") # Provides create_adj_matrix

function kernel0(filenames, scl, EdgesPerVertex)
   n = 2^scl # Total number of vertices
   m = EdgesPerVertex * n # Total number of edges

   # Make sure that we distribute the workload over the workers.
   EdgesPerWorker = m ÷ nworkers()
   surplus = m % nworkers()

   @assert length(filenames) == nworkers()

   lastWorker = maximum(workers())
   @sync begin
      for (id, filename) in zip(workers(), filenames)
         nEdges = EdgesPerWorker
         nEdges += ifelse(id == lastWorker, surplus, 0)
         @async remotecall_wait(kronGraph500, id, filename, scl, nEdges)
      end
   end
   return n
end

function kernel1(filenames, path)
   info("Read data")
   @time edges = DArray(dread(filenames)) # DArrayt construction will wait on the futures

   info("Sort edges")
   @time sorted_edges = sort(edges, by = first)
   close(edges)

   info("Write edges")
   filenames = collect(joinpath(path, "chunk_$i.tsv") for i in 1:nworkers())
   @time dwrite(filenames, sorted_edges)
   filenames
end

function kernel2(filenames, N)
   info("Read data and turn it into a sparse matrix")
   @time begin
      rrefs = dread(filenames)
      adj_matrix = create_adj_matrix(rrefs, N)
      rrefs = nothing
   end

   @assert size(adj_matrix) == (N, N)
   info("Pruning and scaling")
   @time begin
      din = sum(adj_matrix, 1)                  # Compute in degree
      adj_matrix[find(din == maximum(din))]=0   # Eliminate the super-node.
      adj_matrix[find(din == 1)]=0              # Eliminate the leaf-node.
      # dout = sum(adj_matrix, 2)                 # Compute out degree
      # is = find(dout)                           # Find vertices with outgoing edges (dout > 0).
      # DoutInvD = zeros(size(adj_matrix, 1))     # Create diagonal weight matrix.
      # DoutInvD[is] = 1./dout[is]
      # scale!(DoutInvD, adj_matrix)              # Apply weight matrix.
   end

   adj_matrix
end

function run(path, scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)
   info("Number of workers: ", nworkers())

   filenames = collect(joinpath(path, "$i.tsv") for i in 1:nworkers())

   info("Executing kernel 0")
   @time N = kernel0(filenames, scl, EdgesPerVertex)

   # Shuffle the filenames so that we minimise cache effect
   # TODO ideally we would like to make sure that no processor reads in
   # its own file.
   shuffle!(filenames)

   info("Executing kernel 1")
   # The filenames changed and for the darray construction we need to retain order...
   @time filenames = kernel1(filenames, path)

   info("Executing kernel 2")
   @time adj_matrix = kernel2(filenames, N)
end

function dread(filenames)
   map(zip(filenames, workers())) do iter
      filename, id = iter
      remotecall(read_edges, id, filename)
   end
end

function dwrite(filenames, edges)
   @sync for (id, filename) in zip(workers(), filenames)
      @async remotecall_wait(id, filename) do filename
         write_edges(filename, localpart(edges))
      end
   end
end

# Two helper functions to make sort work on tuples
Base.typemin{T}(::Type{Tuple{T,T}}) = (typemin(T),typemin(T))
Base.typemax{T}(::Type{Tuple{T,T}}) = (typemax(T),typemax(T))

end
