"""
Creates a sparse adjacency matrix from an list of edges stored as Array of Tuple{Int, Int}
on remote nodes. It takes a list of futures and the maximum number of vertices.

# Assumptions
- Edge arrays are sorted after the first elemnt of the tuple
- The fvectors is sorted such that last(fvectors[i]) < first(fvectors[i+1])
- The edges vectors are not overlapping
"""
function create_adj_matrix(fvectors, N)
   # List of participating workers
   pworkers = map(f -> f.where, fvectors)
   # The first worker will need to make sure that it starts at 1
   firstWorker = first(pworkers)
   # The last worker will include up to N elements
   lastWorker = last(pworkers)

   # Turn the Vector of Tuple{Int, Int} into the arrays I, J, V
   # and collect the minimum and maximum
   rrefs = map(fvectors) do fvec
      remotecall(fvec.where) do
         edges = fetch(fvec)
         # maybe we should do a sort! here to ensure consistency

         # SparseMatrixCSC only accepts I, J, V for construction
         I = Vector{Int64}(length(edges))
         J = Vector{Int64}(length(edges))
         V = Vector{Int64}(length(edges))

         min_i = ifelse(myid() == firstWorker, 1, first(edges)[1])
         max_i = last(edges)[1]

         for ind in eachindex(edges, I, J, V)
            i, j = edges[ind]
            i = i - (min_i - 1) # localparts of sparse matrix need to start at 1
            I[ind] = i
            J[ind] = j
            V[ind] = 1
         end
         (I, J, V, max_i, min_i)
      end
   end

   # Collect the maxima
   max_is = map(rrefs) do rref
      remotecall_fetch(r -> fetch(r)[4], rref.where, rref)
   end

   # The last worker needs to extend the number of elements it knows about
   surplus = N - maximum(max_is)

   # Collect the minima
   min_is = map(rrefs) do rref
      remotecall_fetch(r -> fetch(r)[5], rref.where, rref)
   end

   # Consistency check
   for i in 1:(length(rrefs)-1)
      @assert max_is[i] < min_is[i+1]
   end

   # Construct sparse array
   lparts = map(rrefs) do rref
      remotecall(rref.where, rref) do ref
         (I, J, V, max_i, min_i) = fetch(ref)
         max_i = ifelse(myid() == lastWorker, max_i + surplus, max_i) - (min_i - 1)
         sparse(I, J, V, max_i, N)
      end
   end
   return DArray(reshape(lparts, (length(lparts), 1)))
end
