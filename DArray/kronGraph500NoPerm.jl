function kronGraph500(filename, scl, nEdges)
    edges = kronGraph500NoPerm(scl, nEdges)
    write_edges(filename, edges)
end

# change loop order to reduce loads and stores
function kronGraph500NoPerm(scl, nEdges)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij = Vector{Tuple{Int64, Int64}}(nEdges)   # Initialize index arrays.
    ab = a + b                                 # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    randbuf = Vector{Float64}(2scl)
    for i = 1:nEdges
        rand!(randbuf)
        ij1 = one(Int64)
        ij2 = one(Int64)
        @inbounds for ib = 1:scl                   # Loop over each scale.
            sc = 1 << (ib - 2)
            ii_bit  = randbuf[ib] > ab
            jj_bit  = randbuf[ib+scl] > ifelse(ii_bit, c_norm, a_norm)
            ij1 += sc * ii_bit
            ij2 += sc * jj_bit
        end
        ij[i] = (ij1, ij2)
    end

    return ij
end

# FASTEST VERSION
# With threading
function kronGraph500NoPermWithThreads(scl, EdgesPerVertex)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij = Vector{Tuple{Int64, Int64}}(m)   # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    @assert m % Threads.nthreads() == 0
    Base.Threads.@threads for l = 1:Threads.nthreads()
        # ccall(:jl_, Void, (Any,), Threads.threadid())
        r = Base.MersenneTwister(Threads.threadid())
        randbuf = Vector{Float64}(2scl)
        range = 1:div(m, Threads.nthreads())
        lrange = length(range)*(Threads.threadid() - 1) + range
        for i = lrange
            # ccall(:jl_, Void, (Any,), i)
            # tl_randbuf = randbuf[Threads.threadid()]
            rand!(r, randbuf)
            ij1 = one(Int)
            ij2 = one(Int)
            @inbounds for ib = 1:scl                   # Loop over each scale.
                sc = 1 << (ib - 2)
                ii_bit  = randbuf[ib] > ab
                jj_bit  = randbuf[ib+scl] > ifelse(ii_bit, c_norm, a_norm)
                # ii_bit  = rand() > ab
                # jj_bit  = rand() > ifelse(ii_bit, c_norm, a_norm)
                ij1 += sc * ii_bit
                ij2 += sc * jj_bit
            end
            ij[i] = (ij1, ij2)
        end
    end

    return ij
end
