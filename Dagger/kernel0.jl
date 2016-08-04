using Dagger

@everywhere begin
    """
    Given the scale, generate a single vertex pair representing an edge

    This is the Graph 500 algorithm for generating a random edge
    """
    function rand_edge_gen(scl)
        a = 0.57; b = 0.19; c = 0.19; d = 1-(a+b+c)
        ab = a+b
        c_norm = c/(1-(a+b))
        a_norm = a/ab
      
        function (_)
            f = 1
            t = 1
            for ib=1:scl
                k = 1 << (ib-1)
                if rand() > ab
                    f += k
                    if rand() > c_norm
                        t += k
                    end
                elseif rand() > a_norm
                    t += k
                end
            end
            f,t
        end
    end
end

"""
Generate avg_connections*(2^scl) edges serially
"""
function generate_serial(scl, avg_connections)
  N = 2^scl # no of vertex
  M = round(Int, avg_connections .* N) # no of edges
  map(rand_edge_gen(scl), 1:M)
end

"""
Generate avg_connections*(2^scl) edges in parallel
"""
function generate_par(scl, avg_connections, nparts)
  N = 2.^scl # no of vertex
  M = round(Int, avg_connections .* N) # no of edges
  p = BlockPartition(ceil(Int, M/nparts))
  map(rand_edge_gen(scl), Distribute(p, 1:M))
end

@everywhere function write_file(io, data)
    for (x,y) in data
        show(io, x)
        write(io, '\t')
        show(io, y)
        write(io, '\n')
    end
end

function write_files(data, nparts, path)
    files = [string("data-$i") for i in 1:nparts]
    p = BlockPartition(1)
    dist_files = Distribute(p, files)
    mappart(data, dist_files) do d, f
        open(joinpath(path, f[1] * ".tsv"), "w") do io
            write_file(io, d)
        end
        d
    end
end

function kernel0(scl, avg_connections, nparts, ctx=Context(), path="data")
    X = generate_par(scl, avg_connections, nparts)
    compute(ctx, write_files(X, nparts, path))
end
