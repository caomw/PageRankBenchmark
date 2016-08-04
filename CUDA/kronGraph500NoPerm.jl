using CUDAdrv, CUDAnative


#
# Utility functions
#

@target ptx function xorshift64(key::UInt64)
    key = (~key) + (key << 21)
    key = key $ (key >> 24)
    key = (key + (key << 3)) + (key << 8)
    key = key $ (key >> 14)
    key = (key + (key << 2)) + (key << 4)
    key = key $ (key >> 28)
    key = key + (key << 31)
    return key
end

"Return a 64-bit floating point number between 0 and 1"
@target ptx function gpurand(seed)
    input = reinterpret(UInt64, Float64(seed))
    val = xorshift64(input)::UInt64
    output = reinterpret(Float64, val / typemax(UInt64))
    return output
end


#
# Basic version with shared memory parallel reduction
#

function kronGraph500NoPerm(scl, EdgesPerVertex)
    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij1_d, ij2_d = CuArray(Int, m), CuArray(Int, m) # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    m_x = min(65535, m)
    m_y = ceil(Int, m/m_x)
    @cuda ((m_x, m_y), scl, scl*4*sizeof(Int)) kronGraph500NoPerm_kernel(m, ab, a_norm, c_norm, ij1_d, ij2_d)

    ij1 = Array(ij1_d)
    free(ij1_d)
    ij2 = Array(ij2_d)
    free(ij2_d)

    return ij1, ij2
end

@target ptx function kronGraph500NoPerm_kernel{T}(m, ab, a_norm, c_norm, ij1::CuDeviceArray{T}, ij2::CuDeviceArray{T})
    i = (blockIdx().y-1) * gridDim().x + blockIdx().x

    ib = threadIdx().x
    scl = blockDim().x

    if i <= m
        # get references to shmem
        temp = @cuDynamicSharedMem(T, 4*scl)
        ij1_buf = view(temp, 1:2*scl)
        ij2_buf = view(temp, 2*scl+1:4*scl)

        seed64 = Int64(ib) << 32 + i
        a = gpurand(seed64)
        b = gpurand(a)

        # calculate per-thread values
        k = 1 << (ib - 1)
        ii_bit  = a > ab
        jj_bit  = b > ifelse(ii_bit, c_norm, a_norm)
        ij1_buf[ib] = k * ii_bit
        ij2_buf[ib] = k * jj_bit
        sync_threads()

        # parallel reduction
        pin, pout = 1, 0
        offset = 1
        while offset < scl
            pout = 1 - pout
            pin = 1 - pin
            if ib > offset
                ij1_buf[pout * scl + ib] = ij1_buf[pin * scl + ib] + ij1_buf[pin * scl + ib - offset]
                ij2_buf[pout * scl + ib] = ij2_buf[pin * scl + ib] + ij2_buf[pin * scl + ib - offset]
            else
                ij1_buf[pout * scl + ib] = ij1_buf[pin * scl + ib]
                ij2_buf[pout * scl + ib] = ij2_buf[pin * scl + ib]
            end
            sync_threads()
            offset = offset * 2
        end
        ij1_buf[pin * scl + ib] = ij1_buf[pout * scl + ib]
        ij2_buf[pin * scl + ib] = ij2_buf[pout * scl + ib]
        sync_threads()

        # write back (also includes "initialization" with one)
        if ib == scl
            ij1[i] = one(T) + ij1_buf[ib]
            ij2[i] = one(T) + ij2_buf[ib]
        end
    end

    return nothing
end


#
# Version using intra-warp shuffle and intra-block shmem
#

function kronGraph500NoPerm_shuffle(scl, EdgesPerVertex)
    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij1_d, ij2_d = CuArray(Int, m), CuArray(Int, m) # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    m_x = min(65535, m)
    m_y = ceil(Int, m/m_x)
    @cuda ((m_x, m_y), nearest_warpsize(scl)) kronGraph500NoPerm_shuffle_kernel(m, scl, ab, a_norm, c_norm, ij1_d, ij2_d)

    ij1 = Array(ij1_d)
    free(ij1_d)
    ij2 = Array(ij2_d)
    free(ij2_d)

    return ij1, ij2
end

@target ptx function reduce_warp{T,F<:Function}(val::T, op::F)
    offset = warpsize ÷ 2
    while offset > 0
        val = op(val, shfl_down(val, offset))
        offset ÷= 2
    end
    return val
end

@target ptx function reduce_block{T,F<:Function}(val::T, op::F)
    shared = @cuStaticSharedMem(T, 32)

    wid, lane = fldmod1(threadIdx().x, warpsize)

    val = reduce_warp(val, op)

    if lane == 1
        shared[wid] = val
    end

    sync_threads()

    # read from shared memory only if that warp existed
    val = (threadIdx().x <= fld(blockDim().x, warpsize)) ? shared[lane] : zero(T)

    if wid == 1
        # final reduce within first warp
        val = reduce_warp(val, op)
    end

    return val
end

@target ptx function kronGraph500NoPerm_shuffle_kernel{T}(m, scl, ab, a_norm, c_norm, ij1::CuDeviceArray{T}, ij2::CuDeviceArray{T})
    i = (blockIdx().y-1) * gridDim().x + blockIdx().x

    ib = threadIdx().x

    if i <= m && ib <= scl
        seed64 = Int64(ib) << 32 + i
        a = gpurand(seed64)
        b = gpurand(a)

        # calculate per-thread values
        k = 1 << (ib - 1)
        ii_bit  = a > ab
        jj_bit  = b > ifelse(ii_bit, c_norm, a_norm)
        ij1_val = k * ii_bit
        ij2_val = k * jj_bit

        # parallel reduction
        ij1_val = reduce_block(ij1_val, +)
        ij2_val = reduce_block(ij2_val, +)

        # write back (also includes "initialization" with one)
        if ib == 1
            ij1[i] = one(T) + ij1_val
            ij2[i] = one(T) + ij2_val
        end
    end

    return nothing
end
