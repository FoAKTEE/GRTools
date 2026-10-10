# Iterative radix-2 decimation in time. No dependencies beyond Base, no global
# caches, and one call-local twiddle table reused across every batch of an axis.

"""
    fft_radix2!(A::AbstractArray{Complex{T}}; dims=1:ndims(A), inverse=false) -> A

In-place radix-2 Fourier transform along `dims` (an integer or iterable of distinct
dimensions). For a length-`N` axis, the forward convention is
`Y[k+1] = sum(A[j+1] * exp(-2π*im*j*k/N) for j in 0:N-1)`.
`inverse=true` reverses the sign and divides by `N` on each transformed axis.
Bins are in ordinary, unshifted order; `dims=()` is an identity operation.

Every transformed length must be a positive power of two (including one), and
indices must be one-based. Invalid dimensions/lengths throw `ArgumentError`
before mutation. Untransformed batch dimensions may have arbitrary lengths.
Writable strided views are supported; elements of `A` must not alias each other.

`T <: AbstractFloat` is preserved. BigFloat arithmetic inherits the ambient
precision, as elsewhere in NumericsCore; this function never changes precision.
Set that precision before starting concurrent work and keep it fixed until all
such work finishes. Calls on disjoint storage are thread-safe: all twiddles are
call-local, with no shared plans, caches, or mutable scratch buffers.

Work is `O(length(A) * sum(log2(size(A,d)) for d in dims))`. Each axis allocates
one half-length `Complex{T}` twiddle table, reused across its batches, with no
full-size data copy. BigFloat arithmetic additionally allocates its scalar
results. Twiddles are evaluated directly with `sincospi`, avoiding recurrence
drift. No arbitrary-length fallback or packed real transform is provided.
"""
function fft_radix2!(A::AbstractArray{Complex{T}};
                     dims=1:ndims(A), inverse::Bool=false) where {T<:AbstractFloat}
    Base.require_one_based_indexing(A)
    transform_dims = dims isa Integer ? (dims,) : Tuple(dims)
    for (i, d) in enumerate(transform_dims)
        d isa Integer && 1 <= d <= ndims(A) ||
            throw(ArgumentError("FFT dimensions must be integers in 1:$(ndims(A)); got $d"))
        any(j -> transform_dims[j] == d, 1:i-1) &&
            throw(ArgumentError("FFT dimensions must be distinct; repeated dimension $d"))
        ispow2(size(A, d)) ||
            throw(ArgumentError("FFT length must be a positive power of two; dimension $d has length $(size(A, d))"))
    end
    isempty(A) && return A
    for d in transform_dims
        _fft_radix2_axis!(A, Int(d), inverse)
    end
    return A
end

"""
    fft_radix2(A::AbstractArray{<:Number}; dims=1:ndims(A), inverse=false)

Return a complex floating-point copy transformed by [`fft_radix2!`](@ref).
Real input returns the full complex spectrum; integers use Float64 by default,
and floating input preserves its real component type (including BigFloat).
`A` is unchanged. The sign, normalization, dimension and precision conventions
are those of `fft_radix2!`.
"""
function fft_radix2(A::AbstractArray{<:Number}; dims=1:ndims(A), inverse::Bool=false)
    return fft_radix2!(complex.(float.(A)); dims, inverse)
end

function _fft_radix2_axis!(A::AbstractArray{Complex{T},D}, d::Int,
                          inverse::Bool) where {T,D}
    n = size(A, d)
    n == 1 && return A
    twiddles = Vector{Complex{T}}(undef, n ÷ 2)
    direction = inverse ? 1 : -1
    @inbounds for k in 0:n÷2-1
        s, c = sincospi(T(direction * 2k) / T(n))
        twiddles[k+1] = complex(c, s)
    end
    step = CartesianIndex(ntuple(i -> Int(i == d), D))
    origins = CartesianIndices(ntuple(i -> i == d ? Base.OneTo(1) : axes(A, i), D))
    for origin in origins
        _fft_radix2_line!(A, origin, step, twiddles, n, inverse)
    end
    return A
end

function _fft_radix2_line!(A::AbstractArray{Complex{T}}, origin, step,
                          twiddles, n::Int, inverse::Bool) where {T}
    # Bit-reversal permutation, using zero-based offsets within this axis.
    j = 0
    @inbounds for i in 0:n-1
        if i < j
            a, b = origin + i * step, origin + j * step
            A[a], A[b] = A[b], A[a]
        end
        bit = n >> 1
        while bit != 0 && (j & bit) != 0
            j ⊻= bit
            bit >>= 1
        end
        j ⊻= bit
    end
    width = 2
    @inbounds while width <= n
        half = width ÷ 2
        twiddle_step = n ÷ width
        for first in 0:width:n-1, k in 0:half-1
            a = origin + (first + k) * step
            b = a + half * step
            u = A[a]
            v = twiddles[k * twiddle_step + 1] * A[b]
            A[a], A[b] = u + v, u - v
        end
        width == n && break # avoids overflow after the largest possible axis
        width *= 2
    end
    if inverse
        scale = inv(T(n))
        @inbounds for i in 0:n-1
            index = origin + i * step
            A[index] *= scale
        end
    end
    return A
end
