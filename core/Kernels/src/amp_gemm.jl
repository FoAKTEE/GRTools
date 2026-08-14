# amp_gemm.jl — ROMAN dense neural-network forward pass (GEMM + leaky-ReLU).
#
# Canonical home for `roman_forward`/`leaky_relu` (the ROMAN amplitude surrogate's dense NN). BLAS-guarded
# (NumericsCore.with_blas_threads) inside any threaded sweep. Grounded by a small hand-computed reference
# GEMM + leaky_relu identity. Built tests-first in a subsequent round.
#
# Consolidates: KerrBackwards·romanamplitude·{roman_amplitudes,leaky_relu}; FEW·romannet·leaky_relu.
#
# Kernels owns ONLY the pure GEMM given (weights, biases, reduced_basis, transform_factor_inv, x). HDF5
# weight loading + ModeIndexTable/special_index_map conj-partner selection stay in
# EMRIWaveforms.AmplitudeModels.

export leaky_relu, roman_forward

"""
    leaky_relu(x)

ROMAN-network activation: identity for `x ≥ 0`, slope `0.2` otherwise (the fixed `0.2` leak of
`matmul.cu::LeakyReLU`). Type-generic in `x`.
"""
leaky_relu(x::Real) = x >= 0 ? x : 0.2 * x

"""
    roman_forward(weights, biases, reduced_basis, transform_factor_inv, x; threaded=false) -> Vector

Pure dense-NN forward pass of the ROMAN Teukolsky-amplitude surrogate (the de-dup of KerrBackwards
`roman_amplitudes` and FEW `romannet`, HDF5/table stripped). Each layer applies `v → Wᵢ·v + bᵢ` with
`leaky_relu` on all but the last layer; the last layer's `2·bi` outputs (`bi = length/2`) fold into `bi`
complex coefficients `nnⱼ = (vⱼ + i v_{bi+j})·transform_factor_inv`, projected through `reduced_basis`
(`num_teuk_modes × bi`) to the mode amplitudes.

`weights[i]` is the `(n_out, n_in)` matrix (the HDF5 `Wᵀ`); `x` the input feature vector. One work-item
per output mode-row + a serial dot over `bi` ⇒ the serial and bounded-threaded projections are
bit-identical (design §8 bit-equality gate); `threaded=true` chunks the rows ≤20, BLAS-guarded.
"""
function roman_forward(weights::AbstractVector, biases::AbstractVector,
                       reduced_basis::AbstractMatrix, transform_factor_inv,
                       x::AbstractVector; threaded::Bool = false)
    nlayers = length(weights)
    nlayers == length(biases) ||
        throw(ArgumentError("weights and biases must have equal length"))
    v = x
    @inbounds for i in 1:nlayers
        v = weights[i] * v .+ biases[i]
        i < nlayers && (v = leaky_relu.(v))
    end
    iseven(length(v)) ||
        throw(ArgumentError("last layer output length $(length(v)) must be even (2·break_index)"))
    bi = length(v) ÷ 2
    size(reduced_basis, 2) == bi ||
        throw(ArgumentError("reduced_basis must have $bi columns (break_index), got $(size(reduced_basis,2))"))

    nn = [complex(v[j], v[bi + j]) * transform_factor_inv for j in 1:bi]
    nrows = size(reduced_basis, 1)
    Tout = typeof(zero(eltype(reduced_basis)) * zero(eltype(nn)))
    out = Vector{Tout}(undef, nrows)
    with_kernel_threads(nrows; threaded = threaded) do r
        acc = zero(Tout)
        @inbounds for c in 1:bi
            acc += reduced_basis[r, c] * nn[c]
        end
        @inbounds out[r] = acc
    end
    return out
end
