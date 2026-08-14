# SphereCollocation — Gauss-Legendre collocation surface for the spin-weighted angular sector
# (TeukEvolution.jl Sphere.jl port). swal is a convention adapter over the ₛY seed (no Jacobi
# dep); the sign law (−1)^(l+m+max(|s|−|m|,0)) is gated vs the upstream swal_vals fixture over
# the full (|s|≤4,|m|≤6,l≤12) grid — do not "simplify" it (a shorter form fit from l∈{2,3}
# fails 524/1083 points). {T}-generic; upstream was Float64-pinned.

# ---------------------------------------------------------------------------------------------------
# mode-count truncation buffer
# ---------------------------------------------------------------------------------------------------

"""
    num_l(ny; max_s=2, max_m=6) -> Int

Number of retained `l` angular modes for `ny` Gauss-Legendre collocation points:
`num_l = ny − 2·max_s − 2·max_m` (drops the top `2(max_s+max_m)` collocation-corrupted modes).
`max_s`/`max_m` are the production spin/azimuthal ceilings (upstream `Sphere.jl` hard-codes
`max_s=2, max_m=6`; kept here as defaulted keywords). Requires `ny > 2(max_s+max_m)+4`.
"""
function num_l(ny::Integer; max_s::Integer=2, max_m::Integer=6)::Int
    @assert ny > (2 * max_s + 2 * max_m) + 4
    return ny - 2 * max_s - 2 * max_m
end

# ---------------------------------------------------------------------------------------------------
# pointwise swal value
# ---------------------------------------------------------------------------------------------------

"""
    swal(spin, m_ang, l_ang, y) -> T

Spin-weighted associated-Legendre value in the TeukEvolution collocation convention,
at `y = −cosθ ∈ [−1,1]`:

    swal(s,m,l,y) = √(2π) · (−1)^(l + m + max(|s|−|m|,0)) · Re ₛY_{lm}(acos(−y), 0)

`√(2π)` converts the full-sphere ₛY normalization to the upstream `∫_{−1}^{1}dy` one; the sign
reconciles the Goldberg-vs-Jacobi conventions (gated vs the upstream fixture — see file header).
Returns `0` when `l < max(|m|,|s|)`. Type-generic in `y`.
"""
function swal(spin::Integer, m_ang::Integer, l_ang::Integer, y::T) where {T<:Real}
    F = float(T)
    # upstream n = l − max(|m|,|s|); n<0 ⇒ the harmonic does not exist ⇒ 0.
    l_ang < max(abs(m_ang), abs(spin)) && return zero(F)
    θ = acos(-F(y))                                                    # y = −cosθ chart
    e = l_ang + m_ang + max(abs(spin) - abs(m_ang), 0)                 # (−1)^(l+m+max(|s|−|m|,0))
    sgn = iseven(e) ? 1 : -1
    Y = SpinWeightedSphericalHarmonicY(spin, l_ang, m_ang, θ, zero(θ))
    return sqrt(2 * F(pi)) * sgn * real(Y)
end

# ---------------------------------------------------------------------------------------------------
# Gauss-Legendre collocation grids (y=−cosθ, cosθ, sinθ) + inner product
# ---------------------------------------------------------------------------------------------------

"""
    Y_vals(::Type{T}=Float64, ny) -> Vector{T}

The `ny` Gauss-Legendre collocation points `y = −cosθ` on `[−1,1]` (ascending), from
`NumericsCore.gl_nodes_weights`. Type-generic (BigFloat nodes for `T=BigFloat`).
"""
Y_vals(::Type{T}, ny::Integer) where {T<:AbstractFloat} = gl_nodes_weights(T, ny)[1]
Y_vals(ny::Integer) = Y_vals(Float64, ny)

"""
    cos_vals(::Type{T}=Float64, ny) -> Vector{T}

`cosθ = −y` at the `ny` Gauss-Legendre points.
"""
cos_vals(::Type{T}, ny::Integer) where {T<:AbstractFloat} = -gl_nodes_weights(T, ny)[1]
cos_vals(ny::Integer) = cos_vals(Float64, ny)

"""
    sin_vals(::Type{T}=Float64, ny) -> Vector{T}

`sinθ = √(1−y)·√(1+y)` at the `ny` Gauss-Legendre points.
"""
function sin_vals(::Type{T}, ny::Integer) where {T<:AbstractFloat}
    yv = gl_nodes_weights(T, ny)[1]
    return [sqrt(one(T) - pt) * sqrt(one(T) + pt) for pt in yv]
end
sin_vals(ny::Integer) = sin_vals(Float64, ny)

"""
    inner_product(v1, v2) -> T

Gauss-Legendre discretization of `∫_{−1}^{1} v1·conj(v2) dy`, i.e. `Σ_j w_j v1[j] conj(v2[j])`
with the `length(v1)`-point GL weights. Type/precision follow the vector element type.
"""
function inner_product(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T<:Number}
    n = length(v1)
    n == length(v2) || throw(DimensionMismatch("inner_product: v1,v2 length mismatch"))
    _, wv = gl_nodes_weights(real(float(T)), n)
    s = zero(T)
    @inbounds for j in 1:n
        s += wv[j] * v1[j] * conj(v2[j])
    end
    return s
end

# ---------------------------------------------------------------------------------------------------
# internal: swal collocation table (reused across matrix builders — same values as the upstream
# triple-loop swal recomputation, precomputed once per spin for efficiency)
# ---------------------------------------------------------------------------------------------------

# V[j,k] = swal(spin, m_ang, lvals[k], yv[j]) — ny × length(lvals)
function _swal_table(yv::AbstractVector{T}, spin::Integer, m_ang::Integer,
                     lvals::AbstractVector{<:Integer}) where {T<:Real}
    ny = length(yv)
    V = Matrix{T}(undef, ny, length(lvals))
    @inbounds for (k, l) in enumerate(lvals), j in 1:ny
        V[j, k] = swal(spin, m_ang, l, yv[j])
    end
    return V
end

# ---------------------------------------------------------------------------------------------------
# collocation operator matrices ( v[i]·M[i,j] -> v[j] convention, applied by angular_matrix_mult! )
# ---------------------------------------------------------------------------------------------------

"""
    swal_vals(::Type{T}=Float64, ny, spin, m_ang) -> Matrix{T}

`ny × num_l(ny)` table of `swal(spin, m_ang, l, yⱼ)` over the retained `l = lmin … lmin+num_l−1`
(`lmin = max(|spin|,|m_ang|)`) at the GL points. The spectral-to-collocation basis columns.
"""
function swal_vals(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer) where {T<:AbstractFloat}
    yv = gl_nodes_weights(T, ny)[1]
    nl = num_l(ny)
    lmin = max(abs(spin), abs(m_ang))
    return _swal_table(yv, spin, m_ang, lmin:(lmin + nl - 1))
end
swal_vals(ny::Integer, spin::Integer, m_ang::Integer) = swal_vals(Float64, ny, spin, m_ang)

"""
    swal_laplacian_matrix(::Type{T}=Float64, ny, spin, m_ang) -> Matrix{T}

Collocation matrix of the spin-weighted spherical Laplacian: spectral factor `−(l−s)(l+s+1)`
(the angular eigenvalue) summed over the retained `l`, projected through the GL weights.
"""
function swal_laplacian_matrix(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer) where {T<:AbstractFloat}
    yv, wv = gl_nodes_weights(T, ny)
    nl = num_l(ny)
    lmin = max(abs(spin), abs(m_ang))
    lvals = lmin:(lmin + nl - 1)
    V = _swal_table(yv, spin, m_ang, lvals)
    lap = zeros(T, ny, ny)
    @inbounds for j in 1:ny, i in 1:ny
        acc = zero(T)
        for (k, l) in enumerate(lvals)
            acc -= (l - spin) * (l + spin + 1) * (V[i, k] * V[j, k])
        end
        lap[j, i] = acc * wv[j]
    end
    return lap
end
swal_laplacian_matrix(ny::Integer, spin::Integer, m_ang::Integer) =
    swal_laplacian_matrix(Float64, ny, spin, m_ang)

"""
    swal_raising_matrix(::Type{T}=Float64, ny, spin, m_ang) -> Matrix{T}

Collocation matrix of the spin-raising operator ð: spectral factor `+√((l−s)(l+s+1))` coupling the
spin-`s` basis (row `j`) to the spin-`s+1` basis (row `i`). `lmin = max(|s|,|s+1|,|m|)`.
"""
function swal_raising_matrix(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer) where {T<:AbstractFloat}
    yv, wv = gl_nodes_weights(T, ny)
    nl = num_l(ny)
    lmin = max(abs(spin), abs(spin + 1), abs(m_ang))
    lvals = lmin:(lmin + nl - 1)
    Vs = _swal_table(yv, spin, m_ang, lvals)        # spin s at y_j
    Vsp = _swal_table(yv, spin + 1, m_ang, lvals)   # spin s+1 at y_i
    raise = zeros(T, ny, ny)
    @inbounds for i in 1:ny, j in 1:ny
        acc = zero(T)
        for (k, l) in enumerate(lvals)
            acc += sqrt(T((l - spin) * (l + spin + 1))) * wv[j] * Vs[j, k] * Vsp[i, k]
        end
        raise[j, i] = acc
    end
    return raise
end
swal_raising_matrix(ny::Integer, spin::Integer, m_ang::Integer) =
    swal_raising_matrix(Float64, ny, spin, m_ang)

"""
    swal_lowering_matrix(::Type{T}=Float64, ny, spin, m_ang) -> Matrix{T}

Collocation matrix of the spin-lowering operator ð′: spectral factor `−√((l+s)(l−s+1))` coupling the
spin-`s` basis (row `j`) to the spin-`s−1` basis (row `i`). `lmin = max(|s|,|s−1|,|m|)`.
"""
function swal_lowering_matrix(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer) where {T<:AbstractFloat}
    yv, wv = gl_nodes_weights(T, ny)
    nl = num_l(ny)
    lmin = max(abs(spin), abs(spin - 1), abs(m_ang))
    lvals = lmin:(lmin + nl - 1)
    Vs = _swal_table(yv, spin, m_ang, lvals)        # spin s at y_j
    Vsm = _swal_table(yv, spin - 1, m_ang, lvals)   # spin s−1 at y_i
    lower = zeros(T, ny, ny)
    @inbounds for i in 1:ny, j in 1:ny
        acc = zero(T)
        for (k, l) in enumerate(lvals)
            acc += -sqrt(T((l + spin) * (l - spin + 1))) * wv[j] * Vs[j, k] * Vsm[i, k]
        end
        lower[j, i] = acc
    end
    return lower
end
swal_lowering_matrix(ny::Integer, spin::Integer, m_ang::Integer) =
    swal_lowering_matrix(Float64, ny, spin, m_ang)

"""
    swal_filter_matrix(::Type{T}=Float64, ny, spin, m_ang) -> Matrix{T}

Collocation matrix of the exponential angular low-pass filter, spectral weight
`exp(−30·(l/(nl−1))^10)` (`l` the absolute angular number, `nl = num_l(ny)`), projected through the
GL weights. Applied on the left: `v[i]·M[i,j] -> v[j]`.
"""
function swal_filter_matrix(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer) where {T<:AbstractFloat}
    yv, wv = gl_nodes_weights(T, ny)
    nl = num_l(ny)
    lmin = max(abs(spin), abs(m_ang))
    lvals = lmin:(lmin + nl - 1)
    V = _swal_table(yv, spin, m_ang, lvals)
    filt = zeros(T, ny, ny)
    @inbounds for j in 1:ny, i in 1:ny
        acc = zero(T)
        for (k, l) in enumerate(lvals)
            acc += exp(-30 * (T(l) / T(nl - 1))^10) * (V[i, k] * V[j, k])
        end
        filt[j, i] = acc * wv[j]
    end
    return filt
end
swal_filter_matrix(ny::Integer, spin::Integer, m_ang::Integer) =
    swal_filter_matrix(Float64, ny, spin, m_ang)

"""
    swal_killtop_matrix(::Type{T}=Float64, ny, spin, m_ang, N) -> Matrix{T}

Collocation matrix of the `killtopN` hard truncation: projector onto the lowest `nl−N` retained
`l`-modes (drops the top `N`), projected through the GL weights. Requires `nl > N+1`.
Applied on the left: `v[i]·M[i,j] -> v[j]`.
"""
function swal_killtop_matrix(::Type{T}, ny::Integer, spin::Integer, m_ang::Integer,
                             N::Integer) where {T<:AbstractFloat}
    yv, wv = gl_nodes_weights(T, ny)
    nl = num_l(ny)
    @assert nl > N + 1
    lmin = max(abs(spin), abs(m_ang))
    lvals = lmin:(lmin + nl - 1 - N)
    V = _swal_table(yv, spin, m_ang, lvals)
    filt = zeros(T, ny, ny)
    @inbounds for j in 1:ny, i in 1:ny
        acc = zero(T)
        for k in eachindex(lvals)
            acc += V[i, k] * V[j, k]
        end
        filt[j, i] = acc * wv[j]
    end
    return filt
end
swal_killtop_matrix(ny::Integer, spin::Integer, m_ang::Integer, N::Integer) =
    swal_killtop_matrix(Float64, ny, spin, m_ang, N)

# ---------------------------------------------------------------------------------------------------
# 2D field angular apply
# ---------------------------------------------------------------------------------------------------

"""
    angular_matrix_mult!(f_m, f, mat) -> nothing

Right-multiply a 2D `(nx,ny)` field by an `ny×ny` collocation matrix in the angular direction:
`f_m[i,j] = Σ_k f[i,k]·mat[k,j]` (the `v[i]·M[i,j] -> v[j]` apply convention). In-place into `f_m`;
`f` and `f_m` may not alias.
"""
function angular_matrix_mult!(f_m::AbstractArray{<:Number}, f::AbstractArray{<:Number},
                              mat::AbstractMatrix{<:Number})
    nx, ny = size(f)
    @inbounds for j in 1:ny, i in 1:nx
        acc = zero(eltype(f_m))
        for k in 1:ny
            acc += f[i, k] * mat[k, j]
        end
        f_m[i, j] = acc
    end
    return nothing
end
