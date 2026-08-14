"""
    EMRIWaveforms.AmplitudeModels.KerrAmpSpline

Clamped tensor-product cubic B-spline evaluator for the `ZNAmps` Kerr eccentric-equatorial amplitude
grid — the real primitive the `ZNAmps_l10_m10_n55_DS2Outer.h5` file demands (it stores B-spline
COEFFICIENTS over clamped knot vectors, not grid nodes). The evaluator is numerics and lives here; the
HDF5 loader BODY moved to `Data`'s HDF5 weakdep-ext (S4-b, design §7.2/§10: ALL grid reader bodies in
`DataHDF5Ext`) — `read_kerrampgrid`/`kerr_mode_table`/`KerrAmpGrid` are re-exported from `Data` below.

The evaluator mirrors the FEW `interp2D` backend (`few/amplitude/ampinterp2d.py` `AmpInterp2D.__call__`,
which calls the scipy `fpbisp`/`bisplev` machinery over `w_knots`/`u_knots` + flattened coefficients);
`deboor_cubic` is the 1D FITPACK `fpbspl`+`splev` De Boor recursion, `bisplev_clamped` the tensor product.

Surface:
- `deboor_cubic(t, c, x)`         — 1D clamped cubic B-spline point evaluation (De Boor / FITPACK).
- `bisplev_clamped(tw, tu, C, w, u)` — 2D tensor-product cubic B-spline (scipy `bisplev` equivalent).
- `read_kerrampgrid(file; modes)` — re-export of `Data.read_kerrampgrid` (HDF5 body in `DataHDF5Ext`).
- `kerr_mode_table(; lmax, nmax)` — re-export of `Data.kerr_mode_table` (FEW `KerrEccentricEquatorial`
  `(l,m,n)` ↔ grid-index ordering).
"""
module KerrAmpSpline

import Data

export deboor_cubic, bisplev_clamped, read_kerrampgrid, KerrAmpGrid, kerr_mode_table

# =====================================================================================
# 1D clamped cubic B-spline (FITPACK fpbspl + splev). Type-generic over the knot/coeff/point
# element types (BigFloat-inherited): all temporaries promote.
# =====================================================================================

# Locate the knot span: clamped knots ⇒ valid interval index l ∈ [k+1, n−k−1] with t[l] ≤ x ≤ t[l+1].
@inline function _find_interval(t::AbstractVector, k::Int, x)
    n = length(t)
    nk1 = n - k - 1
    l = k + 1
    @inbounds while l < nk1 && x >= t[l + 1]
        l += 1
    end
    return l
end

# FITPACK fpbspl: the k+1 nonzero B-spline basis values B_{l−k},…,B_l at x (in-place into h).
function _fpbspl!(h::AbstractVector{T}, hh::AbstractVector{T}, t::AbstractVector, k::Int, x::T,
                  l::Int) where {T}
    @inbounds begin
        h[1] = one(T)
        for j in 1:k
            for i in 1:j
                hh[i] = h[i]
            end
            h[1] = zero(T)
            for i in 1:j
                li = l + i
                lj = li - j
                if t[li] == t[lj]
                    h[i + 1] = zero(T)
                else
                    f = hh[i] / (t[li] - t[lj])
                    h[i] = h[i] + f * (t[li] - x)
                    h[i + 1] = f * (x - t[lj])
                end
            end
        end
    end
    return h
end

"""
    deboor_cubic(t, c, x; k=3) -> scalar

Evaluate the 1D clamped cubic B-spline with knot vector `t` and coefficients `c`
(`length(c) == length(t) − k − 1`) at `x` via the FITPACK `fpbspl`/`splev` De Boor recursion.
On a clamped knot vector the B-splines partition unity, so constant coefficients reproduce the constant
and Greville-abscissa coefficients reproduce `x` exactly.
"""
function deboor_cubic(t::AbstractVector, c::AbstractVector, x; k::Int = 3)
    T = promote_type(eltype(t), eltype(c), typeof(x))
    h = Vector{T}(undef, k + 1)
    hh = Vector{T}(undef, k)
    l = _find_interval(t, k, x)
    _fpbspl!(h, hh, t, k, T(x), l)
    s = zero(T)
    @inbounds for i in 1:(k + 1)
        s += c[l - k - 1 + i] * h[i]
    end
    return s
end

"""
    bisplev_clamped(tw, tu, C, w, u; k=3) -> scalar

Evaluate the tensor-product cubic B-spline surface with knot vectors `tw` (first/`w` axis) and `tu`
(second/`u` axis) and coefficient matrix `C` of size `(length(tw)−k−1, length(tu)−k−1)` at `(w,u)`:
`Σ_{i,j} C[i,j] B^w_i(w) B^u_j(u)`. Mirrors the FEW `interp2D` / scipy `bisplev` backend.
"""
function bisplev_clamped(tw::AbstractVector, tu::AbstractVector, C::AbstractMatrix, w, u; k::Int = 3)
    T = promote_type(eltype(tw), eltype(tu), eltype(C), typeof(w), typeof(u))
    hw = Vector{T}(undef, k + 1); hhw = Vector{T}(undef, k)
    hu = Vector{T}(undef, k + 1); hhu = Vector{T}(undef, k)
    lw = _find_interval(tw, k, w); _fpbspl!(hw, hhw, tw, k, T(w), lw)
    lu = _find_interval(tu, k, u); _fpbspl!(hu, hhu, tu, k, T(u), lu)
    s = zero(T)
    @inbounds for a in 1:(k + 1)
        iw = lw - k - 1 + a
        ca = zero(T)
        for b in 1:(k + 1)
            ju = lu - k - 1 + b
            ca += C[iw, ju] * hu[b]
        end
        s += ca * hw[a]
    end
    return s
end

# C-order (row-major) flat → Julia matrix: cmat[i,j] = flat[(i−1)*nu + j], i over w (outer), j over u.
# (The ZNAmps coeffs are scipy RectBivariateSpline tck[2] flattened, w outer / u inner.)
function _flat_to_cmat(flat::AbstractVector, nw::Int, nu::Int)
    return permutedims(reshape(flat, nu, nw))
end

# =====================================================================================
# ZNAmps loader + FEW KerrEccentricEquatorial mode ordering — canonical home is the Data spine
# (S4-b, design §7.2/§10: ALL grid reader BODIES live in Data's HDF5 weakdep-ext; EMRIWaveforms'
# direct HDF5 dep is the ext trigger). Thin re-exports keep the historical KerrAmpSpline surface.
# =====================================================================================

"`KerrAmpGrid{T}` — raw ZNAmps coefficient-slab container; canonical home [`Data.KerrAmpGrid`](@ref)."
const KerrAmpGrid = Data.KerrAmpGrid

"""
    kerr_mode_table(; lmax=10, nmax=55) -> (modes, index)

Re-export of [`Data.kerr_mode_table`](@ref): the `(l,m,n)` mode list in ZNAmps grid-axis order
(FEW `KerrEccentricEquatorial`) and the inverse `Dict((l,m,n) => 1-based grid index)`.
"""
kerr_mode_table(; kwargs...) = Data.kerr_mode_table(; kwargs...)

"""
    read_kerrampgrid(file; modes) -> KerrAmpGrid

Re-export of [`Data.read_kerrampgrid`](@ref) (HDF5 body in `DataHDF5Ext`, activated by
EMRIWaveforms' direct HDF5 dep): the `regionA`/`regionB` clamped knot vectors + `CoeffsRegion{A,B}`
B-spline coefficient slabs for the requested `(l,m,n)` `modes`.
"""
read_kerrampgrid(file::AbstractString; kwargs...) = Data.read_kerrampgrid(file; kwargs...)

end # module KerrAmpSpline
