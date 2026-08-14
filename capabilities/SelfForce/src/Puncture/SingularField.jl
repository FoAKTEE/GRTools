# Singular-field l-modes hS[i,l,m] sampled on the worldtube Δr-grid + the worldline left/right split.
# Port of Punctures/singular_field.jl (h1S-exact-order-eps2.wl L415-475, L550-555, Wardell et al.
# arXiv:2508.00087).
#
# On the grid the singular field IS the mode-mixed puncture mode `h[i][l,m]` evaluated at each grid
# displacement Δr (upstream `hS[i,l,m] = h[i][l,m,mmax,r0,M,ΔrgridP]`); two source modes carry no
# puncture content and are set identically zero (`hS[2,0,0]`, `hS[9,1,0]`). The grid carries the
# worldline at Δr=0, which is a NaN of the puncture integrand — its one-sided limits are recovered by
# `split_worldline` from off-worldline samples.

"""
    SINGULAR_ZERO_MODES

Modes `(i,l,m)` whose singular field is identically zero in the source (`hS[2,0,0]` even t-r static
monopole, `hS[9,1,0]` odd static; h1S-...wl L460-475) — no puncture content.
"""
const SINGULAR_ZERO_MODES = ((2, 0, 0), (9, 1, 0))

"""
    hS_mode(i, l, m, r0, M, Δrgrid; mmax=10)

Singular-field mode `hS[i,l,m]` sampled on `Δrgrid` (radial displacements Δr = r − r0): a complex
vector equal to the mode-mixed puncture mode `h[i][l,m]` at each Δr, except the identically-zero
`SINGULAR_ZERO_MODES` which short-circuit to exact zeros. `mmax` truncates the mp mixing sum.

The puncture integrand is singular at Δr=0, so a grid point exactly on the worldline yields NaN;
sample off the worldline and recover the one-sided limits via `split_worldline`.
"""
function hS_mode(i::Integer, l::Integer, m::Integer, r0::Real, M::Real,
                 Δrgrid::AbstractVector; mmax::Integer = 10)
    T = float(promote_type(typeof(r0), typeof(M), eltype(Δrgrid)))
    if (Int(i), Int(l), Int(m)) in SINGULAR_ZERO_MODES
        return zeros(Complex{T}, length(Δrgrid))
    end
    return Complex{T}[h_mode(i, l, m, r0, M, Δr; mmax = mmax) for Δr in Δrgrid]
end

"""
    worldline_index(Δrgrid; atol=0)

Index of the worldline (Δr = 0) in `Δrgrid` (upstream `r0iP = Position[ΔrgridP, 0]`). With the
default `atol=0` an exact zero is required, as the grid carries the worldline point explicitly.
"""
function worldline_index(Δrgrid::AbstractVector; atol::Real = 0)
    for (k, Δr) in enumerate(Δrgrid)
        abs(Δr) <= atol && return k
    end
    throw(ArgumentError("worldline_index: no Δr=0 point in grid"))
end

"""
    split_worldline(values, r0i)

Split a per-grid-point vector into `(left, right)` at the worldline index `r0i`, following the source
convention `[[;;r0iP]]` / `[[r0iP+1;;]]`: the left limit includes the worldline point, the right limit
excludes it.
"""
function split_worldline(values::AbstractVector, r0i::Integer)
    return (values[1:r0i], values[r0i+1:end])
end
