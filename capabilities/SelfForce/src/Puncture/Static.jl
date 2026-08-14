# Static h0 first-order field: the l-mode static-sector coefficients h0[i,l,p,r0,M] for a circular
# Schwarzschild geodesic (closed forms in the complete elliptic integrals K,E with parameter
# m = M/(r0−2M), rational in (l,r0,M) with factorial ratios). Port of Punctures/static.jl, REDESIGNED:
#
#   • elliptic K/E COLLAPSED onto the one canonical type-generic `NumericsCore.ellipK`/`ellipE`
#     (Mathematica m=k² convention) — the legacy `_ellipk_agm`/`_ellipe_agm` AGM fork is DROPPED (AHA:
#     one canonical elliptic primitive). `ellipK`/`ellipE` are re-exported here (identity, no copy).
#   • the factorial `_fact(x)=Γ(x+1)` COLLAPSED onto `BigSpecial.fact_big` (AHA with the existing
#     BigSpecial Gamma — one factorial primitive). `fact_big` re-exported.
#   • the worldline Schwarzschild connection CONSUMED from `Spacetime.christoffel`(_dr/_dtheta) (a=0
#     method) — re-exported identity, the puncture worldline self-field never re-owns the connection.
#
# Type-generic {T}: the closed form inherits the float precision of r0/M (pass BigFloat for extended
# precision; the high-p families subtract two large K/E terms and need BigFloat to beat the Float64
# catastrophic cancellation). Opens the Puncture namespace with its first grounded physics primitive;
# the α-integration / β=π/2 mode-mixing / regsum / subtraction / worldline self-field follow as later
# modules.

import NumericsCore
import Spacetime
using BigSpecial: fact_big

# ── consumed elliptic K/E (re-export identity; collapse the legacy AGM fork) ──────────────────────
const ellipK = NumericsCore.ellipK
const ellipE = NumericsCore.ellipE

# the coefficient table calls `_ellipk(m)`/`_ellipe(m)` expecting a REAL result of the working type:
# NumericsCore returns Complex (principal branch); take the real part for the real-m static sector.
_ellipk(m::T) where {T<:AbstractFloat} = real(ellipK(m))
_ellipe(m::T) where {T<:AbstractFloat} = real(ellipE(m))

# Mathematica factorial n! = Γ(n+1) (used as factorial RATIOS, finite only for l above the family's
# threshold). Collapsed onto BigSpecial.fact_big (type-generic Gamma) — one factorial primitive.
_fact(x) = fact_big(x)

# ── consumed worldline Schwarzschild connection (re-export identity; a=0 method) ──────────────────
const christoffel = Spacetime.christoffel
const christoffel_dr = Spacetime.christoffel_dr
const christoffel_dtheta = Spacetime.christoffel_dtheta

# auto-generated closed-form coefficient table (_h0_<i>_<p>(l,r0,M) + the (i,p)→fn _H0_TABLE dict),
# vendored verbatim from Punctures/h0_coefficients.jl; resolves _ellipk/_ellipe/_fact/Pi/sq above.
include("h0_coefficients.jl")

"""
    h0_static(i, l, p, r0, M) -> T

Static-sector multipole coefficient `h0[i,l,p,r0,M]` for a circular Schwarzschild geodesic at areal
radius `r0` about a mass `M`. `i` selects the tensor/parity component (1..10), `p` the radial index
of its `(i,p)` family. The result inherits the float precision of `r0`,`M` (`T = float(promote_type(...))`;
pass BigFloat inputs for extended precision). The closed form is finite only for `l` at or above the
family threshold; an absent `(i,p)` throws `ArgumentError`.
"""
function h0_static(i::Integer, l::Integer, p::Integer, r0::Real, M::Real)
    f = get(_H0_TABLE, (Int(i), Int(p)), nothing)
    f === nothing &&
        throw(ArgumentError("h0_static: no static coefficient for (i=$i, p=$p)"))
    T = float(promote_type(typeof(r0), typeof(M)))
    return T(f(T(l), T(r0), T(M)))
end
