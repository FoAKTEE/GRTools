# Wigner-D β=π/2 mode mixing of the orbit-frame puncture modes into the BL harmonic basis. Port of
# Punctures/mode_mixing.jl (h1P-exact-order-eps2.wl, Wardell et al. arXiv:2508.00087), REDESIGNED:
#
#   • the reduced half-π d-matrix `_wignerd_halfpi` and the fixed-angle rotation `_wignerD_fixed`
#     (Mathematica `WignerD[{j,m,mp}, π, π/2, π/2] = (-1)^m i^{mp} d^j_{m,mp}(π/2)`) are CONSUMED from
#     the canonical `Wigner` home (`wignerd_halfpi`/`wignerD_fixed` — AHA: one Wigner-D primitive,
#     BigFloat-inherited). They are re-exported here (identity) for the internal callers and the
#     collapse cross-check; no local d-matrix is re-implemented.
#
# The mixed mode is `h[i][l,m] = Σ_{mp} D^l_{m,mp}(π,π/2,π/2) hP[i][l,mp]` with a per-component
# (sign, mp_start) pattern selecting the even/odd mp combination — the `_MIX_PATTERN` table below.

import Wigner

# ── consumed Wigner-D β=π/2 (re-export identity; collapse the legacy local d-matrix fork) ──────────
const _wignerd_halfpi = Wigner.wignerd_halfpi
const _wignerD_fixed = Wigner.wignerD_fixed

# Per-component mixing pattern (i => (sign, mp_start)): sign=+1 uses the (D[mp]+D[-mp]) combination,
# −1 uses (D[mp]−D[-mp]); mp_start=0 adds a standalone D[m,0] term, then steps by 2. Scalar-type
# i∈{1,3,5,6,7} use + over even mp≥0; i∈{2,4} use − over odd mp; i=8 uses + over odd mp; i∈{9,10}
# use − over even mp≥2 (h1P-...wl mixing block).
const _MIX_PATTERN = Dict(
    1 => (1, 0), 3 => (1, 0), 5 => (1, 0), 6 => (1, 0), 7 => (1, 0),
    2 => (-1, 1), 4 => (-1, 1), 8 => (1, 1),
    9 => (-1, 2), 10 => (-1, 2),
)

"""
    h_mode(i, l, m, r0, M, Δr; mmax=l)

Mixed puncture mode `h[i][l,m,r0,M,Δr]` obtained by Wigner-D rotation of the orbit-frame modes
`hP[i][l,mp]` over `mp` up to `min(l, mmax)` at the fixed Euler angles `(π, π/2, π/2)`. Complex-valued;
inherits the float precision of `r0`,`M`,`Δr`.
"""
function h_mode(i::Integer, l::Integer, m::Integer, r0::Real, M::Real, Δr::Real;
                mmax::Integer = l)
    T = float(promote_type(typeof(r0), typeof(M), typeof(Δr)))
    sign, start = _MIX_PATTERN[i]
    top = min(l, mmax)
    acc = zero(Complex{T})
    if start == 0
        acc += _wignerD_fixed(T, l, m, 0) * hP(i, l, 0, r0, M, Δr)
        start = 2
    end
    for mp in start:2:top
        w = _wignerD_fixed(T, l, m, mp) + sign * _wignerD_fixed(T, l, m, -mp)
        acc += w * hP(i, l, mp, r0, M, Δr)
    end
    return acc
end
