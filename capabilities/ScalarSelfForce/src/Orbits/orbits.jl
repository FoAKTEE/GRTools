# Orbits — CircularOrbit{T}: the orbit-level adapter+cache over Spacetime for the m-mode
# elliptic scalar-self-force pipeline (ESSF node `orbits-adapters`; arXiv:2606.04998
# Sec. II.A, main.tex:176-209 + puncture scales main.tex:445-455). PURE ADAPTER: every
# physical quantity comes from an exported Spacetime function (constants_of_motion,
# kerr_circular_ut_uphi, horizons, kerr_isco, delta_kerr) — no geodesic physics is
# reimplemented here; the only formulas OWNED by this file are the two puncture-scale
# constants (ϱ² prefactors, z_c²), transcribed from main.tex:452-453 and gated against
# those literals in test/test_orbits.jl.

"""
    CircularOrbit(a, r0) -> CircularOrbit{T}

All orbit-level constants of a circular equatorial Kerr geodesic at Boyer–Lindquist radius
`r0` (G = c = 1, M = 1), assembled once from Spacetime's exported closed forms and cached
for the downstream m-mode pipeline (arXiv:2606.04998 Sec. II.A). `T = float(promote_type(...))`
of the arguments, so BigFloat precision is inherited (never down-cast).

Sign convention (Spacetime / main.tex:196-198): the motion is ALWAYS in +φ (`uphi > 0`,
`Omega > 0`); `a > 0` is prograde, `a < 0` retrograde (`prograde == (a >= 0)`).

Fields:
  - `M`   — background mass, `one(T)` (fixed package-wide).
  - `a`   — Kerr spin, `|a| < 1`.
  - `r0`  — orbital Boyer–Lindquist radius (exterior, outside the photon orbit; radii
            inside the ISCO are valid UNSTABLE circular geodesics).
  - `prograde` — `a >= 0` (sign convention above).
  - `E`, `Lz` — specific energy / axial angular momentum, `constants_of_motion(a, r0, 0, 1)`
            (equatorial `x = +1`; retrograde is encoded by `a < 0`, main.tex:176-190).
  - `Omega` — orbital frequency `uphi/ut = 1/(a + √(r0³/M))` (main.tex:196). Deliberately
            NOT `kerr_frequencies`: its circ-eq Mino path DomainErrors at exactly
            `r0 = kerr_isco(a)` (roundoff; ESSF error ledger 2026-07-28), while this
            closed form is regular through the ISCO.
  - `ut`, `uphi` — four-velocity components, `kerr_circular_ut_uphi(a, r0)` (main.tex:192-196).
  - `rp`, `rm` — horizons r± = M ± √(M² − a²), `horizons(a)`.
  - `r_isco` — `kerr_isco(a)` (Bardeen closed form; main.tex:207-209 condition).
  - `rho2_dr2`, `rho2_dtheta2` — puncture-scale ϱ² prefactors (main.tex:452):
            `ϱ² = rho2_dr2·Δr² + rho2_dtheta2·Δθ²` with `rho2_dr2 = r0²/Δ(r0)`
            (note r0² − 2M r0 + a² ≡ Δ(a, r0)) and `rho2_dtheta2 = r0²`.
  - `zc2` — puncture azimuthal scale (main.tex:453): `z_c² = 4(Lz² + r0² + a² + 2M a²/r0)`.

Invalid inputs surface as the natural `DomainError`s of the underlying closed forms
(`|a| > 1` in `horizons`; `r0` at/inside the photon orbit in `kerr_circular_ut_uphi`).
"""
struct CircularOrbit{T<:AbstractFloat}
    M::T
    a::T
    r0::T
    prograde::Bool
    E::T
    Lz::T
    Omega::T
    ut::T
    uphi::T
    rp::T
    rm::T
    r_isco::T
    rho2_dr2::T
    rho2_dtheta2::T
    zc2::T
end

function CircularOrbit(a::Real, r0::Real)
    T = float(promote_type(typeof(a), typeof(r0)))
    aT, r0T = T(a), T(r0)
    com = constants_of_motion(aT, r0T, zero(T), one(T))    # equatorial x = +1
    ut, uphi = kerr_circular_ut_uphi(aT, r0T)
    rp, rm = horizons(aT)
    # puncture scales, main.tex:452-453 (M = 1); r0² − 2M r0 + a² ≡ Δ(a, r0)
    rho2_dr2 = r0T^2 / delta_kerr(aT, r0T)
    zc2 = 4 * (com.Lz^2 + r0T^2 + aT^2 + 2 * aT^2 / r0T)
    CircularOrbit{T}(one(T), aT, r0T, aT >= 0, com.E, com.Lz, uphi / ut, ut, uphi,
                     rp, rm, kerr_isco(aT), rho2_dr2, r0T^2, zc2)
end
