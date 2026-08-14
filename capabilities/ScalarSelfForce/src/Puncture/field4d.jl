# Puncture — 4D Detweiler-Whiting puncture field Φ^P (ESSF node `puncture-field-4d`;
# arXiv:2606.04998 Sec. II.C, eq:puncture + eq:s2, main.tex:437-457). The A_ijk coefficient
# closed forms are a MECHANICAL transcription of kerr-circular.c effsource_set_particle
# (generated file field4d_coeffs.jl — see its header and generator/transcribe_field4d.py);
# the ϱ²/z_c² scales are NOT recomputed here: they live on CircularOrbit (rho2_dr2,
# rho2_dtheta2, zc2 — the same main.tex:452-453 literals the C bakes into alpha20/alpha02/
# beta, kerr-circular.c:844-847; single knowledge home, gated in test_orbits.jl).

include("field4d_coeffs.jl")   # generated: _puncture_Aijk (50 A_ijk closed forms)

"""
    Puncture(orbit::CircularOrbit{T}) -> Puncture{T}

The fourth-order 4D Detweiler-Whiting puncture field for a circular equatorial Kerr
geodesic (arXiv:2606.04998 eq:puncture, main.tex:440-457):

    Φ^P(Δr, Δθ, Δφ) = R⁻⁷ Σ_{ijk} A_ijk Δr^i Δθ^j Q^k,     6 ≤ i+j+k ≤ 9,

with `Q = sin(Δφ/2)`, `R² = ϱ² + z_c² Q²` (eq:s2), `ϱ² = r_p²/(r_p²−2M r_p+a²)·Δr² +
r_p²·Δθ²` and `z_c² = 4(L² + r_p² + a² + 2M a²/r_p)` — the scales are consumed directly
from `orbit.rho2_dr2/rho2_dtheta2/zc2`. On circular equatorial orbits A_ijk = 0 whenever
`j` or `k` is odd (main.tex:457), so exactly the 50 even-j, even-k coefficients are stored
(fields named after the C statics A006…A900 of kerr-circular.c:18, evaluated by the
transcribed `_puncture_Aijk` at construction). Type-generic: a `CircularOrbit{BigFloat}`
yields full-precision BigFloat coefficients.

Evaluate with [`puncture_phi`](@ref). The m-mode integrals of this field (₃P/2F1 route)
and its analytic derivatives are separate nodes (`mmode-puncture-integrals`,
`puncture-derivatives`).
"""
struct Puncture{T<:AbstractFloat}
    orbit::CircularOrbit{T}
    # the 50 A_ijk, in the C static order (kerr-circular.c:18) == sorted name order ==
    # the field order of the _puncture_Aijk NamedTuple (constructor splats it).
    A006::T; A008::T; A024::T; A026::T; A042::T; A044::T; A060::T; A062::T; A080::T
    A106::T; A108::T; A124::T; A126::T; A142::T; A144::T; A160::T; A162::T; A180::T
    A204::T; A206::T; A222::T; A224::T; A240::T; A242::T; A260::T
    A304::T; A306::T; A322::T; A324::T; A340::T; A342::T; A360::T
    A402::T; A404::T; A420::T; A422::T; A440::T
    A502::T; A504::T; A520::T; A522::T; A540::T
    A600::T; A602::T; A620::T
    A700::T; A702::T; A720::T
    A800::T
    A900::T
end

function Puncture(orbit::CircularOrbit{T}) where {T}
    # C call chain: effsource_init(M, a); effsource_set_particle(xp, E, L, ur_p=0)
    A = _puncture_Aijk(orbit.M, orbit.a, orbit.r0, orbit.E, orbit.Lz)
    return Puncture{T}(orbit, A...)
end

"""
    puncture_phi(p::Puncture, dr, dtheta, dphi) -> Φ^P

The 4D puncture field at coordinate offsets `dr = r − r_p`, `dtheta = θ − π/2`,
`dphi = φ − φ_p` from the particle (eq:puncture; transcription of `effsource_PhiS`,
kerr-circular.c:75-120). Real-valued; diverges as 1/λ on approach to the worldline and is
even in `dtheta` and `dphi` (only even j, k occur). Arguments promote with the table type,
so BigFloat offsets through a BigFloat table stay full-precision.
"""
function puncture_phi(p::Puncture, dr::Real, dtheta::Real, dphi::Real)
    (; A006, A008, A024, A026, A042, A044, A060, A062, A080,
       A106, A108, A124, A126, A142, A144, A160, A162, A180,
       A204, A206, A222, A224, A240, A242, A260,
       A304, A306, A322, A324, A340, A342, A360,
       A402, A404, A420, A422, A440,
       A502, A504, A520, A522, A540,
       A600, A602, A620, A700, A702, A720, A800, A900) = p

    # local powers (kerr-circular.c:91-104)
    dr2 = dr * dr
    dr3 = dr2 * dr
    dr4 = dr2 * dr2
    dr6 = dr3 * dr3
    dr8 = dr4 * dr4
    dtheta2 = dtheta * dtheta
    dtheta4 = dtheta2 * dtheta2
    dtheta8 = dtheta4 * dtheta4
    sindphi = sin(dphi / 2)          # C: sin(0.5*dphi) — this is Q of eq:puncture
    sindphi2 = sindphi * sindphi
    sindphi4 = sindphi2 * sindphi2
    sindphi8 = sindphi4 * sindphi4

    # A = Σ A_ijk Δr^i Δθ^j Q^k, verbatim from kerr-circular.c:106-113
    A = (dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 +
         (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) +
            (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2) +
         (dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 +
            dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2))*sindphi2 +
         (A008 + A108*dr)*sindphi8 + sindphi4*(dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) +
            (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4 +
            (A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2)*sindphi2))

    # R² = ϱ² + z_c² Q² (eq:s2; kerr-circular.c:115-116 alpha/rho2) — scales from the orbit
    o = p.orbit
    alpha = o.rho2_dr2 * dr2 + o.rho2_dtheta2 * dtheta2   # ϱ², main.tex:452
    rho2 = alpha + o.zc2 * sindphi2                       # R², main.tex:447

    # Φ^P = A / R⁷ (C: A/pow(rho2, 3.5); written ^3·√ to stay exact under type promotion)
    return A / (rho2^3 * sqrt(rho2))
end
