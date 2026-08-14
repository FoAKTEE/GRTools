# Kerr equatorial thin-disc background (Dyson & D'Orazio arXiv:2601.19123, M=1, BL coords).
# The κ²-from-geodesic-SHEAR radial epicyclic frequency is PHYSICS — independent of the metric
# perturbation 𝒲 — and is the foundation the resonance/reduced-ODE leaves build on. Type-generic
# {T} so BigFloat precision (256..2048-bit) is inherited. From DiscModes/relativistic_disc.jl.
# Phase7-D dedup: delta_r / risco_disc are thin disc-struct delegators to the canonical Kerr
# primitives Spacetime.delta_kerr / Spacetime.kerr_isco (identical M=1 convention & sign branch;
# verified bit-exact at a ∈ {-0.9, 0, 0.6, 0.9, 0.99₂₅₆} before collapsing).

import Spacetime: delta_kerr, kerr_isco

"""
    RelativisticDisc(; a, h0, r_I, eden_I)

Kerr equatorial disc background (Dyson & D'Orazio 2026). `a`=spin (|a|<1, M=1), `h0`=constant
dimensionless aspect ratio (A9), `r_I`/`eden_I`=reference radius/energy-density for the background
ē_{(0,1)}=ē_I (r/r_I)^{-2}.
"""
Base.@kwdef struct RelativisticDisc{T<:Real}
    a::T      = 0.0
    h0::T     = 0.1
    r_I::T    = 1.0
    eden_I::T = 1.0
end

"Kerr horizon function Δ_r = r² − 2r + a² (M=1, equatorial) — delegates to `Spacetime.delta_kerr`."
delta_r(d::RelativisticDisc, r) = delta_kerr(d.a, r)

"Equatorial geometric circular frequency Ω^k = 1/(a + r^{3/2}) (eq:backgroundVel)."
disc_OmegaK(d::RelativisticDisc, r) = 1 / (d.a + r^(3//2))

"""
    kappar2(d, r)

Kerr radial epicyclic frequency κ_r² from the geodesic shear (Torok & Stuchlik 2005 closed form).
PHYSICS — independent of the perturbation Wronskian. Matches `KerrGeoFrequencies(a,r,0,1).Omegar²`.
"""
kappar2(d::RelativisticDisc, r) =
    (-3 * d.a^2 * sqrt(r) + 8 * d.a * r + (r - 6) * r^(3//2)) /
    (r^(5//2) * (d.a + r^(3//2))^2)

"""
    csound(d, r)

Sound speed from the constant-aspect-ratio thermodynamics (inverted eq:hdefinition): for fixed
dimensionless height `h0`, c_s(r) = h0 √(num/den) with the Kerr metric coefficients below.
"""
function csound(d::RelativisticDisc, r)
    h = d.h0
    num = 4 * d.a * sqrt(r) - 3 * d.a^2 - r^2
    den = 3 * d.a^2 * h^2 + (3 + h^2 - r) * r^2 - 2 * d.a * sqrt(r) * (2 * h^2 + r)
    h * sqrt(num / den)
end

"Background equatorial energy density ē_{(0,1)} = ē_I (r/r_I)^{-2} (eq:background_eofr)."
eden0(d::RelativisticDisc, r) = d.eden_I * (r / d.r_I)^(-2)

"""
    risco_disc(d) -> r_ISCO

Kerr equatorial prograde ISCO radius (Bardeen 1972) — the inner edge / resonance-search inner
bound. a=0 ⇒ 6; the sign branch selects prograde (a>0) vs retrograde (a<0).
Delegates to the canonical `Spacetime.kerr_isco` (identical closed form and sign convention).
"""
risco_disc(d::RelativisticDisc) = kerr_isco(d.a)
