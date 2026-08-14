# Analytic Wronskian (eq:wronskain) + Lindblad/corotation resonance structure (RT-A4), from
# DiscModes/relativistic_disc.jl. The reduced-ODE Wronskian is the CANONICAL cold-residual form
# (the F6_QRES pressure-shift fork is collapsed away — AHA): its numerator is exactly the Lindblad
# residual m²(Ω−Ω_p)²−κ², so the same zeros drive the Wronskian, the WKB turning points, and the
# resonance locations. Root-finding reuses NumericsCore.bisect (the hand-rolled _bisect is dropped).
# Type-generic {T}; BigFloat-inherited.

using ..Modes: RelativisticDisc, disc_OmegaK, kappar2, delta_r, risco_disc
import NumericsCore: bisect

"""
    lindblad_residual(d, r, r_p, m)

Lindblad resonance residual m²(Ω^k_r − Ω^k_{r_p})² − κ_r² — its zeros are the inner/outer Lindblad
resonances (eq:wronskain zeros). This residual IS the Wronskian numerator.
"""
lindblad_residual(d::RelativisticDisc, r, r_p, m) =
    m^2 * (disc_OmegaK(d, r) - disc_OmegaK(d, r_p))^2 - kappar2(d, r)

"Corotation residual Ω^k_r − Ω^k_{r_p} (zero at r=r_p)."
corotation_residual(d::RelativisticDisc, r, r_p) = disc_OmegaK(d, r) - disc_OmegaK(d, r_p)

"""
    wronskian(d, r, r_p, m; Wp=1)

Reduced-ODE Wronskian (eq:wronskain), CANONICAL cold-residual form:

    𝒲(r) = W_p · [m²(Ω^k_r−Ω^k_{r_p})² − κ_r²] · √r / [(2a+(r−3)√r)(Ω^k_r)²(Ω^k_{r_p})² Δ_r].

The bracketed factor is exactly `lindblad_residual`, so 𝒲 vanishes precisely at the Lindblad
resonances. (F6_QRES pressure-shift fork collapsed — this is the validated analytic Wronskian.)
"""
function wronskian(d::RelativisticDisc, r, r_p, m; Wp = one(float(r)))
    Or = disc_OmegaK(d, r); Op = disc_OmegaK(d, r_p)
    num = lindblad_residual(d, r, r_p, m)
    Wp * num / ((2 * d.a + (r - 3) * sqrt(r)) * Or^2 * Op^2 * delta_r(d, r)) * sqrt(r)
end

"""
    lindblad_locations(d, r_p, m) -> (r_ILR, r_OLR)

Inner/outer Lindblad resonance radii for a perturber at `r_p` (mode `m`): the roots of
m(Ω^k_r − Ω^k_{r_p}) = ±√κ_r² (the relativistic Lindblad condition). Returns `NaN` for a root not
bracketed in the disc. (Same condition as Hirata 2010 — the link to the 1010.0759 results.)
"""
lindblad_locations(d::RelativisticDisc, r_p, m) = _lindblad_locations_impl(d, r_p, m)

function _lindblad_locations_impl(d::RelativisticDisc{T}, r_p, m) where {T}
    Op = disc_OmegaK(d, r_p)
    fILR(r) = m * (disc_OmegaK(d, r) - Op) - sqrt(kappar2(d, r))   # inner: m(Ω−Ω_p) = +κ
    fOLR(r) = m * (disc_OmegaK(d, r) - Op) + sqrt(kappar2(d, r))   # outer: m(Ω−Ω_p) = −κ
    ε = T(1e-6)
    risco = risco_disc(d)
    lo_i, hi_i = risco + ε, r_p - ε
    lo_o, hi_o = r_p + ε, 12 * r_p
    rILR = (fILR(lo_i) * fILR(hi_i) <= 0) ? bisect(fILR, lo_i, hi_i) : T(NaN)
    rOLR = (fOLR(lo_o) * fOLR(hi_o) <= 0) ? bisect(fOLR, lo_o, hi_o) : T(NaN)
    rILR, rOLR
end
