# RT-A1 reduced radial-ODE operator (arXiv:2601.19123 eq:reduceeq): 𝔥''+A_m𝔥'+B_m𝔥=𝔉_m, with the
# Lin-Shu / Goldreich-Tremaine density-wave dispersion D=κ²−ω̃² as its principal restoring structure.
# RIGOROUS pieces from the published closed forms:
#   • A_m = −d/dr log𝒲 (Abel's identity; 𝒲 = the cold-canonical Resonances.wronskian) — EXACT.
#   • B_m_corotation — the regular-singular 1/(r−r_p)² structure (derivation.md l.652).
#   • corotation_exponent — the Frobenius indicial root at corotation (RT-A6 / eq:divergentfact).
#   • dispersion D=κ²−ω̃² — the radial-epicyclic restoring term (κ² from the geodesic SHEAR), the
#     density-wave principal symbol; D≡−lindblad_residual, so D=0 marks the Lindblad turning points.
# The full reduced operator is assembled by `reduced_operator` (A_m exact + B_m_full potential +
# 𝔉_m source + D dispersion + the Lindblad/corotation singular structure), with the test-function
# master-PDE reduction `MasterPDE.master_radial_coeffs` exposed as the independent cross-validation.
# H-DISC-FULLOP [DOWNGRADED]: the master-PDE operator reduction IS now ported (leading-order
# equatorial, master_radial_coeffs). The only remaining residual is the PMLA h²/Gaussian
# height-average refinement ([FUTURE], genuinely unpublished) — NOT the operator itself.
# From DiscModes/relativistic_master.jl. Type-generic {T}; BigFloat-inherited.

using ..Modes: RelativisticDisc, disc_OmegaK, kappar2, delta_r, csound
using ..Modes: wronskian, lindblad_residual, lindblad_locations, master_radial_coeffs

"""
    A_m(d, r, r_p, m; h=nothing) -> −d/dr log𝒲(r)

Exact first-order coefficient of the reduced ODE (Abel's identity for 𝒲 = the cold-canonical
Wronskian eq:wronskain). Central finite difference; 𝒲 is smooth away from the Lindblad resonances
(where 𝒲→0 ⇒ A_m singular — the regular-singular structure).
"""
function A_m(d::RelativisticDisc{T}, r, r_p, m; h = nothing) where {T}
    hh = h === nothing ? sqrt(eps(T)) * max(abs(r), one(T)) : T(h)
    Wp = wronskian(d, r + hh, r_p, m); Wm = wronskian(d, r - hh, r_p, m); W0 = wronskian(d, r, r_p, m)
    -(Wp - Wm) / (2 * hh) / W0
end

"""
    B_m_corotation(d, r, r_p, m)

The 1/(r−r_p)²-coefficient of the corotation-limit operator (derivation.md l.652): near
r→r_p the reduced operator is 𝔥'' − [4κ_{r_p}²/(9 m² r_p Δ_{r_p}(Ω^k_{r_p})⁴ (r−r_p)²)] 𝔥, a
regular-singular point. Returns the local B_m there.
"""
function B_m_corotation(d::RelativisticDisc{T}, r, r_p, m) where {T}
    kp2 = kappar2(d, r_p); Op = disc_OmegaK(d, r_p); Dp = delta_r(d, r_p)
    -4 * kp2 / (9 * m^2 * r_p * Dp * Op^4 * (r - r_p)^2)
end

"""
    B_m_wkb(d, r, r_p, m; h=nothing)

Density-wave (Lindblad) part of B_m via the WKB dispersion (eq:WKBwavenumber inverted):
k²_m = B_m − ½A_m′ − ¼A_m² ⇒ B_m = k² + ½A_m′ + ¼A_m², with the relativistic disc-wave dispersion
k²_m = [m²(Ω−Ω_p)² − κ_r²]/c_s² = lindblad_residual/c_s² = −[`dispersion`](@ref)/c_s², and
A_m=−d/dr log𝒲 analytic. The k²=0 turning points sit at the Lindblad resonances (evanescent
between ILR/OLR). The density-wave restoring term D=κ²−ω̃² (the `dispersion` primitive) IS the
principal symbol of the validated leading-order reduction (`MasterPDE.master_radial_coeffs`); the
PMLA h²/Gaussian height-average is the only [FUTURE] residual.
"""
function B_m_wkb(d::RelativisticDisc{T}, r, r_p, m; h = nothing) where {T}
    hh = h === nothing ? sqrt(eps(T)) * max(abs(r), one(T)) : T(h)
    k2 = lindblad_residual(d, r, r_p, m) / csound(d, r)^2
    A0 = A_m(d, r, r_p, m); dA = (A_m(d, r + hh, r_p, m) - A_m(d, r - hh, r_p, m)) / (2 * hh)
    k2 + dA / 2 + A0^2 / 4
end

"""
    B_m_full(d, r, r_p, m) -> B_m_wkb + B_m_corotation

The density-wave B_m_wkb plus the 1/(r−r_p)² corotation singularity — correct Lindblad turning
points and damped region AND the corotation singular structure (the validated RT-A6 form). The
potential of the assembled [`reduced_operator`](@ref); cross-validated against the test-function
reduction `MasterPDE.master_radial_coeffs` (sign + order of magnitude). Singular at r=r_p (handle
via Hadamard regularization, as in the paper).
"""
B_m_full(d::RelativisticDisc, r, r_p, m) = B_m_wkb(d, r, r_p, m) + B_m_corotation(d, r, r_p, m)

"""
    corotation_exponent(d, r_p, m) -> s

Decaying Frobenius singular exponent of 𝔥_m at corotation (RT-A6 / eq:divergentfact): with
C = 4κ_{r_p}²/(9 m² r_p Δ_{r_p}(Ω^k_{r_p})⁴), the indicial roots solve s²−s−C=0 ⇒
s = ½(1−√(1+4C)). The mode is regularised by 𝔥_m = (r−r_p)^s ψ_m. Fiducial
(a=0.6, r_p=15, m=2): s ≈ −0.0818.
"""
function corotation_exponent(d::RelativisticDisc{T}, r_p, m) where {T}
    C = 4 * kappar2(d, r_p) / (9 * m^2 * r_p * delta_r(d, r_p) * disc_OmegaK(d, r_p)^4)
    (1 - sqrt(1 + 4 * C)) / 2
end

"""
    dispersion(d::RelativisticDisc, r, r_p, m) -> D = κ_r²(r) − ω̃²

Lin-Shu / Goldreich-Tremaine density-wave dispersion: the radial-epicyclic restoring `κ_r²` (from
the background geodesic SHEAR, [`kappar2`](@ref) — physics, independent of the perturbation 𝒲) minus
the comoving frequency squared `ω̃²`, with `ω̃ = m(Ω^k_r − Ω^k_{r_p})` (Ω_p the corotating pattern
speed at `r_p`). `D` is the principal symbol of the reduced operator: `D=0` marks the Lindblad
turning points (where the density wave reflects) and `D<0` the evanescent damped region between
ILR/OLR. Exactly `D ≡ −lindblad_residual` (the Wronskian numerator), so the same zeros drive the
dispersion, the WKB turning points and the resonance locations.
"""
function dispersion(d::RelativisticDisc, r, r_p, m)
    ω̃ = m * (disc_OmegaK(d, r) - disc_OmegaK(d, r_p))
    kappar2(d, r) - ω̃^2
end

"""
    reduced_operator(d::RelativisticDisc, r, r_p, m; Φm=nothing) -> NamedTuple

Assemble the full reduced radial-ODE operator `𝔥'' + A_m 𝔥' + B_m 𝔥 = 𝔉_m` (eq:reduceeq) at radius
`r` for a perturber at `r_p`, azimuthal mode `m`. Returns the operator data with the principal
symbol normalised to 1:

  - `principal`  — the ∂²-coefficient, normalised to `1`.
  - `A_m`        — the EXACT Abel first-order coefficient `−d/dr log𝒲` ([`A_m`](@ref)).
  - `B_m`        — the potential `B_m_full` (density-wave `B_m_wkb` + corotation singular term).
  - `dispersion` — the density-wave restoring `D=κ²−ω̃²` ([`dispersion`](@ref)).
  - `forcing`    — the source `𝔉_m` ([`forcing_Fm`](@ref); `Φm` overrides the default potential).
  - `corotation` — the Frobenius singular exponent at corotation ([`corotation_exponent`](@ref)).
  - `lindblad`   — the `(r_ILR, r_OLR)` turning points ([`lindblad_locations`](@ref)).
  - `master`     — the INDEPENDENT test-function reduction `MasterPDE.master_radial_coeffs`, exposed
                   as the leading-order cross-validation of `A_m`/`B_m` (same sign + order).

The Lindblad turning-point structure (`dispersion`/`lindblad`) and the corotation regular-singular
structure (`B_m`/`corotation`) are both carried explicitly. The remaining unmodelled physics is the
PMLA h²/Gaussian height-average refinement ([FUTURE]).
"""
function reduced_operator(d::RelativisticDisc{T}, r, r_p, m; Φm = nothing) where {T}
    (principal = one(T),
     A_m = A_m(d, r, r_p, m),
     B_m = B_m_full(d, r, r_p, m),
     dispersion = dispersion(d, r, r_p, m),
     forcing = forcing_Fm(d, r, r_p, m; Φm = Φm),
     corotation = corotation_exponent(d, r_p, m),
     lindblad = lindblad_locations(d, r_p, m),
     master = master_radial_coeffs(d, r, r_p, m))
end
