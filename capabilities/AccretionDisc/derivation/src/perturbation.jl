# perturbation.jl — the linearised system of arXiv:2608.17003 Sec. II.B in the Fourier domain.
#
# Route R3 stage 2, part 1.  Everything is obtained by evaluating the SAME `fluid.jl` operators on
# λ-jets (jets.jl) built from
#
#     G_{μν} = g_{μν} + λ ĥ_{μν}(r) e^{−iωt+imφ}          (14), (45a), COVARIANT ĥ (ERRATUM E-M2)
#     U^α    = u^α    + λ δû^α(r) e^{−iωt+imφ}            (45a)
#     p → p + λ δp̂ e^{…},   ε → ε + λ δε̂ e^{…}
#
# so that the perturbed Christoffels, the perturbed projector q^{μν} and the perturbed
# u^ν∇_νu^μ are produced by the arithmetic — nothing is hand-derived.  The O(λ) parts of the Euler
# equation (3c), the energy equation (3b) and the normalisation are then Eqs. (15c), (15b) and the
# first step of the reduction scheme (M:604).
#
# Equatorial reflection symmetry of the background and of a circular equatorial orbit kills the
# z-mixing components ĥ_tz = ĥ_rz = ĥ_φz = 0 (route-R2 convention, kept here); ĥ_tr and ĥ_rφ ARE kept
# symbolic so that the "raw" source can be reported — they carry an explicit i and drop out for a
# circular orbit.

"""
    HMODE

The covariant Fourier metric-perturbation amplitude ĥ_{μν}(r) in (t,r,φ,z).
"""
const HMODE = Num[htt     htr      htphi     0
                  htr     hrr      hrphi     0
                  htphi   hrphi    hphiphi   0
                  0       0        0         hzz]

"δû^α = (δû^t, δû^r, δû^φ, 0) — the Eulerian velocity amplitude (no z-motion in the plane)."
const DUUP = Num[dut, dur, duphi, 0]

"G_{μν} = g_{μν} + λ ĥ_{μν} as a jet."
const GJET = [Jet(METRIC[μ, ν], HMODE[μ, ν]) for μ in 1:4, ν in 1:4]

"G^{μν} = g^{μν} − λ g^{μα}ĥ_{αβ}g^{βν} (the exact O(λ) inverse) as a jet."
const GINVJET = let dginv = -METRIC_INV * HMODE * METRIC_INV
    [Jet(METRIC_INV[μ, ν], dginv[μ, ν]) for μ in 1:4, ν in 1:4]
end

"U^α as a jet."
const UJET = [Jet(UUP[α], DUUP[α]) for α in 1:4]
"p as a jet."
const PJET = Jet(p, dp)
"ε as a jet."
const EJET = Jet(eps, deps)

"∂_α G_{μν}."
const DGJET = metric_partials(GJET, ∂jet)
"Γ^λ_{μν}[G] to O(λ) — the perturbed connection, never written by hand."
const GAMMAJET = christoffel_from_metric(GINVJET, DGJET)

"ℰ^μ (3c) as a jet; its O(λ) part is the linearised Euler equation (15c) with the index raised."
const EULER_JET = euler_up(GINVJET, GAMMAJET, UJET, PJET, EJET, ∂jet)
"Energy equation (3b) as a jet; its O(λ) part is (15b)."
const ENERGY_JET = energy_eq(GAMMAJET, UJET, PJET, EJET, ∂jet)
"g_{μν}u^μu^ν + 1 as a jet; its O(λ) part is 2u_αδû^α + ĥ_{αβ}u^αu^β (M:604)."
const NORM_JET = norm_constraint(GJET, UJET)

"δℰ^μ — the O(λ) Euler equation, index up."
const DEULER_UP = pert.(EULER_JET)
"δℰ_α = g_{αμ}δℰ^μ + ĥ_{αμ}ℰ^μ|₀ — the O(λ) Euler equation with the index DOWN (15c)."
const DEULER_LOW = [sum(METRIC[α, μ] * DEULER_UP[μ] + HMODE[α, μ] * bg(EULER_JET[μ]) for μ in 1:4)
                    for α in 1:4]
"O(λ) energy equation (15b)."
const DENERGY = pert(ENERGY_JET)
"O(λ) normalisation constraint."
const DNORM = pert(NORM_JET)

"""
    DUPHI_LOW

The COVARIANT azimuthal velocity perturbation δû_φ = g_{φμ}δû^μ + ĥ_{φμ}u^μ — the quantity the paper
prints in (50a) and consumes in the AMF (72a).  ERRATUM E-M1: (50a) as printed is the ORTHONORMAL
component δu_φ̂ = δu_φ/r, i.e. the Newtonian limit of this expression equals r × (50a).
"""
const DUPHI_LOW = sum(METRIC[IPHI, μ] * DUUP[μ] + HMODE[IPHI, μ] * UUP[μ] for μ in 1:4)

"""
    XI_R

Radial Lagrangian displacement ξ^r = i δû^r/(u^t ω_d), obtained by inverting (16): with u^r = 0 the
projector q^r{}_β = δ^r{}_β, so δû^r = 𝓛_uξ^r|_amplitude = u^t(−iω + imΩ)ξ^r = −i u^t ω_d ξ^r.
(Certified in `stage2_certificates` from (16) itself, not assumed.)
"""
const XI_R = im * dur / (UT * (omega - m * Ω_expr))

"""
    lie_u_xi_r(ξr)

𝓛_uξ^r for an amplitude ξ^r with the ansatz (45a): u^μ∂_μξ^r − ξ^μ∂_μu^r = u^t(−iω + imΩ)ξ^r
(u^r = 0 and the background is stationary/axisymmetric).
"""
lie_u_xi_r(ξr) = UT * (-im * omega + im * m * Ω_expr) * ξr
