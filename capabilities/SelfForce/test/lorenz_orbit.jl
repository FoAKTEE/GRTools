# Grounded tests for SelfForce.LorenzField circular-orbit parameters + BLS source coefficients.
#
# GROUNDING — fully ANALYTIC (Schwarzschild closed forms) + core/Spacetime characterization:
#
# (1) orbit_energy/orbit_ang_mom/orbit_phi_freq are RE-EXPORTS of Spacetime.schw_orbit_* (committed
#     green). We assert object IDENTITY (===) and ground the VALUES against the libSch.c circular
#     closed forms E=(1-2/p)/√(1-3/p), L=√p/√(1-3/p), Ωφ=p^(-3/2), plus the circular-orbit identity
#     E² = f(p)(1 + L²/p²) at the radial turning point (contra_ur=0). All M=1.
#
# (2) R1..R10 source coefficients (χ=π/2, SOURCE_TOTAL of ref-code/h1Lorenz/src/sources.c): the
#     circular turning-point closed forms are self-contained analytic literals. We assert
#       R5_source ≡ R9_source ≡ 0  (∝ source_sin / contra_ur, both vanish for e=0),
#     and the explicit -32π·P̄_lm·(…) forms at sample (l,m,p), with the normalized-Legendre factor
#     P̄ cross-checked against NumericsCore.assoc_legendre_P(;normalized=true) (e.g. P̄_22(π/2)=
#     √(15/(32π))) and the θ-derivative factor against NumericsCore.legendre_P_thetaderiv. R7=
#     (l(l+1)-2m²)·R6 is checked as the exact algebraic relation.

using Test
using SelfForce.LorenzField: OrbitalParams, orbit_energy, orbit_ang_mom, orbit_phi_freq,
                             R1_source, R3_source, R5_source, R6_source, R7_source,
                             R9_source, R10_source
using NumericsCore: assoc_legendre_P, legendre_P_thetaderiv
import Spacetime

@testset "orbit constants: Spacetime re-export identity + closed forms" begin
    @test orbit_energy === Spacetime.schw_orbit_energy
    @test orbit_ang_mom === Spacetime.schw_orbit_ang_mom
    @test orbit_phi_freq === Spacetime.schw_orbit_phi_freq
    for p in (6.5, 8.0, 10.0, 20.0, 50.0)
        E = orbit_energy(p); L = orbit_ang_mom(p); Ω = orbit_phi_freq(p)
        @test E ≈ (1 - 2 / p) / sqrt(1 - 3 / p) rtol = 1e-14
        @test L ≈ sqrt(p) / sqrt(1 - 3 / p)     rtol = 1e-14
        @test Ω ≈ p^(-1.5)                       rtol = 1e-14
        f = 1 - 2 / p
        # circular turning-point identity: contra_ur² = E² - f(1 + L²/p²) = 0
        @test E^2 ≈ f * (1 + L^2 / p^2) rtol = 1e-13
    end
end

@testset "OrbitalParams populates E/L/Ω at p (type-generic)" begin
    p = 8.0
    o = OrbitalParams(p, [p], 0, 0, 0)
    @test o.p == p
    @test o.E ≈ orbit_energy(p) rtol = 1e-14
    @test o.L ≈ orbit_ang_mom(p) rtol = 1e-14
    @test o.Omega_phi ≈ orbit_phi_freq(p) rtol = 1e-14
    setprecision(BigFloat, 256) do
        ob = OrbitalParams(big"8.0", [big"8.0"], 0, 0, 0)
        @test ob.E isa BigFloat
        @test ob.E ≈ (1 - 2 / big"8.0") / sqrt(1 - 3 / big"8.0") rtol = 1e-60
    end
end

@testset "R-sources vanish for circular orbit (R5≡R9≡0) + Legendre cross-check" begin
    for (l, m, p) in ((2, 2, 8.0), (3, 2, 10.0), (2, 0, 12.0))
        o = OrbitalParams(p, [p], 0, 0, 0)
        @test R5_source(l, m, o) == 0
        @test R9_source(l, m, o) == 0
    end
    # P̄_22(π/2) = √(15/(32π)); core normalized Legendre is the source's P̄ factor
    @test assoc_legendre_P(2, 2, π / 2; normalized = true) ≈ sqrt(15 / (32π)) rtol = 1e-13
    @test assoc_legendre_P(2, 1, π / 2; normalized = true) == 0.0   # (l+m) odd parity
end

@testset "R1/R3/R6 explicit χ=π/2 closed forms at (l,m,p)" begin
    for (l, m, p) in ((2, 2, 8.0), (4, 2, 10.0), (3, 1, 12.0))
        o = OrbitalParams(p, [p], 0, 0, 0)
        E = orbit_energy(p); L = orbit_ang_mom(p); fp = 1 - 2 / p
        P̄ = assoc_legendre_P(l, m, π / 2; normalized = true)
        c8 = -4 * 8 * π * P̄
        R1exp = c8 * fp^2 * (2 * E^2 * p^2 - fp * p^2 - L^2 * fp) / (E * p^3)
        R3exp = c8 * fp^2 * (p^2 + L^2) / (E * p^3)
        R6exp = (-4 * 8 * π * L^2 / E * P̄) * fp^2 / p^3
        @test real(R1_source(l, m, o)) ≈ R1exp rtol = 1e-13
        @test imag(R1_source(l, m, o)) == 0
        @test real(R3_source(l, m, o)) ≈ R3exp rtol = 1e-13
        @test real(R6_source(l, m, o)) ≈ R6exp rtol = 1e-13
    end
end

@testset "R7 = (l(l+1)-2m²)·R6 + R10 carries θ-derivative·I" begin
    for (l, m, p) in ((2, 2, 8.0), (3, 2, 10.0), (4, 1, 9.0))
        o = OrbitalParams(p, [p], 0, 0, 0)
        @test R7_source(l, m, o) ≈ (l * (l + 1) - 2 * m^2) * R6_source(l, m, o) rtol = 1e-13
        E = orbit_energy(p); L = orbit_ang_mom(p); fp = 1 - 2 / p
        Pθ = legendre_P_thetaderiv(l, m, π / 2; normalized = true)
        R10exp = im * (-4 * 16 * π * m * L^2 / E * Pθ) * fp^2 / p^3
        @test R10_source(l, m, o) ≈ R10exp rtol = 1e-13
        @test real(R10_source(l, m, o)) == 0   # purely imaginary (explicit I)
    end
end

@testset "R-sources BigFloat precision path" begin
    setprecision(BigFloat, 256) do
        o = OrbitalParams(big"8.0", [big"8.0"], 0, 0, 0)
        r1 = R1_source(2, 2, o)
        @test r1 isa Complex{BigFloat}
        E = orbit_energy(big"8.0"); L = orbit_ang_mom(big"8.0"); fp = 1 - 2 / big"8.0"
        P̄ = assoc_legendre_P(2, 2, big(π) / 2; normalized = true)
        R1exp = (-4 * 8 * big(π) * P̄) * fp^2 * (2 * E^2 * 64 - fp * 64 - L^2 * fp) / (E * big"8.0"^3)
        @test real(r1) ≈ R1exp rtol = 1e-60
    end
end
