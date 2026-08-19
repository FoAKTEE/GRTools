# SpiralWaves/test_background.jl — EGD::sw-background (arXiv:2608.17003 Sec. II.A, main.tex:121-205).
# Tests FIRST for `background.jl`: the `DiscBackground{T}` equilibrium (power-law ε, p, c_s (12a-c),
# Ω•(a•) (13), Γ1 (22), the relativistic hydrostatic balance (11) solved EXACTLY for v and Ω, the
# four-velocity (8)-(10), κ², r_ISCO, F_{J,0}/F_{J,1} (76)). GROUNDING:
#   • (11) residual: p'/(ε+p) + (ν' − v²ψ' + e^{ψ−ν} v ω̃')/(1−v²) ≤ 1e-13 at 60 deterministic
#     pseudo-random (χ, a•, h0, 𝔭, 𝔮, r) points — the metric potentials and their derivatives come
#     from Spacetime.kerr_equatorial_potentials(_dr) (EGD::spacetime-kerr-equatorial-potentials,
#     solid) and p'/(ε+p) is recomputed HERE from the power laws (12a-c).
#   • p' → 0 (𝔭 = 𝔮 = 0): Ω(r) == the circular geodesic frequency Modes.disc_OmegaK == 1/(χ + r^{3/2})
#     to 1e-14, u^t == Spacetime.kerr_circular_dtdtau to 1e-14 (prograde AND χ<0).
#   • a = 0 Newtonian limit: Ω² r³ − 1 = −(𝔭+𝔮) h0² (1 + O(h0², 1/r)) — scales as h0².
#   • BigFloat(256) inherits: residual ≤ 1e-60, every output is BigFloat; @inferred type stability.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.Modes: disc_OmegaK, kappar2, RelativisticDisc
using Spacetime: kerr_equatorial_potentials, kerr_equatorial_potentials_dr, kerr_circular_dtdtau,
                 kerr_isco, kerr_circular_ut_uphi

# deterministic pseudo-random stream (no Random stdlib in the test target): LCG in [0,1)
function _lcg_stream(seed::UInt64, n::Int)
    out = Vector{Float64}(undef, n); s = seed
    for i in 1:n
        s = s * 0x5851f42d4c957f2d + 0x14057b7ef767814f
        out[i] = (s >> 11) / 2.0^53
    end
    out
end

# independent (11) residual: the hydrostatic balance with p'/(ε+p) from the power laws (12a-c):
# p ∝ r^{-(𝔭+𝔮)} ⇒ p' = −(𝔭+𝔮) p / r, ε + p = ε (1 + c_s²)
function _balance_residual(bg, r)
    T = typeof(r)
    ε = eps_profile(bg, r); p = pressure(bg, r)
    dp = -(bg.𝔭 + bg.𝔮) * p / r
    lhs = dp / (ε + p)
    _, ψ, _, _ = kerr_equatorial_potentials(bg.χ, r)
    ν = kerr_equatorial_potentials(bg.χ, r).ν
    d = kerr_equatorial_potentials_dr(bg.χ, r)
    rp = rotation_profile(bg, r)
    v = rp.v
    rhs = -(d.dν - v^2 * d.dψ + exp(ψ - ν) * v * d.dω̃) / (1 - v^2)
    (lhs - rhs) / max(abs(d.dν), one(T))
end

@testset "SpiralWaves.background (arXiv:2608.17003 Sec. II.A)" begin

    @testset "DiscBackground construction, defaults, promotion, validation" begin
        bg = DiscBackground(χ = 0.9, aEMRI = 15.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
        @test bg isa DiscBackground{Float64}
        @test bg.b == 0.6 && bg.εp == 1e-5 && bg.σε == 1.5 && bg.ε0 == 1.0      # M:735, M:884 nominal
        @test bg.χ == 0.9 && bg.aEMRI == 15.0 && bg.h0 == 0.05 && bg.𝔭 == 1.5 && bg.𝔮 == 1.0
        bgi = DiscBackground(χ = 0, aEMRI = 15, h0 = 0.05)                        # Int promote, 𝔭=𝔮=0 default
        @test bgi isa DiscBackground{Float64} && bgi.𝔭 == 0 && bgi.𝔮 == 0
        bgb = DiscBackground(χ = big"0.9", aEMRI = 15, h0 = 0.05)
        @test bgb isa DiscBackground{BigFloat} && bgb.σε isa BigFloat && bgb.σε == big"1.5"
        bgb2 = DiscBackground{BigFloat}(χ = 0.9, aEMRI = 15, h0 = 0.05)
        @test bgb2 isa DiscBackground{BigFloat}
        @test_throws DomainError DiscBackground(χ = 0.0, aEMRI = 5.0, h0 = 0.05)    # inside the ISCO
        @test_throws DomainError DiscBackground(χ = 1.2, aEMRI = 15.0, h0 = 0.05)   # |χ| > 1
        @test_throws DomainError DiscBackground(χ = 0.0, aEMRI = 15.0, h0 = 0.0)    # h0 ≤ 0
        @test_throws DomainError DiscBackground(χ = 0.0, aEMRI = 15.0, h0 = 0.05, b = -0.1)
        @test r_isco(bg) == kerr_isco(0.9)
        @test softening_length(bg) == 15.0 * 0.05 * 0.6                          # (59) ε_soft = a• h0 b
    end

    @testset "power laws (12a-c), Ω• (13), Γ1 (22)" begin
        bg = DiscBackground(χ = 0.6, aEMRI = 20.0, h0 = 0.1, 𝔭 = 1.5, 𝔮 = 1.0, ε0 = 2.0)
        Om = Omega_emri(bg)
        @test Om == 1 / (20.0^1.5 + 0.6)                                            # (13), M = 1
        @test Om == disc_OmegaK(bg, 20.0)                                          # one name, one concept
        @test disc_OmegaK(bg, 8.0) == disc_OmegaK(RelativisticDisc(a = 0.6), 8.0)
        @test eps_profile(bg, 20.0) == 2.0                                         # ε(a•) = ε0
        @test eps_profile(bg, 40.0) ≈ 2.0 * 2.0^(-1.5) rtol = 1e-15
        @test eden0(bg, 40.0) == eps_profile(bg, 40.0)                             # generic name → same concept
        @test csound(bg, 20.0) ≈ 0.1 * 20.0 * Om rtol = 1e-15                      # c_s(a•) = h0 a• Ω•
        @test csound(bg, 80.0) ≈ csound(bg, 20.0) * 4.0^(-0.5) rtol = 1e-15        # (r/a•)^{-𝔮/2}
        @test pressure(bg, 30.0) ≈ csound(bg, 30.0)^2 * eps_profile(bg, 30.0) rtol = 1e-15
        @test gamma1(bg, 30.0) ≈ (eps_profile(bg, 30.0) + pressure(bg, 30.0)) / pressure(bg, 30.0) *
                                 csound(bg, 30.0)^2 rtol = 1e-15
        @test gamma1(bg, 30.0) ≈ 1 + csound(bg, 30.0)^2 rtol = 1e-15                # exact consequence of (12b)
        # p'/(ε+p) closed form vs finite difference of pressure (Float64 FD, 1e-7 is the FD floor)
        r = 25.0; h = 1e-4
        fd = (pressure(bg, r + h) - pressure(bg, r - h)) / (2h) / (eps_profile(bg, r) + pressure(bg, r))
        @test dlogp_dr_over_enthalpy(bg, r) ≈ fd rtol = 1e-7
        # normalisations (76): F_{J,0} = ε0 a•⁴ Ω•² q²/h0³, F_{J,1} = F_{J,0}/a•
        @test F_J0(bg; q = 1e-4) ≈ 2.0 * 20.0^4 * Om^2 * 1e-8 / 0.1^3 rtol = 1e-15
        @test F_J1(bg; q = 1e-4) ≈ F_J0(bg; q = 1e-4) / 20.0 rtol = 1e-15
        @test F_J0(bg) == F_J0(bg; q = 1)
    end

    @testset "(11) residual ≤ 1e-13 at 60 pseudo-random (χ, a•, h0, 𝔭, 𝔮, r)" begin
        u = _lcg_stream(0x2608170300000001, 6 * 60)
        worst = 0.0
        for i in 1:60
            χ  = -0.95 + 1.9 * u[6i-5]
            ris = kerr_isco(χ)
            ap = ris * (1.2 + 8.0 * u[6i-4])        # a• ∈ [1.2, 9.2] r_ISCO
            h0 = 0.02 + 0.18 * u[6i-3]
            𝔭  = 3.0 * u[6i-2]
            𝔮  = 2.0 * u[6i-1]
            r  = ris * (1.02 + 6.0 * u[6i])
            bg = DiscBackground(χ = χ, aEMRI = ap, h0 = h0, 𝔭 = 𝔭, 𝔮 = 𝔮)
            res = abs(_balance_residual(bg, r))
            worst = max(worst, res)
            rp = rotation_profile(bg, r)
            @test 0 < rp.v < 1
            @test rp.ut > 0 && rp.uφ > 0 && rp.ut_low < 0
            # u·u = −1 from (8)-(10): u^t u_t + u^φ u_φ = −1
            @test rp.ut * rp.ut_low + rp.Ω * rp.ut * rp.uφ ≈ -1 atol = 1e-13
            # v = (Ω − ω̃) e^{ψ−ν} (9)
            pot = kerr_equatorial_potentials(χ, r)
            @test rp.v ≈ (rp.Ω - pot.ω̃) * exp(pot.ψ - pot.ν) rtol = 1e-13
            @test Omega(bg, r) == rp.Ω
        end
        @test worst ≤ 1e-13
        println("  (11) worst scaled residual over 60 points: ", worst)
    end

    @testset "p' → 0 (𝔭 = 𝔮 = 0): Ω == Kerr circular geodesic frequency, u^t == dt/dτ" begin
        for χ in (-0.9, -0.5, 0.0, 0.6, 0.9, 0.998), rr in (1.05, 1.5, 3.0, 10.0)
            ris = kerr_isco(χ)
            bg = DiscBackground(χ = χ, aEMRI = 4ris, h0 = 0.1, 𝔭 = 0.0, 𝔮 = 0.0)
            r = rr * ris
            @test dlogp_dr_over_enthalpy(bg, r) == 0
            rp = rotation_profile(bg, r)
            @test rp.Ω ≈ disc_OmegaK(bg, r) rtol = 1e-14
            @test rp.Ω ≈ 1 / (χ + r^1.5) rtol = 1e-14
            @test rp.ut ≈ kerr_circular_dtdtau(χ, r) rtol = 1e-14
            @test rp.Ω * rp.ut ≈ kerr_circular_ut_uphi(χ, r).uphi rtol = 1e-14
        end
    end

    @testset "a = 0 Newtonian limit: Ω² r³ − 1 = −(𝔭+𝔮) h0² (1 + O(h0², 1/r))" begin
        ap = 1.0e3
        for (𝔭, 𝔮) in ((1.5, 0.0), (2.0, 1.0), (0.0, 1.0))
            δ(h0) = Omega(DiscBackground(χ = 0.0, aEMRI = ap, h0 = h0, 𝔭 = 𝔭, 𝔮 = 𝔮), ap)^2 * ap^3 - 1
            d1 = δ(0.05); d2 = δ(0.025)
            @test d1 < 0                                   # pressure support ⇒ sub-Keplerian
            @test d1 / d2 ≈ 4 rtol = 2e-2                  # O(h0²) scaling
            @test d1 ≈ -(𝔭 + 𝔮) * 0.05^2 rtol = 1e-2       # leading Newtonian coefficient
            @test abs(d1 + (𝔭 + 𝔮) * 0.05^2) ≤ 5e-3 * abs(d1)
        end
        # 𝔭 = 𝔮 = 0 ⇒ exactly Keplerian in BL coordinate time (a = 0)
        @test Omega(DiscBackground(χ = 0.0, aEMRI = ap, h0 = 0.05), ap)^2 * ap^3 ≈ 1 rtol = 1e-14
    end

    @testset "κ² (geodesic radial epicyclic frequency) and ω_d" begin
        bg = DiscBackground(χ = 0.7, aEMRI = 12.0, h0 = 0.05, 𝔭 = 1.5)
        @test kappa2(bg, 9.0) == kappar2(RelativisticDisc(a = 0.7), 9.0)          # same concept, delegated
        @test kappa2(bg, 9.0) == kappar2(bg, 9.0)
        bg0 = DiscBackground(χ = 0.0, aEMRI = 1.0e4, h0 = 0.05)
        @test kappa2(bg0, 1.0e4) * 1.0e12 ≈ 1 - 6e-4 rtol = 1e-12                  # (r−6)/r⁴ at a = 0
        # ω_d = ω − mΩ = m(Ω• − Ω(r)) (47), (51): vanishes at r = a• when 𝔭 = 𝔮 = 0
        bgk = DiscBackground(χ = 0.7, aEMRI = 12.0, h0 = 0.05)
        @test abs(omega_d(bgk, 3, 12.0)) ≤ 1e-15
        @test omega_d(bgk, 3, 20.0) > 0 && omega_d(bgk, 3, 8.0) < 0
        @test omega_d(bgk, 3, 20.0) ≈ 3 * (Omega_emri(bgk) - Omega(bgk, 20.0)) rtol = 1e-15
    end

    @testset "BigFloat(256) inherits: (11) residual ≤ 1e-60" begin
        setprecision(BigFloat, 256) do
            bg = DiscBackground(χ = big"0.9", aEMRI = big"15", h0 = big"0.05", 𝔭 = big"1.5", 𝔮 = big"1")
            @test bg isa DiscBackground{BigFloat}
            for r in (big"3", big"7.5", big"15", big"40")
                rp = rotation_profile(bg, r)
                @test rp.Ω isa BigFloat && rp.v isa BigFloat && rp.ut isa BigFloat
                @test abs(_balance_residual(bg, r)) ≤ big"1e-60"
            end
            @test csound(bg, big"15") isa BigFloat && gamma1(bg, big"15") isa BigFloat
            @test Omega_emri(bg) isa BigFloat && F_J0(bg) isa BigFloat
            @test kappa2(bg, big"15") isa BigFloat && omega_d(bg, 2, big"15") isa BigFloat
            # pressureless branch in BigFloat: geodesic to 1e-70
            bg0 = DiscBackground(χ = big"-0.9", aEMRI = big"20", h0 = big"0.1")
            @test abs(rotation_profile(bg0, big"15").Ω - 1 / (big"-0.9" + big"15"^big"1.5")) ≤ big"1e-70"
        end
    end

    @testset "type stability (@inferred)" begin
        bg = DiscBackground(χ = 0.5, aEMRI = 15.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
        @test @inferred(rotation_profile(bg, 12.0)) isa NamedTuple
        @test @inferred(Omega(bg, 12.0)) isa Float64
        @test @inferred(csound(bg, 12.0)) isa Float64
        @test @inferred(pressure(bg, 12.0)) isa Float64
        @test @inferred(eps_profile(bg, 12.0)) isa Float64
        @test @inferred(gamma1(bg, 12.0)) isa Float64
        @test @inferred(kappa2(bg, 12.0)) isa Float64
        @test @inferred(omega_d(bg, 4, 12.0)) isa Float64
        @test @inferred(Omega_emri(bg)) isa Float64
        @test @inferred(F_J0(bg)) isa Float64
        @test @inferred(r_isco(bg)) isa Float64
        bgb = DiscBackground(χ = big"0.5", aEMRI = 15, h0 = 0.05, 𝔭 = 1.5)
        @test @inferred(rotation_profile(bgb, big"12")) isa NamedTuple
        @test @inferred(csound(bgb, big"12")) isa BigFloat
    end
end
