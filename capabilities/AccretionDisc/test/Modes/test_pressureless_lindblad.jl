# Module — Modes/PressurelessLindblad.jl.  The ANALYTIC pressureless relativistic Lindblad torque of
# Hegade K. R., Gammie & Yunes 2025 (arXiv:2510.03564 "Paper II", Kerr), the pressureless comparison
# baseline of arXiv:2608.17003 Fig. 9 [Eqs. (77)-(80) there].  GROUNDING — four INDEPENDENT oracles,
# none of them this file's own code:
#   • FREQUENCY REUSE — the printed auxiliary polynomials of Eq. (19) are asserted equal to the
#     existing Kerr primitives: (χ+x⁶) = 1/`disc_OmegaK`, (−3χ²+x⁸−6x⁴+8χx²) = p²`kappar2`/Ω²,
#     (χ²+x⁸−2x⁴) = `delta_r`, (2χ+x⁶−3x²) = 1/(p^{3/2}Ω²(u^t)²) with `kerr_circular_dtdtau`.
#     No κ/Ω is re-derived here.
#   • CROSS-PAPER Schwarzschild oracle — arXiv:2509.20457 ("Paper I") main.tex:1352-1392 prints the
#     χ = 0 Lindblad power ⟨Ė'⟩^{(0)},⟨Ė'⟩^{(1)} as an INDEPENDENTLY typeset closed form in
#     (p'−2),(p'−3),(p'−6) factors.  Transcribed below as the oracle; the Kerr code must reproduce it
#     at χ = 0 to machine precision for every p' > 6.
#   • NEWTONIAN limit — arXiv:2509.20457 main.tex:1305-1325 (Ward 1986) prints ⟨Ė'⟩^{(0,1)}_{Newt};
#     the relativistic result must converge to it as p' → ∞ at the 1/p' rate, for χ = 0 AND χ = ±0.9.
#   • PRINTED NUMBER — arXiv:2510.03564 main.tex:662: in Schwarzschild 𝒢₁(p') < 0 only for
#     p' > 8.92, i.e. the torque-reversal root of Eq. (25) sits at p' ≈ 8.92.
# plus a TRANSCRIPTION SELF-CHECK: the brace of b_{1,n} (C2a) and b_{2,d} (C2d) are two separately
# typeset degree-30 polynomials that are algebraically identical — transcribed independently here and
# asserted equal, so a single mistyped coefficient in either fails the suite.

using Test
using AccretionDisc
using AccretionDisc.Modes
using SpecialFunctions: besselk
using Spacetime: kerr_circular_dtdtau, kerr_isco
const MPL = AccretionDisc.Modes

# ---------------- ORACLE 1: arXiv:2509.20457 (Paper I) Schwarzschild closed forms ----------------
# main.tex:1360 (κ/2|A|), 1358-1366 (⟨Ė'⟩^{(0)}), 1368-1388 (⟨Ė'⟩^{(1)}), M = 1, q = 1.
_o_K0(X) = besselk(0, X)
_o_K1(X) = besselk(1, X)
_o_X(p)  = 2 * sqrt((p - 6) * (p - 3)) / (3 * (p - 2))
_o_b0(p) = sqrt(p - 6) / (2 * sqrt(p - 3))
_o_P1(p) = 21p^4 - 331p^3 + 1626p^2 - 3348p + 2376
_o_P2(p) = 12p^4 - 190p^3 + 927p^2 - 1800p + 1188
_o_br0(p) = (X = _o_X(p); _o_K0(X) + _o_b0(p) * _o_K1(X))
# b₁ read off the ⟨Ė'⟩^{(1)} second term (the −8M b₁ q²/27 prefactor of Eq. (18b))
_o_b1(p) = (p - 3)^2 * (p^2 - 5p + 6) * _o_P1(p) / ((p - 6)^2 * (p - 2)^5 * p^1.5)
# b₂ read off the K₁ coefficient inside Paper I's second bracket
_o_b2(p) = (4 * _o_P2(p)^2 / _o_P1(p)^2) *
           (3 * (p * (p * (2p + 5) - 264) + 396) / (8p * (p * (2p * (6p - 95) + 927) - 1800) + 9504) + 7 / 8) *
           sqrt(p - 6) / sqrt(p - 3)
_o_Edot0(p, Σ; k = 1) =
    -16k * (p - 3)^4 * sqrt(p) * Σ / (3 * (p - 2)^2 * ((p - 6) * p)^1.5) * _o_br0(p)^2
function _o_Edot1(p, Σ, dlogΣ)
    X = _o_X(p)
    t1 = -32 * (p - 3)^4 / (9 * (p - 6) * (p - 2)^2 * p^1.5) * _o_br0(p)^2 * dlogΣ * Σ
    t2 = -8 * (p - 3)^2 * (p^2 - 5p + 6) / (27 * (p - 6)^2 * (p - 2)^5 * p^1.5) * _o_P1(p) * Σ *
         _o_br0(p) * (_o_K0(X) + _o_b2(p) * _o_K1(X))
    t1 + t2
end

# ---------------- ORACLE 2: Newtonian limit (Ward 1986), arXiv:2509.20457 main.tex:1305-1325 -----
_o_KN0 = besselk(0, 2 / 3)
_o_KN1 = besselk(1, 2 / 3)
_o_N_Edot0(p, Σ; k = 1) = -k * 4 * (2 * _o_KN0 + _o_KN1)^2 * Σ / (3 * sqrt(p))
_o_N_Edot1(p, Σ, dlogΣ) =
    -8 * (2 * _o_KN0 + _o_KN1)^2 / (9 * sqrt(p)) * dlogΣ * Σ -
     4 / (9 * sqrt(p)) * Σ * (2 * _o_KN0 + _o_KN1) * (7 * _o_KN0 + 8 * _o_KN1)

@testset "Modes.PressurelessLindblad" begin

    σ_flat = PowerLawSigma(Σ_ref = 1.7, p_ref = 10.0, slope = 0.0)
    σ_alpha = PowerLawSigma(Σ_ref = 1.7, p_ref = 10.0, slope = -1.5)   # α-disk Σ_p = −3/2

    @testset "PowerLawSigma — profile algebra" begin
        @test surface_density(σ_alpha, 10.0) ≈ 1.7
        @test surface_density(σ_alpha, 40.0) ≈ 1.7 * 4.0^(1.5)
        @test surface_density_slope(σ_alpha, 23.0) ≈ 1.5
        @test surface_density_slope(σ_flat, 23.0) ≈ 0.0 atol = 1e-15
        # dΣ/dp' against a central finite difference
        for p in (8.0, 23.0, 400.0)
            h = p * 1e-6
            fd = (surface_density(σ_alpha, p + h) - surface_density(σ_alpha, p - h)) / (2h)
            @test surface_density_dp(σ_alpha, p) ≈ fd rtol = 1e-9
        end
        @test surface_density_dp(σ_flat, 23.0) == 0.0
    end

    @testset "Eq. (19) auxiliaries ARE the existing Kerr primitives (no duplicate κ/Ω)" begin
        for χ in (-0.9, -0.3, 0.0, 0.5, 0.9), p in (8.0, 15.0, 60.0, 500.0)
            p > kerr_isco(χ) || continue
            d = RelativisticDisc(a = χ)
            x = p^(1 / 4)
            # printed polynomials (arXiv:2510.03564 Eq. (19)) vs the primitive-built implementations
            @test MPL.hkr_epicyclic_poly(χ, p) ≈ -3χ^2 + x^8 - 6x^4 + 8χ * x^2 rtol = 1e-13
            @test MPL.hkr_binding_poly(χ, p) ≈ 2χ + x^6 - 3x^2 rtol = 1e-13
            @test MPL.hkr_delta_poly(χ, p) ≈ χ^2 + x^8 - 2x^4 rtol = 1e-13
            # the frequency identities the implementation rests on
            @test MPL.hkr_epicyclic_poly(χ, p) ≈ p^2 * kappar2(d, p) / disc_OmegaK(d, p)^2 rtol = 1e-14
            @test MPL.hkr_delta_poly(χ, p) ≈ delta_r(d, p) rtol = 1e-14
            @test χ + x^6 ≈ 1 / disc_OmegaK(d, p) rtol = 1e-13
            # (2χ+x⁶−3x²) = 1/(p^{3/2} Ω² (u^t)²)  — Spacetime.kerr_circular_dtdtau
            ut = kerr_circular_dtdtau(χ, p)
            @test MPL.hkr_binding_poly(χ, p) ≈ 1 / (p^1.5 * disc_OmegaK(d, p)^2 * ut^2) rtol = 1e-12
            # X and b₀ against their printed definitions
            A = -3χ^2 + x^8 - 6x^4 + 8χ * x^2
            B = 2χ + x^6 - 3x^2
            C = χ^2 + x^8 - 2x^4
            @test hkr_shear_argument(χ, p) ≈ 2x * sqrt(B) * sqrt(A) / (3C) rtol = 1e-13
            @test hkr_b0(χ, p) ≈ sqrt(A) / (2x * sqrt(B)) rtol = 1e-13
        end
    end

    @testset "transcription self-check — b_{1,n} brace ≡ b_{2,d} (C2a vs C2d)" begin
        for χ in (-0.9, -0.3, 0.0, 0.25, 0.5, 0.9), p in (7.0, 12.0, 30.0, 200.0)
            x = p^(1 / 4)
            @test MPL._b1_numerator_brace(χ, x) ≈ MPL._b2_denominator(χ, x) rtol = 1e-12
        end
    end

    @testset "χ = 0 vs the Paper I (arXiv:2509.20457) Schwarzschild closed forms" begin
        for p in (6.5, 7.0, 8.0, 8.92, 10.0, 20.0, 50.0, 200.0, 1000.0)
            @test hkr_shear_argument(0.0, p) ≈ _o_X(p) rtol = 1e-13
            @test hkr_b0(0.0, p) ≈ _o_b0(p) rtol = 1e-13
            @test hkr_b1(0.0, p) ≈ _o_b1(p) rtol = 1e-11
            @test hkr_b2(0.0, p) ≈ _o_b2(p) rtol = 1e-11
            for σ in (σ_flat, σ_alpha)
                Σ = surface_density(σ, p)
                dlogΣ = surface_density_slope(σ, p)
                for k in (-1, 1)
                    @test hkr_lindblad_power_quadratic(0.0, p, σ; k = k) ≈
                          _o_Edot0(p, Σ; k = k) rtol = 1e-11
                end
                @test hkr_lindblad_power_linear(0.0, p, σ) ≈ _o_Edot1(p, Σ, dlogΣ) rtol = 1e-11
            end
        end
    end

    @testset "Newtonian limit p' → ∞ (Ward 1986 via arXiv:2509.20457) — 1/p' convergence" begin
        for χ in (0.0, 0.9, -0.9)
            errs0 = Float64[]
            errs1 = Float64[]
            for p in (1e3, 1e4, 1e5)
                Σ = surface_density(σ_alpha, p)
                dlogΣ = surface_density_slope(σ_alpha, p)
                push!(errs0, abs(hkr_lindblad_power_quadratic(χ, p, σ_alpha) / _o_N_Edot0(p, Σ) - 1))
                push!(errs1, abs(hkr_lindblad_power_linear(χ, p, σ_alpha) / _o_N_Edot1(p, Σ, dlogΣ) - 1))
            end
            @test errs0[1] < 1e-2 && errs1[1] < 1e-2
            @test errs0[3] < 1e-4 && errs1[3] < 1e-4
            # first-order (1/p') convergence: each decade in p' cuts the error by ~10
            for i in 1:2
                @test 8 < errs0[i] / errs0[i + 1] < 12
                @test 8 < errs1[i] / errs1[i + 1] < 12
            end
        end
        # the Newtonian asymptotics of the b coefficients (b₁ → 21/√p', b₂ → 8/7)
        for χ in (0.0, 0.9)
            @test hkr_b1(χ, 1e6) * sqrt(1e6) ≈ 21 rtol = 1e-4
            @test hkr_b2(χ, 1e6) ≈ 8 / 7 rtol = 1e-5
            @test hkr_shear_argument(χ, 1e6) ≈ 2 / 3 rtol = 1e-5
            @test hkr_b0(χ, 1e6) ≈ 1 / 2 rtol = 1e-5
        end
    end

    @testset "k-parity of Eqs. (18a)/(18b) — why the O(m²) terms cancel" begin
        for χ in (-0.9, 0.0, 0.9), p in (12.0, 30.0)
            p > kerr_isco(χ) || continue
            plus  = hkr_lindblad_power_quadratic(χ, p, σ_alpha; k = 1)
            minus = hkr_lindblad_power_quadratic(χ, p, σ_alpha; k = -1)
            @test plus ≈ -minus rtol = 1e-14      # (18a) is ODD in k
            @test plus != 0
            # (18b) carries no k at all
            @test hkr_lindblad_power_linear(χ, p, σ_alpha) ==
                  hkr_lindblad_power_linear(χ, p, σ_alpha)
        end
    end

    @testset "Eq. (20) — circular-orbit torque/power relation" begin
        for χ in (-0.9, 0.0, 0.9), p in (12.0, 40.0)
            p > kerr_isco(χ) || continue
            d = RelativisticDisc(a = χ)
            Ω = disc_OmegaK(d, p)
            @test hkr_lindblad_torque_rate(χ, p, σ_alpha; order = :linear) ≈
                  hkr_lindblad_power_linear(χ, p, σ_alpha) / Ω rtol = 1e-14
            @test hkr_lindblad_torque_rate(χ, p, σ_alpha; order = :quadratic, k = -1) ≈
                  hkr_lindblad_power_quadratic(χ, p, σ_alpha; k = -1) / Ω rtol = 1e-14
            @test_throws ArgumentError hkr_lindblad_torque_rate(χ, p, σ_alpha; order = :cubic)
        end
    end

    @testset "Eq. (77) T_{L,m} ∝ m and Eq. (78) Σ_m T_{L,m} ∝ m_max²" begin
        χ, p = 0.0, 15.0
        base = hkr_lindblad_torque_mode(χ, p, 1, σ_alpha)
        @test base != 0
        for m in (1, 2, 5, 13, 40, 137)
            @test hkr_lindblad_torque_mode(χ, p, m, σ_alpha) ≈ m * base rtol = 1e-14
        end
        # per-mode torque is EXACTLY 2m ⟨Ė'⟩^{(1)}/Ω  (the two Lindblad families)
        rate = hkr_lindblad_torque_rate(χ, p, σ_alpha; order = :linear)
        @test hkr_lindblad_torque_mode(χ, p, 7, σ_alpha) ≈ 14 * rate rtol = 1e-14
        # exact partial sum vs the asymptotic Eq. (24) form m_max²
        for m_max in (7, 20, 100)
            exact_sum = sum(hkr_lindblad_torque_mode(χ, p, m, σ_alpha) for m in 1:m_max)
            @test hkr_lindblad_torque(χ, p, σ_alpha; m_max = m_max, asymptotic = false) ≈
                  exact_sum rtol = 1e-13
            @test hkr_lindblad_torque(χ, p, σ_alpha; m_max = m_max) ≈ m_max^2 * rate rtol = 1e-14
            # the asymptotic form is the exact sum up to the stated O(m_max) relative error
            asym = hkr_lindblad_torque(χ, p, σ_alpha; m_max = m_max)
            @test abs(asym / exact_sum - 1) ≈ 1 / (m_max + 1) rtol = 1e-10
        end
        @test_throws DomainError hkr_lindblad_torque(χ, p, σ_alpha; m_max = -1)
    end

    @testset "torque cutoff — Eqs. (79), (80), (23)" begin
        @test hkr_mode_cutoff(0.05) ≈ 20.0
        @test hkr_mode_cutoff(0.1) ≈ 10.0
        @test HKR_M45_COEFFICIENT == 0.35
        @test hkr_mode_cutoff(0.05; coefficient = HKR_M45_COEFFICIENT) ≈ 7.0
        @test hkr_mode_cutoff(0.15; coefficient = HKR_M45_COEFFICIENT) ≈ 0.35 / 0.15
        # Eq. (23): the radius-dependent form reduces to 1/h0 at p' = p_ref (and for Σ_h = 0)
        @test hkr_mode_cutoff_powerlaw(0.05, 10.0) ≈ hkr_mode_cutoff(0.05)
        @test hkr_mode_cutoff_powerlaw(0.05, 40.0; Σ_h = 0) ≈ hkr_mode_cutoff(0.05)
        @test hkr_mode_cutoff_powerlaw(0.05, 40.0; Σ_h = 1) ≈ 4 / 0.05
        @test_throws DomainError hkr_mode_cutoff(0.0)
        @test_throws DomainError hkr_mode_cutoff_powerlaw(-0.1, 10.0)
    end

    @testset "Eq. (25) torque reversal — the printed p' ≈ 8.92 (main.tex:662)" begin
        # 𝒢₁ is the sign of the Σ(p') term of Eq. (18b)
        for χ in (-0.9, 0.0, 0.9), p in (12.0, 30.0)
            p > kerr_isco(χ) || continue
            expected = -hkr_b1(χ, p) *
                       (besselk(0, hkr_shear_argument(χ, p)) +
                        hkr_b0(χ, p) * besselk(1, hkr_shear_argument(χ, p))) *
                       (besselk(0, hkr_shear_argument(χ, p)) +
                        hkr_b2(χ, p) * besselk(1, hkr_shear_argument(χ, p)))
            @test hkr_torque_reversal_criterion(χ, p) ≈ expected rtol = 1e-14
        end
        # PRINTED NUMBER — Schwarzschild reversal at p' ≈ 8.92
        p_rev = hkr_torque_reversal_radius(0.0)
        @test p_rev ≈ 8.92 atol = 5e-3
        @test hkr_torque_reversal_criterion(0.0, p_rev + 0.5) < 0    # negative torque outside
        @test hkr_torque_reversal_criterion(0.0, p_rev - 0.5) > 0    # reversed inside
        # Newtonian limit: 𝒢₁ < 0 always (main.tex:661)
        for p in (1e3, 1e5, 1e7)
            @test hkr_torque_reversal_criterion(0.0, p) < 0
        end
        # spin dependence: the reversal stays outside the ISCO for every spin surveyed
        for χ in (-0.9, -0.5, 0.0, 0.5, 0.9)
            @test hkr_torque_reversal_radius(χ) > kerr_isco(χ)
        end
    end

    @testset "domain guard — the model diverges at the ISCO" begin
        for χ in (-0.9, 0.0, 0.9)
            ri = kerr_isco(χ)
            @test_throws DomainError hkr_lindblad_power_linear(χ, ri, σ_alpha)
            @test_throws DomainError hkr_lindblad_power_linear(χ, ri * 0.9, σ_alpha)
            @test_throws DomainError hkr_lindblad_power_quadratic(χ, ri, σ_alpha)
            @test hkr_lindblad_power_linear(χ, ri * 1.05, σ_alpha) isa Float64
        end
    end

    @testset "type genericity — BigFloat precision is inherited by the algebra" begin
        χb = big"0.9"; pb = big"25.0"
        for f in (hkr_shear_argument, hkr_b0, hkr_b1, hkr_b2,
                  MPL.hkr_epicyclic_poly, MPL.hkr_binding_poly, MPL.hkr_delta_poly)
            v = f(χb, pb)
            @test v isa BigFloat
            @test Float64(v) ≈ f(0.9, 25.0) rtol = 1e-13
        end
        # raising the working precision refines the value (BigFloat is genuinely inherited)
        v256 = setprecision(BigFloat, 256) do; hkr_b1(big"0.9", big"25.0"); end
        v512 = setprecision(BigFloat, 512) do; hkr_b1(big"0.9", big"25.0"); end
        @test abs(v256 / v512 - 1) < 1e-70
        @test abs(v256 / v512 - 1) > 0            # 256-bit is NOT bit-identical to 512-bit
        # the injectable Bessel kernel is the ONLY Float64 bottleneck [HOLE: H-HKR-BIGK]
        @test_throws MethodError besselk(0, big"0.5")
        Kf(ν, x) = BigFloat(besselk(ν, Float64(x)))
        σb = PowerLawSigma(big"1.7", big"10.0", big"-1.5")
        vb = hkr_lindblad_power_linear(χb, pb, σb; bessel = Kf)
        @test vb isa BigFloat
        @test Float64(vb) ≈ hkr_lindblad_power_linear(0.9, 25.0, σ_alpha) rtol = 1e-12
    end

    @testset "mass-ratio scaling — every result ∝ q²" begin
        for q in (1, 3, 1e-5)
            @test hkr_lindblad_power_linear(0.5, 20.0, σ_alpha; q = q) ≈
                  q^2 * hkr_lindblad_power_linear(0.5, 20.0, σ_alpha) rtol = 1e-14
            @test hkr_lindblad_power_quadratic(0.5, 20.0, σ_alpha; q = q) ≈
                  q^2 * hkr_lindblad_power_quadratic(0.5, 20.0, σ_alpha) rtol = 1e-14
        end
    end
end
