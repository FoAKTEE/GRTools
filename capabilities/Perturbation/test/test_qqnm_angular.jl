# test_qqnm_angular.jl — Q4-H Deliverable A gates (the general-a angular reduction).
#   G-Q4H-3  a→0 continuity onto the certified a=0 3-j route (the DECISIVE Deliverable-A gate).
#   + the L†_n ladder identities (a=0 Schwarzschild reduction, Paper I:293,749).
#   + the paper leading-order angular projections ⟨(L₂†S₂₂)²|S₄₄⟩, ⟨S₂₂L₁†L₂†S₂₂|S₄₄⟩ (Paper I:1088).
#   + general-a source-coefficient continuity (qqnm_source_coeffs GL joint projection → a=0 3-j).
using Test
using Printf
using Perturbation
using Harmonics: SpinWeightedSphericalHarmonicY, SpinWeightedSpheroidalHarmonicS, swsh_derivative
import Perturbation: _qqnm_proj_a0_analytic

@testset "Q4-H Deliverable A — general-a angular reduction" begin
    ωL = ComplexF64(0.37367168441804 - 0.08896231568894im)   # 220 (a=0)

    # ── L†_n ladders as first-order θ-operators == the a=0 Schwarzschild reduction (Paper I:293,749) ──
    @testset "L†_n ladder reduction (a=0)" begin
        for (lL, mL) in ((2, 2), (3, 2))
            for θ in (0.4, 1.0, 1.9, 2.6)
                # reconstruct L₂†S, L₁†L₂†S from the structures (𝒜₁=[L₂†S]², 𝒜₂=S·L₁†L₂†S)
                (A1, A2, A3, A4) = qqnm_angular_structures(lL, mL, 0.0 + 0.0im, 0.0, θ)
                S0 = real(SpinWeightedSphericalHarmonicY(-2, lL, mL, θ, 0.0))
                y1 = real(SpinWeightedSphericalHarmonicY(-1, lL, mL, θ, 0.0))
                y0 = real(SpinWeightedSphericalHarmonicY(0, lL, mL, θ, 0.0))
                L2S_expected = -sqrt((lL + 2) * (lL - 1)) * y1
                L1L2S_expected = sqrt((lL + 2) * (lL + 1) * lL * (lL - 1)) * y0
                @test abs(A1 - L2S_expected^2) < 1e-11           # 𝒜₁ = [L₂†S]²
                @test abs(A2 - S0 * L1L2S_expected) < 1e-11      # 𝒜₂ = S·L₁†L₂†S
                @test A3 == 0 && A4 == 0                         # explicit-a structures vanish at a=0
            end
        end
    end

    # ── G-Q4H-3 the a→0 continuity gate (the DECISIVE Deliverable-A gate) ─────────────────────────────
    @testset "G-Q4H-3 a→0 continuity onto the certified 3-j route" begin
        for (lL, lQ, mQ) in ((2, 4, 4), (2, 6, 4), (3, 6, 6))
            (P_A, P_B) = _qqnm_proj_a0_analytic(Float64, lL, lL, lQ, mQ)
            cA1 = (lL + 2) * (lL - 1) * P_B                       # Π₁ a=0 target
            cA2 = sqrt((lL + 2) * (lL + 1) * lL * (lL - 1)) * P_A # Π₂ a=0 target
            # (i) a=0 EXACT: the complex-γ machinery reduces to the 3-j route to the value-tier floor
            Π = qqnm_angular_projections(lL, lL, lQ, mQ, ωL, 0.0; ny=96)
            @test abs(Π[1] - cA1) ≤ 1e-10 * max(abs(cA1), 1.0)
            @test abs(Π[2] - cA2) ≤ 1e-10 * max(abs(cA2), 1.0)
            @test abs(Π[3]) ≤ 1e-12 && abs(Π[4]) ≤ 1e-12          # 𝒜₃,𝒜₄ ∝ a vanish at a=0
        end
        # (ii) finite-a drift is SMOOTH physical O(a): (Π₁(a)−cA1)/a → a finite Richardson constant
        (P_A, P_B) = _qqnm_proj_a0_analytic(Float64, 2, 2, 4, 4)
        cA1 = 4 * 1 * P_B
        drift(a) = real((qqnm_angular_projections(2, 2, 4, 4, ωL, a; ny=96)[1] - cA1) / (a * cA1))
        d2, d3, d4, d6 = drift(1e-2), drift(1e-3), drift(1e-4), drift(1e-6)
        @test abs(d6 - d4) < 1e-3 && abs(d4 - d3) < 1e-3         # Richardson-convergent (finite dΠ/da)
        # (iii) Π₃ = O(a), Π₄ = O(a²) — the explicit-a orders
        Πa = qqnm_angular_projections(2, 2, 4, 4, ωL, 1e-3; ny=96)
        @test abs(Πa[3]) < 5e-3 && abs(Πa[3]) > 1e-8            # O(a)
        @test abs(Πa[4]) < 5e-6                                  # O(a²)
    end

    # ── the paper leading-order angular projections (Paper I:1088 commented block; a=0 3-j) ──────────
    @testset "leading-order projections ⟨(L₂†S₂₂)²|S₄₄⟩, ⟨S₂₂L₁†L₂†S₂₂|S₄₄⟩" begin
        Π = qqnm_angular_projections(2, 2, 4, 4, ωL, 0.0; ny=96)
        # a=0: ⟨(L₂†S₂₂)²|S₄₄⟩ = (l+2)(l−1)·[3jB] = 4·(5/(3√7π)); ⟨S₂₂L₁†L₂†S₂₂|S₄₄⟩ = √((l+2)(l+1)l(l−1))·[3jA]
        @test abs(Π[1] - 4 * (5 / (3 * sqrt(7π)))) < 1e-10
        @test abs(Π[2] - sqrt(4 * 3 * 2 * 1) * (5 / (2 * sqrt(42π)))) < 1e-10
    end

    # ── general-a source-coefficient continuity (GL joint Γ(θ)-coupled projection → a=0 3-j) ─────────
    @testset "qqnm_source_coeffs GL a→0 continuity (W00,W01,W11)" begin
        p = (l=2, m=2, n=0, ω=ωL)
        rs = collect((4.0, 7.0, 12.0))
        co0 = qqnm_source_coeffs(p, p, :pp, 4; a=0.0, r=rs, T=ComplexF64)          # certified 3-j
        coa = qqnm_source_coeffs(p, p, :pp, 4; a=1e-6, r=rs, T=ComplexF64)         # GL joint, a→0
        for fld in (:W00, :W01, :W11)
            v0 = getfield(co0, fld); va = getfield(coa, fld)
            @test maximum(abs.(va .- v0) ./ max.(abs.(v0), 1.0)) < 1e-4            # O(a) drift bound
        end
    end

    # ── visibility: the a→0 continuity numbers (G-Q4H-3 pasted) ──────────────────────────────────────
    println("─ G-Q4H-3 a→0 continuity (lL=2→lQ=4, mQ=4) ─")
    (P_A, P_B) = _qqnm_proj_a0_analytic(Float64, 2, 2, 4, 4)
    cA1 = 4 * P_B; cA2 = sqrt(4 * 3 * 2 * 1) * P_A
    for a in (0.0, 1e-6, 1e-4, 1e-2)
        Π = qqnm_angular_projections(2, 2, 4, 4, ωL, a; ny=96)
        @printf("  a=%.0e  relΔΠ₁=%.3e  relΔΠ₂=%.3e  |Π₃|=%.3e  |Π₄|=%.3e\n",
                a, abs(Π[1] - cA1) / abs(cA1), abs(Π[2] - cA2) / abs(cA2), abs(Π[3]), abs(Π[4]))
    end
    println("  [NOTE] the a=0-EXACT reduction meets ≤1e-10; finite-a residual is the genuine physical")
    println("  O(a) spheroidal drift (Richardson const dΠ₁/da/cA1 ≈ −0.22), NOT numerical error.")
end
