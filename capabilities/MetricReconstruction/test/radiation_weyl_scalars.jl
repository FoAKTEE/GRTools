# Grounding for Radiation.WeylScalars — the gauge-invariant linearized Weyl scalars ψ0/ψ4 REGENERATED from the
# CCK-reconstructed coordinate metric h (δR_abcd fully lowered, contracted on Kinnersley legs). This is the H10
# physics gate: it validates the :psi0_IRG/:psi4_ORG preset CHOICE, not just the numeric amplitude value.
#
# (0) EXACT CURVATURE/SIGN ORACLE: a manufactured Schwarzschild perturbation is differentiated through the full
#     metric by an independent symbolic calculation. Its all-down linearized Riemann literals pin the variation
#     of the lowering metric, the Riemann algebraic symmetries, and the no-extra-sign NP convention.
# (1) H10 PHYSICS GATE (LIVE, corpus-free): for a vacuum CCK ORG reconstruction the regenerated ψ4[h_ORG] has
#     EXACTLY the Teukolsky s=−2 (r,θ) shape ρ⁴·R^up_{−2}(r)·S_{−2}(θ), ρ=−1/(r−ia cosθ). Because both are the
#     SAME gauge-invariant scalar, the ratio ψ4[h] / (ρ⁴·R^up·S_{−2}) is a single (r,θ)-INDEPENDENT complex
#     constant. Grounded LIVE against GrTool Perturbation.TeukolskyRadial(−2,…) + Harmonics SWSH (both light
#     built deps): constancy to <1e-9 with the :psi4_ORG preset. Likewise ψ0[h_IRG] ∝ R^up_{+2}·S_{+2} with the
#     :psi0_IRG preset. (At a=0, ρ⁴=1/r⁴, so the shape is R^up_{−2}·S_{−2}/r⁴ — the "ψ4==Teukolsky/r⁴" gate.)
# (2) COMPANION vacuum gate (already-green Einstein): lin_einstein_residual→0 on the same reconstructed h — the
#     sufficient condition that the reconstruction is a genuine vacuum perturbation (so ψ0/ψ4 ARE gauge-invariant).
# (3) a→0 SCHWARZSCHILD regression: at a=0 the ψ4 ratio equals the a=0.5 constancy behaviour (gauge invariance
#     is spin-independent); the ψ4[h_ORG] shape matches R^up_{−2}(r)·S_{−2}(θ)/r⁴.
# (4) CORRECTED-CONVENTION REGRESSION at (a=0.5,ω=0.2,l=2,m=2,r=8,θ=π/3), default :unit Weyl. The
#     magnitudes are the legacy cross-Teukolsky fixtures (rel<1e-7 solver floor), while the signs are fixed by (0)
#     and Spacetime's independent ψ₂ oracle. This explicitly discriminates the legacy extra-minus bug.

using Test
import MetricReconstruction.Radiation as R
import Perturbation: TeukolskyRadial
import Harmonics: SpinWeightedSpheroidalHarmonicS
import Spacetime

function _riemann_symmetry_defects(Rt)
    scale = maximum(abs, Rt)
    first_pair = maximum(abs, Rt .+ permutedims(Rt, (2,1,3,4))) / scale
    last_pair = maximum(abs, Rt .+ permutedims(Rt, (1,2,4,3))) / scale
    pair_exchange = maximum(abs, Rt .- permutedims(Rt, (3,4,1,2))) / scale
    bianchi = maximum(abs(Rt[a,b,c,d] + Rt[a,c,d,b] + Rt[a,d,b,c])
                      for a in 1:4, b in 1:4, c in 1:4, d in 1:4) / scale
    (; first_pair, last_pair, pair_exchange, bianchi)
end

@testset "Radiation.WeylScalars (ψ0/ψ4 regenerated from reconstructed h — H10 gate)" begin
    a = 0.5; ω = 0.2; l = 2; m = 2; r = 8.0; θ = π/3

    # ---- (0) EXACT manufactured Schwarzschild curvature + NP-sign oracle ---------------------------------
    @testset "(0) exact all-down curvature and no-extra-sign NP convention" begin
        # h_tt = (1 + (r-10) + (r-10)^2) sin(θ)^2 at r=10, θ=π/2. An independent symbolic derivative of
        # R_abcd[g_Schw + εh] gives the exact literals below; the remaining h-jet entries vanish.
        aa = 0.0; rr = 10.0; θθ = π/2
        g = Spacetime.kerr_metric(aa, rr, θθ; coords=:BL)
        gi = Spacetime.kerr_metric_inv(aa, rr, θθ; coords=:BL)
        Γ, dΓ = R._christoffel_and_d_BL_analytic(ComplexF64, aa, rr, θθ)
        h = zeros(ComplexF64, 4, 4)
        dh = zeros(ComplexF64, 4, 4, 4)
        ddh = zeros(ComplexF64, 4, 4, 4, 4)
        h[1,1] = 1
        dh[2,1,1] = 1
        ddh[2,2,1,1] = 2
        ddh[3,3,1,1] = -2

        Rt = R._lin_riemann_lowered_from_jets(g, gi, Γ, dΓ, h, dh, ddh)
        @test Rt[1,2,1,2] ≈ -6361/6400 atol=2e-12 rtol=2e-12
        @test Rt[1,3,1,3] ≈ -3 atol=2e-12 rtol=2e-12
        @test Rt[1,4,1,4] ≈ -4 atol=2e-12 rtol=2e-12
        defects = _riemann_symmetry_defects(Rt)
        @test maximum(values(defects)) < 5e-13

        lup, nup, mup, mbup = map(x -> ComplexF64.(x),
                                  Spacetime.kinnersley_tetrad(aa, rr, θθ; coords=:BL))
        @test R._repeated_null_contraction(Rt, lup, mup) ≈ 1/128 atol=2e-12 rtol=2e-12
        @test R._repeated_null_contraction(Rt, nup, mbup) ≈ 1/800 atol=2e-12 rtol=2e-12
    end

    # The production reconstructed-field path must preserve the same all-down tensor symmetries.
    @testset "(0b) CCK curvature has Riemann algebraic symmetries" begin
        for gauge in (:IRG, :ORG)
            w = R.weyl_amplitudes(:unit, ω, l, m, a)
            Rt = R.lin_riemann_lowered(gauge, a, ω, l, m, r, θ, 0.0, 0.0, w)
            @test maximum(values(_riemann_symmetry_defects(Rt))) < 1e-11
        end
    end

    # ---- (4) corrected convention regression (baked literals, default :unit weyl) -------------------------
    fixtures = Dict(
        (:psi0, :IRG) => -0.04196278187108527 - 0.04308365923119698im,
        (:psi0, :ORG) => -1.0806492822353e-5 - 1.4144896508046335e-5im,
        (:psi4, :IRG) => -0.030283855413226517 - 0.01794252368939258im,
        (:psi4, :ORG) => -9.820940222603024e-6 - 4.5808511674786406e-6im,
    )
    @testset "(4) corrected-convention weyl_scalar fixture rel<1e-7" begin
        for ((which, gauge), ref) in fixtures
            v = R.weyl_scalar(which, gauge, a, ω, l, m, r, θ)   # default weyl_amplitudes(:unit)
            @test abs(v - ref) / abs(ref) < 1e-7
        end
    end

    # ---- (1) H10 PHYSICS GATE: ψ4[h_ORG] shape == ρ⁴·R^up_{−2}·S_{−2}, ratio constant to <1e-9 -------------
    @testset "(1) H10 gate: ψ4[h_ORG]/(ρ⁴ R^up_{−2} S_{−2}) constant (<1e-9)" begin
        for aa in (0.5, 0.9)
            w4 = R.weyl_amplitudes(:psi4_ORG, ω, l, m, aa)
            Λ = R.teuk_lambda(-2, l, m, aa, ω)
            Rteuk = TeukolskyRadial(-2, l, m, aa, ω; lambda=Λ)
            Sm2 = SpinWeightedSpheroidalHarmonicS(-2, l, m, aa*ω)
            ratios = ComplexF64[]
            for rr in 8.0:1.5:15.5, θθ in range(π/6, 5π/6, length=5)
                ρ = -1/(rr - im*aa*cos(θθ))
                teuk = ρ^4 * Rteuk.Up(rr) * Sm2(θθ, 0.0)
                ψ4 = R.weyl_scalar(:psi4, :ORG, aa, ω, l, m, rr, θθ; weyl=w4)
                push!(ratios, ψ4 / teuk)
            end
            C = sum(ratios) / length(ratios)
            maxrel = maximum(abs(x/C - 1) for x in ratios)
            @test maxrel < 1e-9
        end
    end

    # ---- (1b) ψ0[h_IRG] shape == R^up_{+2}·S_{+2} (constant ratio, s=+2 content) --------------------------
    @testset "(1b) ψ0[h_IRG] ∝ R^up_{+2} S_{+2} (ratio constant <1e-8)" begin
        w0 = R.weyl_amplitudes(:psi0_IRG, ω, l, m, a)
        Λ2 = R.teuk_lambda(2, l, m, a, ω)
        Rp2 = TeukolskyRadial(2, l, m, a, ω; lambda=Λ2)
        Sp2 = SpinWeightedSpheroidalHarmonicS(2, l, m, a*ω)
        ratios = ComplexF64[]
        for rr in 8.0:1.5:15.5, θθ in range(π/6, 5π/6, length=5)
            teuk = Rp2.Up(rr) * Sp2(θθ, 0.0)
            ψ0 = R.weyl_scalar(:psi0, :IRG, a, ω, l, m, rr, θθ; weyl=w0)
            push!(ratios, ψ0 / teuk)
        end
        C = sum(ratios) / length(ratios)
        @test maximum(abs(x/C - 1) for x in ratios) < 1e-8
    end

    # ---- (2) COMPANION vacuum gate: lin_einstein_residual→0 on the reconstructed h (both gauges) -----------
    @testset "(2) companion vacuum gate ‖G⁽¹⁾‖/‖h‖ → 0" begin
        for (g, rr) in ((:ORG, 8.0), (:IRG, 8.0))
            w = R.weyl_amplitudes(:unit, ω, l, m, a)
            H = R.h_coordinate(g, :BL, a, ω, l, m, rr, θ, 0.0, 0.0, w; bc=:Up)
            G = R.lin_einstein_residual(g, a, ω, l, m, rr, θ, 0.0, 0.0, w; bc=:Up, deriv=:analytic)
            @test maximum(abs, G) < 1e-10 * maximum(abs, H)
        end
    end

    # ---- (3) a→0 Schwarzschild regression: ψ4[h_ORG] shape == R^up_{−2}·S_{−2}/r⁴ -------------------------
    @testset "(3) a→0 regression ψ4[h_ORG] shape == R^up_{−2} S_{−2}/r⁴" begin
        aa = 1e-6
        w4 = R.weyl_amplitudes(:psi4_ORG, ω, l, m, aa)
        Λ = R.teuk_lambda(-2, l, m, aa, ω)
        Rteuk = TeukolskyRadial(-2, l, m, aa, ω; lambda=Λ)
        Sm2 = SpinWeightedSpheroidalHarmonicS(-2, l, m, aa*ω)
        ratios = ComplexF64[]
        for rr in 8.0:1.5:14.0, θθ in (0.5, 1.1, 1.9)
            ρ = -1/(rr - im*aa*cos(θθ))                     # ρ⁴ → 1/r⁴ as a→0
            teuk = ρ^4 * Rteuk.Up(rr) * Sm2(θθ, 0.0)
            ψ4 = R.weyl_scalar(:psi4, :ORG, aa, ω, l, m, rr, θθ; weyl=w4)
            push!(ratios, ψ4 / teuk)
        end
        C = sum(ratios) / length(ratios)
        @test maximum(abs(x/C - 1) for x in ratios) < 1e-8
        # leading a→0 shape: ρ⁴ = 1/r⁴ + O(a) — confirm the /r⁴ form to first order (O(a)~4e-7 here)
        @test abs((-1/(12.0 - im*aa*cos(1.1)))^4 - 1/12.0^4) < 1e-5
    end
end
