# Characterization of the normalized associated-Legendre + spin-weighted spherical-harmonic angular SEED.
#
# Grounding:
#  - spin_weighted_Y / _y_dtheta: oracle = SpinWeightedSpheroidalHarmonics, BAKED into
#    fixtures/numericscore_swsh_ref.jl (Phase-7-B, captured once while SWSH was still resolvable).
#    SWSH.SpinWeightedSphericalHarmonicY. Matched at full Float64 precision on representative + edge
#    (s,l,m), INCLUDING the l>=_LHIGH_SWSHY(=20) BigFloat path (Int64-overflow regime). The
#    θ-derivative is checked against the baked oracle's _y_dtheta AND vs central finite-diff.
#  - assoc_legendre_P / legendre_P_thetaderiv: READ-literal closed forms reproduced from
#    H1Lorenz/src/math_funcs.jl, cross-checked vs Abramowitz&Stegun analytic forms:
#       P_1^0=cosθ, P_1^1=-sinθ (Condon-Shortley), P_2^0=(3cos²θ-1)/2, P_2^2=3sin²θ.
#  - Ylm_real / Ylm_eq: analytic spherical-harmonic anchors (DiscModes is heavy → no live dep):
#       Y_00=1/√(4π), Y_10=√(3/4π)cosθ, Y_20=√(5/16π)(3cos²θ-1); orthonormality ∫|Y|²=1.
#  - pl_deriv: analytic closed form (Punctures heavy → no live dep): P_l'(1)=l(l+1)/2, P_l(1)=1,
#       P_l'(-1)=(-1)^(l+1) l(l+1)/2.
#  - Analytic: ₀Y₀₀=1/√(4π); Condon-Shortley phase; BigFloat genericity to 2^-200.

using NumericsCore
using NumericsCore: _fact, _legendre_Plm, assoc_legendre_P, legendre_P_thetaderiv,
                    assoc_legendre, Ylm_real, Ylm_eq, pl_deriv, spin_weighted_Y,
                    _y_dtheta, _LHIGH_SWSHY
# Baked from SpinWeightedSpheroidalHarmonics (captured once, Phase-7-B) — no live legacy dep.
include("fixtures/numericscore_swsh_ref.jl")

@testset "legendre" begin

    @testset "_fact (BigInt-safe factorial in T)" begin
        @test _fact(Float64, 0) == 1.0
        @test _fact(Float64, 1) == 1.0
        @test _fact(Float64, 5) == 120.0
        @test _fact(Int, 6) == 720
        # overflow-free at high l in BigFloat
        @test _fact(BigFloat, 25) ≈ BigFloat(factorial(big(25)))
    end

    @testset "_legendre_Plm closed forms (A&S, Condon-Shortley)" begin
        for x in (-0.7, -0.2, 0.0, 0.3, 0.85)
            s = sqrt(1 - x^2)
            @test _legendre_Plm(0, 0, x) ≈ 1.0
            @test _legendre_Plm(1, 0, x) ≈ x
            @test _legendre_Plm(1, 1, x) ≈ -s              # Condon-Shortley -(1-x²)^½
            @test _legendre_Plm(2, 0, x) ≈ (3x^2 - 1) / 2
            @test _legendre_Plm(2, 1, x) ≈ -3 * x * s
            @test _legendre_Plm(2, 2, x) ≈ 3 * (1 - x^2)
        end
        @test _legendre_Plm(3, 4, 0.3) == 0.0              # m>l ⇒ 0
        # assoc_legendre (DiscModes name) is the SAME consolidated primitive
        @test assoc_legendre(2, 1, 0.4) == _legendre_Plm(2, 1, 0.4)
        @test assoc_legendre(5, 3, -0.6) == _legendre_Plm(5, 3, -0.6)
    end

    @testset "assoc_legendre_P at θ (GSL Plm convention + m<0 sign + parity)" begin
        for θ in (0.4, 1.1, 2.3)
            x = cos(θ)
            @test assoc_legendre_P(1, 0, θ) ≈ x
            @test assoc_legendre_P(1, 1, θ) ≈ -sin(θ)
            @test assoc_legendre_P(2, 0, θ) ≈ (3x^2 - 1) / 2
            # m<0: (-1)^|m| times the |m| value
            @test assoc_legendre_P(2, -1, θ) ≈ -assoc_legendre_P(2, 1, θ)
            @test assoc_legendre_P(2, -2, θ) ≈ assoc_legendre_P(2, 2, θ)
        end
        # θ=π/2 parity shortcut: l+m odd ⇒ exactly 0
        @test assoc_legendre_P(2, 1, π / 2) == 0.0
        @test assoc_legendre_P(3, 0, π / 2) == 0.0
        # normalized=true: GSL sphPlm. ₀Y normalization N_l^|m| P_l^|m|. Check N_1^0 P_1^0 at θ=0.4.
        N10 = sqrt(3 / (4π))
        @test assoc_legendre_P(1, 0, 0.4; normalized=true) ≈ N10 * cos(0.4)
        N11 = sqrt(3 / (8π))
        @test assoc_legendre_P(1, 1, 0.4; normalized=true) ≈ N11 * (-sin(0.4))
    end

    @testset "legendre_P_thetaderiv (legacy -dP/d(cosθ); sinθ × it = true ∂_θP) vs FD" begin
        # The legacy math_funcs.c closed form divides by (1-cos²θ) not sinθ, so it returns
        # -dP_l^m/d(cosθ) = (1/sinθ)·dP_l^m/dθ.  Ground that exact convention vs central FD of
        # the colatitude derivative: sinθ · legendre_P_thetaderiv ≈ ∂_θ assoc_legendre_P.
        h = 1e-6
        for (l, m) in ((1, 0), (2, 0), (2, 1), (3, 2), (4, 1))
            for θ in (0.5, 1.3, 2.0)
                fd = (assoc_legendre_P(l, m, θ + h) - assoc_legendre_P(l, m, θ - h)) / (2h)
                @test sin(θ) * legendre_P_thetaderiv(l, m, θ) ≈ fd rtol=1e-5
                fdn = (assoc_legendre_P(l, m, θ + h; normalized=true) -
                       assoc_legendre_P(l, m, θ - h; normalized=true)) / (2h)
                @test sin(θ) * legendre_P_thetaderiv(l, m, θ; normalized=true) ≈ fdn rtol=1e-5
            end
        end
        # exact legacy closed-form literal at l=1,m=0: -(-cosθ·cosθ + 1·1)/(1-cos²θ) = -1
        @test legendre_P_thetaderiv(1, 0, 0.7) ≈ -1.0
    end

    @testset "Ylm_real / Ylm_eq analytic spherical-harmonic anchors" begin
        for θ in (0.3, 1.0, 1.9, 2.7)
            @test Ylm_real(0, 0, θ) ≈ 1 / sqrt(4π)
            @test Ylm_real(1, 0, θ) ≈ sqrt(3 / (4π)) * cos(θ)
            @test Ylm_real(2, 0, θ) ≈ sqrt(5 / (16π)) * (3cos(θ)^2 - 1)
            # |m| invariance of Ylm_real (uses |m|, no sign flip — DiscModes convention)
            @test Ylm_real(3, -2, θ) == Ylm_real(3, 2, θ)
        end
        # Ylm_eq returns (Y, dY, d2Y); check value + finite-diff derivs
        h = 1e-5
        Y0, dY, d2Y = Ylm_eq(2, 1, 1.0)
        @test Y0 ≈ Ylm_real(2, 1, 1.0)
        @test dY ≈ (Ylm_real(2, 1, 1.0 + h) - Ylm_real(2, 1, 1.0 - h)) / (2h) rtol=1e-4
        @test d2Y ≈ (Ylm_real(2, 1, 1.0 + h) - 2Ylm_real(2, 1, 1.0) +
                     Ylm_real(2, 1, 1.0 - h)) / h^2 rtol=1e-3
    end

    @testset "Ylm orthonormality (Gauss-Legendre over θ)" begin
        # ∫|Y_lm|² dΩ = 1.  φ-integral = 2π (|e^{imφ}|=1), so 2π ∫_0^π |Y_lm|² sinθ dθ = 1.
        nodes, weights = let n = 60
            # simple high-res midpoint over θ∈(0,π); enough for low-l polynomials
            xs = [(k - 0.5) * π / n for k in 1:n]
            ws = fill(π / n, n)
            xs, ws
        end
        for (l, m) in ((0, 0), (1, 0), (2, 1), (3, 2))
            integ = sum(weights[i] * Ylm_real(l, m, nodes[i])^2 * sin(nodes[i]) for i in 1:length(nodes))
            @test 2π * integ ≈ 1.0 rtol=2e-3
        end
    end

    @testset "pl_deriv analytic (Punctures, P_l derivatives at ±1)" begin
        for l in 0:6
            @test pl_deriv(Float64, l, 0, 1) ≈ 1.0           # P_l(1)=1
            @test pl_deriv(Float64, l, 0, -1) ≈ (-1.0)^l     # P_l(-1)=(-1)^l
            @test pl_deriv(Float64, l, 1, 1) ≈ l * (l + 1) / 2          # P_l'(1)
            @test pl_deriv(Float64, l, 1, -1) ≈ (-1.0)^(l + 1) * l * (l + 1) / 2
        end
        # BigFloat genericity
        @test pl_deriv(BigFloat, 10, 1, 1) ≈ BigFloat(10 * 11) / 2
    end

    @testset "spin_weighted_Y vs SWSH oracle (Float64 path, l<20) (baked from SWSH, captured once)" begin
        reps = [(0, 0, 0), (-2, 2, 0), (-2, 2, 2), (-2, 2, -1), (0, 3, 2),
                (-1, 4, -3), (2, 5, 4), (-2, 10, 7), (0, 19, 11)]
        for (s, l, m) in reps, θ in (0.37, 1.05, 2.41), φ in (0.0, 0.83)
            ref = REF_Y[(float(s), float(l), float(m), θ, φ)]
            @test spin_weighted_Y(s, l, m, θ, φ) ≈ ref rtol=1e-13 atol=1e-300
        end
        # analytic anchor ₀Y₀₀ = 1/√(4π)
        @test spin_weighted_Y(0, 0, 0, 1.0, 0.5) ≈ complex(1 / sqrt(4π))
    end

    @testset "spin_weighted_Y high-l BigFloat path (l>=_LHIGH_SWSHY) vs SWSH oracle (baked from SWSH, captured once)" begin
        @test _LHIGH_SWSHY == 20
        for (s, l, m) in [(-2, 22, 15), (0, 25, -10), (-2, 30, 30)], θ in (0.6, 1.7)
            ref = REF_Y_highl[(float(s), float(l), float(m), θ, 0.4)]
            @test spin_weighted_Y(s, l, m, θ, 0.4) ≈ ref rtol=1e-12
        end
    end

    @testset "_y_dtheta vs SWSH oracle + finite difference (baked from SWSH, captured once)" begin
        h = 1e-6
        for (s, l, m) in [(-2, 2, 1), (0, 3, 2), (-2, 5, -3), (-2, 22, 14)], θ in (0.7, 1.6)
            ref = REF_ydtheta[(float(s), float(l), float(m), θ)]
            @test _y_dtheta(s, l, m, θ, 1) ≈ ref rtol=1e-12
            # cross-check vs FD of the real θ-part (l<20 only — BigFloat FD too costly here)
            if l < _LHIGH_SWSHY
                yp = real(spin_weighted_Y(s, l, m, θ + h, 0.0))
                ym = real(spin_weighted_Y(s, l, m, θ - h, 0.0))
                @test _y_dtheta(s, l, m, θ, 1) ≈ (yp - ym) / (2h) rtol=1e-5
            end
        end
    end

    @testset "BigFloat precision inheritance (2^-200)" begin
        setprecision(BigFloat, 256) do
            θ = BigFloat(1) / 3
            x = cos(θ)
            @test abs(_legendre_Plm(2, 0, x) - (3x^2 - 1) / 2) < big(2.0)^-200
            @test assoc_legendre_P(1, 1, θ) isa BigFloat
            @test abs(assoc_legendre_P(1, 1, θ) - (-sin(θ))) < big(2.0)^-200
            @test abs(Ylm_real(1, 0, θ) - sqrt(BigFloat(3) / (4 * BigFloat(π))) * cos(θ)) < big(2.0)^-200
            y = spin_weighted_Y(0, 0, 0, θ, BigFloat(0))
            @test abs(real(y) - 1 / sqrt(4 * BigFloat(π))) < big(2.0)^-200
        end
    end
end
