# Module 1 — SphericalHarmonics: γ=0 spin-weighted spherical-Y front door + unified m<0 parity.
#
# Grounding:
#  (A) CHARACTERIZATION vs the LIGHT legacy oracle SpinWeightedSpheroidalHarmonics (deps =
#      LinearAlgebra only) — Harmonics.SpinWeightedSphericalHarmonicY (parity=:none) must equal
#      SpinWeightedSpheroidalHarmonics.SpinWeightedSphericalHarmonicY to FULL precision on
#      representative + edge (s,l,m,θ,φ), including the high-l BigFloat path (l ≥ 20).
#  (B) ANALYTIC identities — orthonormality (∫|Y|²dΩ=1), l=|m| closed form (s=0 ∝ sin^l θ),
#      Goldberg/Condon-Shortley reality conj relation, and the :few (-1)^l (m<0) identity.
#  (C) :reality restricted to s=0 (the Punctures/DiscModes scalar case): grounded vs the
#      (-1)^m·conj identity, which for the seed is value-equivalent to :none. The s≠0 spin-flip
#      and full Punctures cross-check remain a [HOLE] (asserted to error for now).

using Test
using Harmonics
# Legacy oracle SpinWeightedSpheroidalHarmonics.SpinWeightedSphericalHarmonicY BAKED into a
# fixture (captured once, same (s,l,m,θ,φ) grid as the loops below). Same characterization at
# the SAME `==`.
include(joinpath(@__DIR__, "fixtures", "spherical_harmonics_ref.jl"))  # FIXSPH_Y

# composite Simpson over [a,b] with N (even) panels — enough for the smooth |Y|² integrand.
function _simpson(f, a, b, N)
    h = (b - a) / N
    s = f(a) + f(b)
    for i in 1:N-1
        s += (isodd(i) ? 4 : 2) * f(a + i * h)
    end
    s * h / 3
end

@testset "SphericalHarmonics" begin

    # ---- (A) characterization vs the legacy oracle, FULL precision ----------------------------
    @testset "characterization vs SWSH oracle (parity=:none)" begin
        cases = [(-2, 2, 2), (-2, 2, -1), (0, 3, 0), (0, 4, -3), (2, 5, 2), (-1, 3, -2),
                 (0, 7, 7), (-2, 8, -8)]
        for (s, l, m) in cases, θ in (0.3, 1.1, 2.7), φ in (0.0, 0.9)
            ref = FIXSPH_Y[(s, l, m, θ, φ)]
            got = SpinWeightedSphericalHarmonicY(s, l, m, θ, φ)
            @test got == ref                       # bit-identical: thin wrapper over same seed (baked)
        end
    end

    @testset "characterization high-l BigFloat path (l ≥ 20)" begin
        for (s, l, m) in [(0, 25, 10), (-2, 30, -5), (2, 22, 22)], θ in (0.4, 2.0)
            ref = FIXSPH_Y[(s, l, m, θ, 0.7)]
            got = SpinWeightedSphericalHarmonicY(s, l, m, θ, 0.7)
            @test got == ref
        end
        # BigFloat element type inherited from core (precision flows through the wrapper).
        setprecision(BigFloat, 256) do
            v = SpinWeightedSphericalHarmonicY(0, 3, 2, big"0.7", big"0.0")
            @test v isa Complex{BigFloat}
            @test real(v) ≈ Float64(real(v)) atol = 1e-12
        end
    end

    # ---- (B) analytic identities --------------------------------------------------------------
    @testset "orthonormality ∫|Y|²dΩ = 1" begin
        for (s, l, m) in [(0, 2, 1), (-2, 3, 2), (2, 4, -1)]
            I = _simpson(θ -> abs2(SpinWeightedSphericalHarmonicY(s, l, m, θ, 0.0)) * sin(θ),
                         0.0, π, 2000)
            @test I * 2π ≈ 1.0 rtol = 1e-6       # φ-integral gives 2π
        end
    end

    @testset "l=|m|, s=0 closed form ∝ sin^l θ" begin
        for l in (2, 3, 5)
            r = SpinWeightedSphericalHarmonicY(0, l, l, 0.8, 0.0) /
                SpinWeightedSphericalHarmonicY(0, l, l, 1.3, 0.0)
            @test real(r) ≈ (sin(0.8) / sin(1.3))^l rtol = 1e-12
            @test imag(r) ≈ 0.0 atol = 1e-12
        end
    end

    @testset "Goldberg/Condon-Shortley conj reality: conj(_sY_lm)=(-1)^(s+m) _-sY_l,-m" begin
        for (s, l, m) in [(-2, 4, 2), (2, 3, -1), (0, 5, 3)], θ in (0.6, 2.2), φ in (0.0, 1.1)
            lhs = conj(SpinWeightedSphericalHarmonicY(s, l, m, θ, φ))
            rhs = (-1)^(s + m) * SpinWeightedSphericalHarmonicY(-s, l, -m, θ, φ)
            @test lhs ≈ rhs rtol = 1e-12
        end
    end

    # ---- (C) unified parity front door --------------------------------------------------------
    @testset "spin_weighted_spherical_Y :none == bare seed" begin
        for (s, l, m) in [(-2, 3, 1), (0, 4, -2)], θ in (0.5, 1.9)
            @test spin_weighted_spherical_Y(s, l, m, θ, 0.3; convention = :none) ==
                  SpinWeightedSphericalHarmonicY(s, l, m, θ, 0.3)
        end
    end

    @testset ":few applies (-1)^l for m<0 (KerrBackwards ylm_few convention)" begin
        for (s, l, m) in [(-2, 3, -1), (-2, 4, -3), (0, 5, -2)], θ in (0.7, 2.0), φ in (0.0, 0.8)
            bare = SpinWeightedSphericalHarmonicY(s, l, m, θ, φ)
            few = spin_weighted_spherical_Y(s, l, m, θ, φ; convention = :few)
            @test few == (-1)^l * bare
        end
        # m ≥ 0: :few is a no-op (== bare).
        for (s, l, m) in [(-2, 3, 2), (0, 4, 0)], θ in (0.7, 2.0)
            @test spin_weighted_spherical_Y(s, l, m, θ, 0.0; convention = :few) ==
                  SpinWeightedSphericalHarmonicY(s, l, m, θ, 0.0)
        end
    end

    @testset "angular_parity sign helper" begin
        @test angular_parity(3, -1; convention = :few) == -1   # (-1)^3
        @test angular_parity(4, -2; convention = :few) == 1    # (-1)^4
        @test angular_parity(3, 2; convention = :few) == 1     # m≥0 ⇒ no flip
        @test angular_parity(5, 4; convention = :none) == 1
        @test angular_parity(4, -3; convention = :reality) == -1  # (-1)^m = (-1)^(-3)
        @test angular_parity(4, -2; convention = :reality) == 1
    end

    @testset ":reality (s=0) == bare seed via (-1)^m conj identity" begin
        for (l, m) in [(3, -2), (4, -1), (5, -3)], θ in (0.6, 1.8), φ in (0.0, 1.2)
            r = spin_weighted_spherical_Y(0, l, m, θ, φ; convention = :reality)
            @test r ≈ SpinWeightedSphericalHarmonicY(0, l, m, θ, φ) rtol = 1e-12
            # explicit reality relation it is built from: _0Y_{l,-|m|}=(-1)^|m| conj(_0Y_{l,|m|})
            @test r ≈ (-1)^abs(m) * conj(SpinWeightedSphericalHarmonicY(0, l, abs(m), θ, φ)) rtol = 1e-12
        end
    end

    @testset ":reality general-s negative-m spin-flip vs SWSH oracle" begin
        # CHARACTERIZATION: the general spin-weighted negative-m reality relation
        # {}_sY_{l,-|m|} = (-1)^{s+m} conj({}_{-s}Y_{l,|m|}) reproduces the oracle bit-identically
        # (same NumericsCore seed) for s = -2..2 across m<0 rows and representative l.
        for s in (-2, -1, 0, 1, 2)
            for l in (max(2, abs(s)), abs(s) + 3, 6)
                for m in -l:-1
                    for θ in (0.4, 1.3, 2.6), φ in (0.0, 0.9)
                        r = spin_weighted_spherical_Y(s, l, m, θ, φ; convention = :reality)
                        @test r == SpinWeightedSphericalHarmonicY(s, l, m, θ, φ)
                        # (a) reality/Goldberg self-consistency: the exact phase it is built from.
                        @test r == (-1)^(s + m) *
                              conj(SpinWeightedSphericalHarmonicY(-s, l, abs(m), θ, φ))
                    end
                end
            end
        end
        # high-l BigFloat half-angle path (l ≥ 20), s ≠ 0, m < 0.
        for (s, l, m) in [(-2, 25, -10), (2, 22, -22), (1, 21, -3)], θ in (0.4, 2.0)
            r = spin_weighted_spherical_Y(s, l, m, θ, 0.7; convention = :reality)
            @test r == FIXSPH_Y[(s, l, m, θ, 0.7)]
        end
        # (b) orthonormality for an s≠0, m<0 :reality harmonic.
        let s = -2, l = 3, m = -2
            I = _simpson(θ -> abs2(spin_weighted_spherical_Y(s, l, m, θ, 0.0; convention = :reality)) *
                              sin(θ), 0.0, π, 2000)
            @test I * 2π ≈ 1.0 rtol = 1e-6
        end
        # BigFloat inheritance through the spin-flip front door.
        setprecision(BigFloat, 256) do
            v = spin_weighted_spherical_Y(-2, 4, -3, big"0.7", big"0.2"; convention = :reality)
            @test v isa Complex{BigFloat}
        end
    end

    @testset ":few parity flag still applies (-1)^l for m<0 (regression)" begin
        for (s, l, m) in [(-2, 3, -1), (1, 4, -2), (0, 5, -3)], θ in (0.7, 2.0)
            bare = SpinWeightedSphericalHarmonicY(s, l, m, θ, 0.4)
            @test spin_weighted_spherical_Y(s, l, m, θ, 0.4; convention = :few) == (-1)^l * bare
        end
    end

    # ---- legacy convention shims (kept green for oracle re-export) -----------------------------
    @testset "ylm / ylm_few / ylms_few shims" begin
        # FEW ylm: bare seed (no flip); 2-arg form fixes s=-2.
        for (l, m) in [(3, 2), (4, -1)], θ in (0.7, 2.1)
            @test ylm(-2, l, m, θ, 0.4) == SpinWeightedSphericalHarmonicY(-2, l, m, θ, 0.4)
            @test ylm(l, m, θ, 0.4) == ylm(-2, l, m, θ, 0.4)
        end
        # KerrBackwards ylm_few: s=-2, m<0 ⇒ (-1)^l flip.
        for (l, m) in [(3, -1), (4, -3), (5, 2)], θ in (0.7, 2.1)
            @test ylm_few(l, m, θ, 0.4) ==
                  spin_weighted_spherical_Y(-2, l, m, θ, 0.4; convention = :few)
        end
        v = ylms_few([3, 4, 5], [-1, -3, 2], 1.0, 0.4)
        @test v == [ylm_few(3, -1, 1.0, 0.4), ylm_few(4, -3, 1.0, 0.4), ylm_few(5, 2, 1.0, 0.4)]
        @test_throws ArgumentError ylms_few([2, 3], [1], 1.0, 0.0)
    end

    @testset "Ylm_real / Ylm_eq re-export (DiscModes real-Y path)" begin
        @test Ylm_real(3, 2, 0.7) isa Float64
        Y, dY, d2Y = Ylm_eq(3, 2, 0.7)
        @test Y == Ylm_real(3, 2, 0.7)
    end
end
