# Q1a — COMPLEX-γ spheroidal machinery (the module header's declared [FUTURE], now BUILT).
# γ = aω with a COMPLEX QNM frequency ω (the QQNM angular basis at aω_L and 2aω_L). The eigenvalue
# extends as a complex-SYMMETRIC (not Hermitian) problem: the spectral tier uses the dense
# nonsymmetric geev with eigenvector-l-component branch tracking; the Leaver CF tier is
# {T}/Complex-generic and BigFloat-inheriting. The value tier (SphericalExpansion) carries complex
# coeffs (Float64-capped by design §4.2).
#
# NO legacy oracle: these are ANALYTIC / self-consistency gates (real-γ regression is covered by the
# BAKED fixtures in the sibling test files — the entire existing suite must stay bit-unchanged).
# Tolerances are MEASURED (probe_complex.jl), stated with headroom, never tuned.

using Test
using Harmonics

@testset "Spheroidal COMPLEX-γ (Q1a)" begin

    # ---- (a) REAL-γ regression: complex path at zero imaginary reproduces the real path ----------
    # (the full bit-exact regression is the existing fixture suites; here: continuity of the new
    #  complex spectral/Leaver branch onto the real answer.)
    @testset "(a) complex(x+0im) == real(x) continuity" begin
        for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 4, 2), (2, 4, -3)], x in (0.3, 0.7, 1.2)
            re = SpinWeightedSpheroidalEigenvalue(s, l, m, x; method=:SphericalExpansion)
            cx = SpinWeightedSpheroidalEigenvalue(s, l, m, complex(x); method=:SphericalExpansion)
            @test imag(cx) == 0.0 || abs(imag(cx)) < 1e-12 * abs(cx)
            @test real(cx) ≈ re rtol = 1e-12           # geev vs Symmetric eigen: same eigenvalue
            reL = SpinWeightedSpheroidalEigenvalue(s, l, m, x; method=:Leaver)
            cxL = SpinWeightedSpheroidalEigenvalue(s, l, m, complex(x); method=:Leaver)
            @test real(cxL) ≈ reL rtol = 1e-12
        end
    end

    # ---- (b) Leaver CF residual == 0 (machine floor) at BF-256, complex γ, all quadrants ---------
    # cover aω_L (|γ|≈0.3–0.5) AND 2aω_L (|γ|≈0.95–1.2) magnitudes (χ up to 0.95); arg spans all
    # four quadrants (incl. the physical 4th-quadrant QNM γ = a(ω_R − i|ω_I|)).
    @testset "(b) Leaver CF residual ≈ machine floor @ BF-256, complex γ" begin
        setprecision(BigFloat, 256) do
            for (s, l, m) in [(-2, 2, 2), (-2, 3, 3), (-2, 2, 0), (2, 2, 2), (0, 2, 1), (-2, 4, -2)]
                for mag in (big"0.3", big"0.5", big"0.95", big"1.2")   # aω_L and 2aω_L magnitudes
                    for arg in (big"0.4", big"1.9", big"-0.6", big"2.7", big"-2.4")  # all quadrants
                        γ = mag * exp(im * arg)
                        λ, res = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; return_residual=true)
                        @test res isa BigFloat
                        @test res >= 0
                        @test res < eps(BigFloat) * 256      # observed max ≈ 9.2×eps @256 (1.6e-76)
                        @test λ isa Complex{BigFloat}
                    end
                end
            end
        end
    end

    # ---- (c) conjugation symmetry λ(s,l,m,γ̄) == conj(λ(s,l,m,γ)) --------------------------------
    @testset "(c) λ(conj γ) == conj(λ(γ))" begin
        for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 3, 0), (2, 4, 2)]
            for γ in (0.5 + 0.3im, 1.0 - 0.7im, 0.2 + 0.9im, 1.2 + 0.0im)     # spectral (F64)
                a = SpinWeightedSpheroidalEigenvalue(s, l, m, γ)
                b = SpinWeightedSpheroidalEigenvalue(s, l, m, conj(γ))
                @test abs(b - conj(a)) <= 1e-12 * abs(a)     # observed 0
            end
        end
        setprecision(BigFloat, 256) do
            for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 3, 0)]
                for γ in (big"0.5" + big"0.3" * im, big"1.0" - big"0.7" * im, big"1.2" + big"0.0" * im)
                    a = SpinWeightedSpheroidalEigenvalue(s, l, m, γ)          # Leaver (BF256)
                    b = SpinWeightedSpheroidalEigenvalue(s, l, m, conj(γ))
                    @test abs(b - conj(a)) <= eps(BigFloat) * 256 * abs(a)    # observed 0
                end
            end
        end
    end

    # ---- (d) small-γ consistency vs SpheroidalEigenvalueSeries (|γ|<=0.1, complex) ---------------
    @testset "(d) small-γ series == eigenvalue (complex, |γ|<=0.1)" begin
        for (s, l, m) in [(-2, 2, 0), (-2, 3, 2), (0, 2, 1), (2, 4, 4)]
            c8 = SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=8)  # accuracy (real coeffs)
            c6 = SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=6)  # the :SmallGamma default order
            for γ in (0.05 + 0.02im, 0.08 - 0.05im, 0.03 + 0.09im, 0.1 + 0.0im)
                ser = evalpoly(γ, c8)
                spec = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
                leav = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
                small = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SmallGamma)  # default order 6
                @test small == evalpoly(γ, c6)            # :SmallGamma route wiring is EXACT (same order)
                @test ser ≈ spec rtol = 1e-12             # observed ≤2.5e-15
                @test ser ≈ leav rtol = 1e-12
            end
        end
    end

    # ---- (e) spherical limit γ→0 == (l-s)(l+s+1) exact -------------------------------------------
    @testset "(e) γ→0 == (l-s)(l+s+1) exact" begin
        for (s, l, m) in [(s, l, m) for s in (-2, 0, 2)
                          for l in (max(abs(s), 2):max(abs(s), 2)+2) for m in (-l, 0, 1, l) if abs(m) <= l]
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0 + 0.0im) == (l - s) * (l + s + 1)
            @test SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, 0.0 + 0.0im) ==
                  -2 * m * s^2 / (l + l^2) - 2 * m
        end
    end

    # ---- (f) eltype / BF-256 cleanliness --------------------------------------------------------
    @testset "(f) eltype + BF-256 inheritance" begin
        @test SpinWeightedSpheroidalEigenvalue(-2, 2, 2, 0.5 + 0.3im) isa ComplexF64
        @test SpinWeightedSpheroidalEigenvalueDeriv(-2, 2, 2, 0.5 + 0.3im) isa ComplexF64
        @test SpinWeightedSpheroidalEigenvalue(-2, 3, 2, 0.4 + 0.2im; method=:Leaver) isa ComplexF64
        setprecision(BigFloat, 256) do
            s, l, m = -2, 4, 2
            γ = big"0.6" + big"0.4" * im
            λbig = SpinWeightedSpheroidalEigenvalue(s, l, m, γ)          # :Auto ⇒ Leaver for Complex{BigFloat}
            @test λbig isa Complex{BigFloat}
            @test precision(real(λbig)) == 256
            # BF256 Leaver strictly tighter than the ComplexF64 spectral seed, vs a 512-bit truth
            λref = setprecision(BigFloat, 512) do
                SpinWeightedSpheroidalEigenvalue(s, l, m, big"0.6" + big"0.4" * im)
            end
            err_big = abs(λbig - λref)
            err_f64 = abs(SpinWeightedSpheroidalEigenvalue(s, l, m, ComplexF64(γ); method=:SphericalExpansion) - λref)
            @test err_big < err_f64
            @test err_big < 1e-50
        end
    end

    # ---- spectral ↔ Leaver cross-agreement (complex γ, moderate |γ|) -----------------------------
    @testset "spectral ↔ Leaver cross-agreement (complex γ)" begin
        for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 4, 2), (2, 4, -3)]
            for γ in (0.4 + 0.3im, 1.0 - 0.5im, 1.2 + 0.2im)
                a = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:SphericalExpansion)
                b = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
                @test a ≈ b rtol = 1e-9              # observed ≤1.4e-14
            end
        end
    end

    # ---- dλ/dγ (complex) vs central finite difference -------------------------------------------
    @testset "dλ/dγ complex == finite difference" begin
        for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 4, 2)]
            for γ in (0.4 + 0.3im, 1.0 - 0.5im)
                dana = SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, γ)
                h = 1e-6
                dnum = (SpinWeightedSpheroidalEigenvalue(s, l, m, γ + h; method=:SphericalExpansion) -
                        SpinWeightedSpheroidalEigenvalue(s, l, m, γ - h; method=:SphericalExpansion)) / (2h)
                @test dana ≈ dnum rtol = 1e-6        # observed ≤4.5e-10
            end
        end
    end

    # ---- value tier: complex S(θ,φ), γ→0 limit, θ-derivative, real continuity --------------------
    @testset "value tier (complex γ): S(θ,φ) + θ-derivative" begin
        _angles = [(0.4, 0.0), (1.1, 0.9), (2.6, 1.7)]
        for (s, l, m) in [(-2, 2, 2), (-2, 3, -1), (0, 4, 2), (2, 4, -3)]
            S = SpinWeightedSpheroidalHarmonicS(s, l, m, 0.5 + 0.2im)
            @test eltype(S.coeffs) <: Complex
            @test S.eigenvalue isa Complex
            @test sum(abs2, S.coeffs) ≈ 1.0 rtol = 1e-12          # Euclidean-unit-norm coeffs
            for (θ, φ) in _angles
                @test S(θ, φ) isa Complex
                @test swsh_derivative(S, 0, 1, θ, φ) ≈ im * m * S(θ, φ) rtol = 1e-12   # ∂_φ = im·m
            end
        end
        # γ→0 value → {}_sY_lm (the seed)
        for (s, l, m) in [(-2, 2, 2), (0, 3, 0), (2, 4, -1)]
            S0 = SpinWeightedSpheroidalHarmonicS(s, l, m, 0.0 + 0.0im)
            for (θ, φ) in _angles
                @test S0(θ, φ) ≈ SpinWeightedSphericalHarmonicY(s, l, m, θ, φ) rtol = 1e-12
            end
        end
        # complex(x+0im) reproduces the real value/θ-derivative carrier
        for (s, l, m) in [(-2, 3, 2), (0, 4, -2)], x in (0.5, 1.1)
            Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, complex(x))
            Sr = SpinWeightedSpheroidalHarmonicS(s, l, m, x)
            for (θ, φ) in _angles
                @test Sc(θ, φ) ≈ Sr(θ, φ) rtol = 1e-10
                @test swsh_derivative(Sc, 1, 0, θ, φ) ≈ swsh_derivative(Sr, 1, 0, θ, φ) rtol = 1e-10
            end
        end
        # honest [FUTURE] exits: complex Leaver VALUE rep + complex γ-derivative both throw
        @test_throws ArgumentError SpinWeightedSpheroidalHarmonicS(-2, 3, 2, 0.5 + 0.2im; method=:Leaver)
        @test_throws ArgumentError swsh_gamma_derivative(-2, 3, 2, 0.5 + 0.2im, 1, 0, 0.7, 0.3)
    end
end
