# Module — SpheroidalEigenvalues residual fold: resolve the swsh_eigenvalue_big [HOLE].
# Folds KerrLorenzCirc `precision.jl/swsh_eigenvalue_big` (the (λ, residual) BigFloat refinement)
# into the canonical Harmonics Leaver path so there is ONE BigFloat eigenvalue home:
#   _swsh_eigenvalue_leaver_residual(s,l,m,γ)  -> (λ, |Eq(A)|)   [single CF body]
#   _swsh_eigenvalue_leaver(s,l,m,γ)           -> λ              [bare-λ wrapper, AHA]
#   SpinWeightedSpheroidalEigenvalue(...; return_residual=true)  -> (λ, |Eq(A)|)  [front door]
#
# Grounding:
#  (A) λ vs the LIGHT live oracle SpinWeightedSpheroidalHarmonics (LinearAlgebra-only) — the
#      returned λ EQUALS SWSH._swsh_eigenvalue_leaver / SpinWeightedSpheroidalEigenvalue(:Leaver)
#      bit-for-bit (same CF body), incl edges s∈{-2,0,2}, m<0, l=|m|.
#  (B) residual contract grounded CLOSED-FORM — residual = |Eq(A)| at the secant root is the
#      Leaver-equation defect, hence ≈ machine floor: residual < eps(T)*C. C=256 bounds the
#      observed 0–5.25×eps(T) across the cases (secant rtol = eps(T)*8); no oracle needed.
#  (C) 256-vs-512-bit self-convergence — λ@512 vs λ@256 agree to ~eps(BigFloat@256), confirming
#      the Leaver path INHERITS the active setprecision (BigFloat inheritance).
#  (D) ONE-SHOT KLC literal fixture — for the representative (s,l,m,γ)=(-2,4,2,big"0.5") at
#      256 bits, KLC swsh_eigenvalue_big returns (λ, residual) with λ pinned below (captured from
#      the SWSH :Leaver oracle, which the migrated KLC AUTO-ROUTES to — see KLC precision.jl body).
#      Pins the (λ, residual) TUPLE SHAPE + value; no live KLC test-dep (KLC pulls a heavy closure).

using Test
using Harmonics
# Legacy oracle SpinWeightedSpheroidalHarmonics._swsh_eigenvalue_leaver BAKED into a fixture
# (captured once; same args as the loops below). Same characterization at the SAME `==`.
include(joinpath(@__DIR__, "fixtures", "spheroidal_eigenvalue_residual_ref.jl"))  # FIXRES_leaver

@testset "SpheroidalEigenvalues residual fold" begin

    _cases = [(s, l, m) for s in (-2, 0, 2) for l in (max(abs(s), 2):max(abs(s), 2)+2)
              for m in (-l, -1, 0, 1, l) if abs(m) <= l]

    # ---- (A) λ from the residual variant EQUALS the bare-λ leaf & legacy oracle ---------------
    @testset "λ == bare leaf == legacy oracle (full precision)" begin
        for (s, l, m) in _cases, γ in (0.2, 0.9, 3.0)
            λres, _ = Harmonics._swsh_eigenvalue_leaver_residual(s, l, m, γ)
            @test λres == Harmonics._swsh_eigenvalue_leaver(s, l, m, γ)   # single CF body
            @test λres == FIXRES_leaver[(s, l, m, γ)]                     # legacy CF, bit-identical (baked)
        end
    end

    # ---- (B) residual is the Leaver-equation defect at the root ⇒ ≈ machine floor -------------
    @testset "residual ≈ machine floor (closed form: residual < eps(T)*C)" begin
        for (s, l, m) in _cases, γ in (0.2, 0.9, 3.0)
            _, res = Harmonics._swsh_eigenvalue_leaver_residual(s, l, m, γ)
            @test res isa Float64
            @test res >= 0
            @test res < eps(Float64) * 256          # observed 0–5.25×eps; secant rtol = eps*8
        end
    end

    # ---- (A/front door) return_residual kwarg: tuple, λ matches scalar return -----------------
    @testset "return_residual front door tuple == (scalar-λ, floor residual)" begin
        for (s, l, m) in [(-2, 3, 2), (0, 4, -1), (2, 5, 0)], γ in (0.3, 1.2)
            t = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; return_residual=true)
            @test t isa Tuple && length(t) == 2
            λscalar = SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Leaver)
            @test t[1] == λscalar
            @test t[2] < eps(Float64) * 256
        end
        # default scalar return unchanged (no tuple)
        @test SpinWeightedSpheroidalEigenvalue(0, 2, 0, 0.7) isa Real
        # γ=0 front-door: analytic λ, zero residual
        t0 = SpinWeightedSpheroidalEigenvalue(-2, 4, 2, 0.0; return_residual=true)
        @test t0 == (4 * 5 - (-2) * (-1), 0.0)
    end

    # ---- (C) BigFloat inheritance: 256-vs-512-bit self-convergence ----------------------------
    @testset "BigFloat 256-vs-512 self-convergence (precision inheritance)" begin
        setprecision(BigFloat, 256) do
            s, l, m = -2, 4, 2
            λ256, res256 = Harmonics._swsh_eigenvalue_leaver_residual(s, l, m, big"0.5")
            @test λ256 isa BigFloat
            @test res256 < eps(BigFloat) * 256                    # floor tracks the active 256-bit precision
            λ512 = setprecision(BigFloat, 512) do
                Harmonics._swsh_eigenvalue_leaver(s, l, m, big"0.5")
            end
            relerr = abs(λ256 - λ512) / abs(λ512)
            @test relerr < eps(BigFloat) * 16                     # agrees to ~256-bit eigenvalue floor
            @test eps(BigFloat) < 1e-70                           # 256-bit floor ≪ eps(Float64)
        end
    end

    # ---- (D) ONE-SHOT KLC literal fixture: (λ, residual) tuple shape + pinned value -----------
    @testset "KLC swsh_eigenvalue_big fixture (s,l,m,γ)=(-2,4,2,0.5) @256" begin
        # λ literal captured one-shot from the SWSH :Leaver oracle (which the migrated KLC
        # AUTO-ROUTES to per packages/KerrLorenzCirc/src/precision.jl), 256-bit:
        λ_klc = big"15.76206790420670136761911766362038260105787843587709784382231630857785749299388"
        setprecision(BigFloat, 256) do
            t = SpinWeightedSpheroidalEigenvalue(-2, 4, 2, big"0.5"; return_residual=true)
            @test t isa Tuple && length(t) == 2                   # KLC (λ, residual) tuple shape
            λ, res = t
            @test λ isa BigFloat
            @test abs(λ - λ_klc) / abs(λ_klc) < eps(BigFloat) * 16
            @test res < eps(BigFloat) * 256                       # residual = |Eq(A)| ≈ floor
        end
    end
end
