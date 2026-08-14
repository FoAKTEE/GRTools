# test_qnm_fast.jl — Q1c: the `:fast` Kerr-QNM backend (shift-invert angular eigenpair + spin-
# sequence continuation), a native-Julia port of the `qnm` package speed engineering. Grounded vs
# (i) the `:leaver` reference over the campaign (both root-find the SAME `qnm_cf_residual`), (ii) the
# geev spectral λ (the shift-invert correctness anchor — SAME matrix elements), (iii) the baked
# qnm-package 55-row cross-fixture (ω + A_sep) + the a=0 overtone ladder, and (iv) analytic
# robustness (retrograde + m=0 continuations, backend-switch identity, BigFloat tier routing).
# Tolerances MEASURED via scratchpad/qqnm_phase3/q1c_env/{smoke,bench,measure}.jl, stated w/ headroom.

include(joinpath(@__DIR__, "fixtures", "qnm_fast_reference.jl"))   # QNM_FAST_REF (55), QNM_FAST_LADDER (11)

# deterministic LCG stream in [-1,1) — avoids an undeclared `using Random` in the test sandbox.
function _qnm_lcg(seed::UInt64)
    s = Ref(seed)
    return function ()
        s[] = 6364136223846793005 * s[] + 1442695040888963407    # UInt64 wraps
        return Float64(s[] >> 11) * (2.0^-53) * 2 - 1
    end
end

@testset "QNM :fast backend ω_{s,l,m,n}(a) (Q1c)" begin

    # ── (1) PRIORITY-1a: the hand-rolled pentadiagonal banded LU == dense solve ───────────────────
    @testset "pentadiagonal band LU == dense solve" begin
        rng = _qnm_lcg(0x1c1c1c1c2b2b2b2b)
        worst = 0.0
        for _ in 1:60
            n = 5 + (abs(rand(Int)) % 21)            # 5..25 (Int rand is fine; size only)
            M = zeros(ComplexF64, n, n)
            for i in 1:n, j in 1:n
                abs(i - j) <= 2 && (M[i, j] = complex(rng(), rng()))
            end
            b = [complex(rng(), rng()) for _ in 1:n]
            xref = M \ b
            A = copy(M); Perturbation._qnm_penta_lu!(A)
            x = copy(b); Perturbation._qnm_penta_solve!(A, x)
            worst = max(worst, sqrt(sum(abs2, x - xref)) / sqrt(sum(abs2, xref)))
        end
        @test worst <= 1e-11                          # measured 5.2e-15
    end

    # ── (2) G-Q1c-a (angular): shift-invert λ == geev spectral λ (SAME matrix) ────────────────────
    @testset "G-Q1c-a shift-invert λ == geev λ" begin
        worst = 0.0
        for (s, l, m, γ) in ((-2, 2, 2, 0.23 - 0.043im), (-2, 2, 2, 0.60 - 0.058im),
                             (-2, 3, 3, 0.94 - 0.059im), (-2, 4, 4, 1.47 - 0.051im),
                             (-2, 2, -2, 0.27 - 0.079im), (-2, 2, 0, 0.19 - 0.044im),
                             (+2, 2, 2, 0.60 - 0.058im))
            # geev reference (the :leaver angular path): λ = Harmonics spectral eigenvalue at γ=aω.
            λg = Perturbation._qnm_angular_lambda(s, l, m, 1.0, γ)      # a=1 so aω=γ
            λf, _, _, _, _ = Perturbation._qnm_angular_fast(s, l, m, γ, l * (l + 1) - s * (s + 1))
            worst = max(worst, abs(λf - λg))
        end
        @test worst <= 1e-11                          # measured 3.6e-14
    end

    # ── (3) G-Q1c-a: :fast == :leaver over the full campaign table (5 multiplets × 11 χ) ──────────
    @testset "G-Q1c-a :fast == :leaver (55-row campaign)" begin
        χs = vcat(0.0:0.1:0.9, 0.95)
        worst = 0.0
        for (l, m) in ((2, 2), (3, 3), (4, 4), (2, 0), (2, -2))
            ωL = qnm_frequency(-2, l, m, 0, 0.0)
            ωF = qnm_frequency(-2, l, m, 0, 0.0; backend=:fast)
            for χ in χs
                aL = qnm_frequency(-2, l, m, 0, χ; seed=ωL, backend=:leaver)
                aF = qnm_frequency(-2, l, m, 0, χ; seed=ωF, backend=:fast)
                worst = max(worst, abs(aL - aF) / max(abs(aL), 1e-30))
                ωL = aL; ωF = aF
            end
        end
        @test worst <= 1e-12                          # measured 2.0e-15 (both are CF-residual roots)
    end

    # ── (4) G-Q1c-b: :fast vs the baked qnm-package 55-row table (ω AND A_sep) ────────────────────
    @testset "G-Q1c-b :fast vs qnm 55-row table (ω + A_sep)" begin
        worstω = 0.0; worstA = 0.0
        for (l, m) in ((2, 2), (3, 3), (4, 4), (2, 0), (2, -2))
            ωseed = complex(0.0); Aseed = complex(0.0); have = false
            for (rl, rm, chi, ωref, Aref) in QNM_FAST_REF
                (rl == l && rm == m) || continue
                if !have
                    ωseed = ωref; Aseed = Aref; have = true    # a=0 anchor from the fixture (seed only)
                end
                r = Perturbation._qnm_polish_fast(-2, l, m, chi, ωseed, Aseed; ninv=0)
                @test r.residual <= 1e-10                       # a genuine CF root at this χ
                worstω = max(worstω, abs(r.ω - ωref))
                worstA = max(worstA, abs(r.A - Aref))
                ωseed = r.ω; Aseed = r.A                        # χ-continuation seed
            end
        end
        @test worstω <= 1e-10                          # measured 1.6e-12
        @test worstA <= 1e-10                          # measured 3.5e-12
    end

    # ── (5) G-Q1c-b: the a=0 (2,2,n) overtone ladder n=0..7 (:fast) ───────────────────────────────
    @testset "G-Q1c-b a=0 overtone ladder n=0..7" begin
        # per-n tol tracks the ladder's OWN accuracy (cf_err col: ~1e-16 at n≤2, ~1e-10 at n=3..7).
        tol = [1e-11, 1e-11, 1e-10, 1e-9, 1e-9, 1e-9, 1e-9, 1e-8]
        for n in 0:7
            ω, res = qnm_frequency(-2, 2, 2, n, 0.0; seed=QNM_FAST_LADDER[n+1] * (1 + 1e-3),
                                   backend=:fast, return_residual=true)
            @test abs(ω - QNM_FAST_LADDER[n+1]) <= tol[n+1]
            @test res <= 1e-10                          # genuine CF root
        end
    end

    # ── (6) G-Q1c-d: robustness — retrograde (2,-2,0) + m=0 (2,0,0) complete, no mode-jumps ───────
    @testset "G-Q1c-d robustness (retrograde + m=0 continuations)" begin
        for (l, m) in ((2, -2), (2, 0))
            seq = qnm_sequence_fast(-2, l, m, 0; a_max=0.99, build_interp=false)  # throws [BLOCKING] on fail
            @test seq.a[end] >= 0.99 - 1e-9                 # reached a_max
            @test seq.nnodes >= 200                          # a full adaptive sweep
            @test maximum(abs.(diff(seq.omega))) < 5e-3      # χ-continuity: no branch jumps (meas ≤7e-4)
            @test all(isfinite, seq.omega) && all(isfinite, seq.A)
            # re-polish a few interior spins: each is a real CF root
            for a in (0.3, 0.6, 0.9)
                idx = argmin(abs.(seq.a .- a))
                r = Perturbation._qnm_polish_fast(-2, l, m, seq.a[idx], seq.omega[idx], seq.A[idx]; ninv=0)
                @test r.residual <= 1e-10
            end
        end
    end

    # ── (7) G-Q1c-e: the backend switch — same answers, :leaver default, BigFloat tier routing ────
    @testset "G-Q1c-e backend switch + tier split" begin
        for (l, m, a) in ((2, 2, 0.5), (3, 3, 0.7), (2, -2, 0.9))
            sd = qnm_frequency(-2, l, m, 0, a; seed=0.5 - 0.08im)            # eikonal→continuation seed
            ωf = qnm_frequency(-2, l, m, 0, a; seed=sd, backend=:fast)
            ωl = qnm_frequency(-2, l, m, 0, a; seed=sd, backend=:leaver)
            ωd = qnm_frequency(-2, l, m, 0, a; seed=sd)                       # default
            @test abs(ωf - ωl) <= 1e-12
            @test ωd === ωl                                                    # default IS :leaver
        end
        @test_throws ArgumentError qnm_frequency(-2, 2, 2, 0, 0.5; backend=:nonsense)
        # BigFloat backend=:fast routes to the :leaver BF-256 home (precision inherited, not Float64-capped)
        setprecision(BigFloat, 256) do
            f64 = qnm_frequency(-2, 2, 2, 0, 0.5; seed=0.46 - 0.086im)
            ωbf = qnm_frequency(-2, 2, 2, 0, big"0.5"; seed=Complex{BigFloat}(f64), backend=:fast)
            @test ωbf isa Complex{BigFloat}
            @test precision(real(ωbf)) == 256
            @test abs(ComplexF64(ωbf) - f64) <= 1e-10
        end
    end

    # ── (8) adaptive l_max: the band grows with |γ|, tail component below tol ─────────────────────
    @testset "adaptive l_max (PRIORITY-4)" begin
        prev_lmax = 0
        for χ in (0.1, 0.5, 0.9, 0.95)
            ω = qnm_frequency(-2, 2, 2, 0, χ; seed=0.4 - 0.08im, backend=:fast)
            λ, A, tailfrac, lmax, iters = Perturbation._qnm_angular_fast(-2, 2, 2, χ * ω, complex(4.0))
            @test tailfrac <= 1e-9                     # the growth criterion is satisfied
            @test lmax >= prev_lmax                     # monotone in |γ|
            @test iters <= 10                           # RQI converges fast (measured ≤4)
            prev_lmax = lmax
        end
    end
end
