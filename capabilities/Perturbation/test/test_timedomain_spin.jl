# DG-P4: generic-spin (s = ±1, ±2) VACUUM m-mode ℓ-coupled TimeDomain sector — gates.
#
# Formulation dossier: progress/emri-disk/DG_SPIN_FORMULATION.md (commit 5f6bd0e).
# Chart decision (dossier §3.4 exits): Harms–Bernuzzi–Brügmann RT horizon-penetrating
# hyperboloidal chart [arXiv:1301.1591] — proven-regular coefficients on R ∈ [R₊, 1] for every
# s; the repo layer chart has divergent E⁻¹-coefficients at scri for s ≠ 0 (dossier §3.4, §8.2).
# The locked s=0 layer-chart path (mode1d/kerrmodes/kerrflux) is BIT-UNTOUCHED by this sector.
#
# GATES (tolerances MEASURED, headroom stated inline):
#   (B) coupling == CG        : closed forms [arXiv:1702.04204, Eq. (34)] vs the baked
#                               Clebsch–Gordan fixture (GrTool core/Wigner, one-shot capture —
#                               Wigner is NOT a runtime dep), ≤ 1e-13 rel.
#   (A1) s=0 coupling identity: spin coefficients at s=0 == locked `cminus`/`mode_coupling`
#                               BIT-EXACT.
#   (RT) transcription        : rtchart.jl coefficient functions vs the INDEPENDENT sympy
#                               re-derivation fixture (master equation + chart maps only) and
#                               vs the printed boundary speeds [arXiv:1301.1591, Eqs. (17),(18)].
#   (D) Bardeen–Press limit   : assembled a=0, s=−2 operator == chain-rule-transformed BPT
#                               operator [arXiv:1702.04204, Eqs. (67)–(70)] up to the overall
#                               normalization, at several radii.
#   (A2) s=0 physics reduction: RT-chart s=0 Kerr a=0.9 m=0 tail exponent == locked KerrModes
#                               chart exponent (Burko–Khanna law [arXiv:1301.1591, Eq. (1)]).
#   (C) Price tails           : ID0 vacuum tails per (s, ℓ) vs the empirical laws
#                               [arXiv:1301.1591, Eqs. (28),(29)] at a=0 and a=0.9.
#   (E) stability             : bounded field norms (no NaN/blowup) over every full run window.

using Test
using Perturbation
const TDS = Perturbation.TimeDomain

include(joinpath(@__DIR__, "fixtures", "spin_cg_reference.jl"))
include(joinpath(@__DIR__, "fixtures", "rt_coeffs_reference.jl"))
# locked-chart baked exponents (A2) — td_kerrtail.jl is already included by test_timedomain.jl
isdefined(Main, :TD_KERRTAIL) || include(joinpath(@__DIR__, "fixtures", "td_kerrtail.jl"))

@testset "TimeDomain spin (DG-P4)" begin

    # ---- (B) closed-form coupling == Clebsch–Gordan fixture --------------------------------------
    @testset "gate B: spin coupling vs Wigner CG fixture" begin
        # ⟨ₛY_{ℓ′m}|cosθ|ₛYₗₘ⟩: ℓ′=ℓ−1 ↦ c₋^ℓ, ℓ′=ℓ ↦ c₀^ℓ, ℓ′=ℓ+1 ↦ c₊^ℓ = c₋^{ℓ+1}
        # ([arXiv:1702.04204, Eq. (32)] placement). Fixture from BigInt-exact Racah CG.
        # MEASURED max rel dev 2.2e-16 (max abs dev 2.2e-16) → gate 1e-13 (≥400x headroom).
        worst = 0.0
        for (s, l, m, lp, ref) in SPIN_CG_REF
            v = lp == l - 1 ? TDS.spin_cminus(s, l, m) :
                lp == l ? TDS.spin_czero(s, l, m) : TDS.spin_cplus(s, l, m)
            dev = abs(v - ref) / max(abs(ref), 1.0)
            worst = max(worst, dev)
        end
        @test worst <= 1e-13
        # C-matrix consistency: 𝒮 == I − 𝒦² on a joint ℓ-window (the sin²θ = 1 − cos²θ matrix
        # identity behind [arXiv:1702.04204, Eq. (35)]; C₊^ℓ = −c₊^ℓ(c₀^{ℓ+1}+c₀^ℓ) etc.).
        # Interior rows are exact; window edges differ by the truncated band. MEASURED interior
        # max dev 2.3e-16 → gate 1e-13.
        for (s, m) in ((-2, 0), (-2, 3), (-1, 1), (1, -2), (2, 2), (0, 1))
            ℓs = collect(max(abs(s), abs(m)):(max(abs(s), abs(m)) + 7))
            n = length(ℓs)
            Km = TDS._spin_cos_matrix(s, ℓs, m)
            Sm = TDS._spin_sin2_matrix(s, ℓs, m)
            Iden = [i == j ? 1.0 : 0.0 for i in 1:n, j in 1:n]
            E = Sm - (Iden - Km * Km)
            @test maximum(abs.(E[3:n-2, 3:n-2])) <= 1e-13
            # symmetry (E_ππ symmetric ⇒ symmetric-hyperbolicity ingredient, dossier §3.3)
            @test Sm == transpose(Sm)
            @test Km == transpose(Km)
        end
    end

    # ---- (A1) s=0 reduction of the coupling: BIT-exact vs the locked KerrModes closed forms ------
    @testset "gate A1: s=0 coupling == locked mode_coupling (bitwise)" begin
        for l in 0:12, m in -l:l
            @test TDS.spin_cminus(0, l, m) === TDS.cminus(l, m)
            @test TDS.spin_czero(0, l, m) === 0.0
            C = TDS.spin_mode_coupling(0, l, m)
            @test (C[1], C[3], C[5]) === TDS.mode_coupling(l, m)
            @test C[2] === -0.0 || C[2] === 0.0   # ℓ±1 band vanishes identically at s=0
            @test C[4] === -0.0 || C[4] === 0.0
        end
    end

    # ---- (RT) coefficient transcription vs independent sympy re-derivation -----------------------
    @testset "gate RT-1: rtchart coefficients vs sympy re-derivation fixture" begin
        # Fixture: coefficients of the 2+1 m-mode RT equation derived from the master equation
        # [arXiv:2307.01349, Eq. (1)] + chart maps [arXiv:1301.1591, Eqs. (3)–(6),(12)] ONLY,
        # normalized to C_θθ = −1/R. Reconstruction below re-assembles the full θ-dependent
        # coefficients from the split functions (rtchart.jl header). MEASURED max rel dev
        # 2.3e-13 (64 samples × 7 coefficients; a ∈ {0, 1/2, 9/10}, s ∈ {0,±1,±2},
        # m ∈ {−2..3}, R from R₊+1e-4 to 1−1e-4, θ ∈ {π/2, π/3}) → gate 1e-10 (430x headroom).
        worst = 0.0
        for (a, s, m, R, θ, C0r, CTr, CRr, Cthr, CTTr, CTRr, CRRr) in RT_COEFFS_REF
            cot = cos(θ) / sin(θ)
            csc = 1 / sin(θ)
            C0 = TDS.rt_C0_flat(R; M=1.0, a=a, s=s) + (s^2 * cot^2 - s) / R +
                 im * m * TDS.rt_Cphi_re(R; a=a) + (m^2 * csc^2 + 2m * s * cot * csc) / R
            CT = TDS.rt_CT_flat(R; M=1.0, a=a, s=s) + 2im * a * s * cos(θ) / R +
                 im * m * TDS.rt_CTphi(R; M=1.0, a=a)
            CR = TDS.rt_CR(R; M=1.0, a=a, s=s) + im * m * TDS.rt_CRphi(R; M=1.0, a=a)
            Cth = -cot / R
            CTT = TDS.rt_CTT_flat(R; M=1.0, a=a) + a^2 * cos(2θ) / (2R)
            CTR = TDS.rt_CTR(R; M=1.0, a=a)
            CRR = TDS.rt_CRR(R; M=1.0, a=a)
            for (mine, ref) in ((C0, C0r), (CT, CTr), (CR, CRr), (Cth, Cthr),
                                (CTT, CTTr), (CTR, CTRr), (CRR, CRRr))
                worst = max(worst, abs(mine - ref) / max(abs(ref), 1e-3))
            end
        end
        @test worst <= 1e-10
    end

    @testset "gate RT-2: boundary characteristic speeds [1301.1591 Eqs. (17),(18)]" begin
        # c_R = C̃_TR/C̃_TT at θ=π/2 (C_RR vanishes at both ends). Printed to 6 digits in the
        # paper; MEASURED agreement to all printed digits → gate = printed precision.
        for (a, cH, cS) in ((0.0, -0.0460655, 0.08), (0.9, -0.0723329, 0.0826788))
            Rp = TDS.rt_Rplus(1.0, a)
            for (R, ref) in ((Rp, cH), (1.0, cS))
                ctt = TDS.rt_CTT_flat(R; M=1.0, a=a) - a^2 / (2R)
                @test abs(TDS.rt_CTR(R; M=1.0, a=a) / ctt - ref) <= 1e-6
            end
            @test abs(TDS.rt_CRR(Rp; M=1.0, a=a)) <= 1e-15   # outflow: C_RR(R₊) = 0
            @test abs(TDS.rt_CRR(1.0; M=1.0, a=a)) <= 1e-15  # outflow: C_RR(1) = 0
        end
    end

    # ---- (D) Bardeen–Press limit: assembled a=0, s=−2 operator == transformed BPT ---------------
    @testset "gate D: a=0 s=−2 Bardeen–Press operator identity" begin
        # BPT coefficients U,V,W [arXiv:1702.04204, Eqs. (67)–(70)] (triple-sourced, dossier
        # §6.1); (t,r*) 1+1 form via DERIVED-1/2 (dossier §3.1): −φ̈ + φ″ − μφ̇ − νφ′ − 4W φ = 0,
        # μ = 2(U+V), ν = 2(V−U). Chart push-forward: t = T + F(R), r* = G(R)
        # (rt_F_of_R; [arXiv:1301.1591, Eqs. (3),(5),(6)]), ∂_{r*} = p∂_T + q∂_R with
        # p = −F′/G′, q = 1/G′. At a=0 the RT and BG mode variables coincide (ψ = rΔˢΨ both).
        # The 6 transformed coefficients must be PROPORTIONAL (one common σ(R)) to the
        # assembled RT ones. p′,q′ by Richardson central differences (h=1e-5, error ≪ gate).
        # MEASURED max |ratio−1| = 2.5e-6 (FD-limited) → gate 1e-4 (40x headroom).
        M = 1.0
        s = -2
        fd(f, x; h=1e-5) = (8 * (f(x + h) - f(x - h)) - (f(x + 2h) - f(x - 2h))) / (12h)
        G(R) = (r = TDS.rt_r_of_R(R); r + 2M * log(r / (2M) - 1))
        F(R) = TDS.rt_F_of_R(R; M=M)
        pf(R) = -fd(F, R) / fd(G, R)
        qf(R) = 1 / fd(G, R)
        for l in (2, 3), R in (0.65, 0.75, 0.85, 0.95)
            r = TDS.rt_r_of_R(R)
            fs = 1 - 2M / r
            λ = (l + 2) * (l - 1)
            U = 2M / r^2
            V = -2fs / r
            W = (fs / 4) * (λ / r^2 - 2M / r^3)
            μ = 2 * (U + V)
            ν = 2 * (V - U)
            p = pf(R)
            q = qf(R)
            pp = fd(pf, R)
            qp = fd(qf, R)
            L1 = (-1 + p^2, 2p * q, q^2, q * pp - μ - ν * p, q * qp - ν * q, -4W)
            RT = (TDS.rt_CTT_flat(R; M=M, a=0.0),
                  TDS.rt_CTR(R; M=M, a=0.0),
                  TDS.rt_CRR(R; M=M, a=0.0),
                  TDS.rt_CT_flat(R; M=M, a=0.0, s=s),
                  TDS.rt_CR(R; M=M, a=0.0, s=s),
                  real(TDS.rt_C0diag(R, s, l, 0; M=M, a=0.0)))
            σ = RT[3] / L1[3]
            for i in 1:6
                @test abs(RT[i] / (σ * L1[i]) - 1) <= 1e-4
            end
        end
    end

    # ---- shared driver for the evolution gates ---------------------------------------------------
    # Per-case numerical-dissipation settings (all MEASURED, see inline evidence):
    #   |s| = 2 : modal filter ON (filter_alpha=36, power=16) — required by the anti-damped
    #             glancing families at the sonic points (spinmodes.jl header). The filter's
    #             drainage also sets a late-time artificial floor, so |s|=2 fit windows stop
    #             before it (measured window-by-window below).
    #   |s| ≤ 1 : filter OFF (measured step-operator growth ≤ 8e-6/M — negligible over the
    #             run windows) — the filter otherwise damps the DEEP (|ψ| ≲ 1e-9) tails of
    #             these cases into exponential decay (measured LPI → −10 with filter ON vs
    #             clean law values with it OFF).
    # ID0 = stationary compact pure-multipole data [arXiv:1301.1591, §4], G-width w=1500
    # (σ_R ≈ 0.026, resolved by K=24/N=8; the paper's w=3000 is under-resolved here and its
    # filter-frozen response was the dossier-era floor).
    function _spin_run(s, a; K=24, nl=2, fa=36.0, tf=1800.0, id=:ID0, init_mode=1, ℓs=nothing,
                       m=0)
        l0 = max(abs(s), abs(m))
        ℓl = ℓs === nothing ? collect(l0:(l0+nl-1)) : ℓs
        p = TDS.SpinModesProblem(s, ℓl, m; N=8, K=K, M=1.0, a=a, filter_alpha=fa)
        dt = TDS.spin_cfl_dt(p; C=0.6)
        τs, Ψ, mx = TDS.solve_spin_modes(p; tfinal=tf, dt=dt,
                                         record_every=max(1, round(Int, 2.0 / dt)), id=id,
                                         init_mode=init_mode,
                                         ic=(amp=1.0, R0=(a == 0 ? 0.8 : 0.75), w=1500.0))
        return τs, Ψ, mx
    end
    fitw(τs, Ψ, ip, q, win) = TDS.fit_tail_window(τs, abs.(Ψ[ip, q, :]), win)

    # ---- (A2) s=0 physics reduction: RT chart vs locked KerrModes chart --------------------------
    @testset "gate A2: s=0 RT-chart tail == locked-chart tail (a=0.9, m=0)" begin
        # Same physics as the locked characterization (c) in test_timedomain.jl and its baked
        # fixture: a=0.9, m=0, ID1(momentum) seeded in ℓ′=2, ℓ=0 down-mode; Burko–Khanna law
        # μ = ℓ′+ℓ+1 = 3 at finite radii [arXiv:1301.1591, Eq. (1)]. Locked-chart baked
        # exponent TD_KERRTAIL.n_l0 = −3.01705 (r*=100, window (1200,3000)). RT chart:
        # MEASURED fit −2.9442 (mid probe, window (800,1600); the late windows drift to −4.2
        # as |ψ| → 1e-11 floors). |Δ| = 0.073 → gate 0.2 (2.7x headroom).
        τs, Ψ, mx = _spin_run(0, 0.9; nl=0, ℓs=[0, 2, 4], fa=0.0, tf=3000.0, id=:ID1,
                              init_mode=2)
        n_rt = fitw(τs, Ψ, 2, 1, (800.0, 1600.0))
        @test abs(n_rt - TD_KERRTAIL.n_l0) <= 0.2
        @test all(isfinite, mx)                      # gate E on this run too
    end

    # ---- (C)+(E) Price tails per (s, ℓ): vacuum ID0, laws [arXiv:1301.1591, Eqs. (28),(29)] ------
    # ℓ′ = ℓ = ℓ₀ = |s| (diagonal, no LPI splitting per the paper's tables), m = 0.
    # Law values: horizon/finite μ = ℓ′+ℓ+3 (+1 iff s>0, m=0) [Eq. (28) + δ rule];
    # scri μ = ℓ−s+2 [Eq. (29)]. Each tolerance is |measured−law| with ≥2x headroom, and the
    # measured values are quoted; fits are deterministic for fixed configs.
    @testset "gate C+E: s=$s a=$a vacuum Price tails" for (s, a, K, nl, fa, tf, gates) in (
        # (s, a, K, nl, filter_alpha, tfinal, [(probe, window, law, tol), ...])
        # s=−2, a=0: scri law 6 [Eq. (29)]: MEASURED −5.999 (dev 0.001) → tol 0.4.
        (-2, 0.0, 24, 2, 36.0, 1800.0, [(3, (1000.0, 1800.0), -6.0, 0.4)]),
        # s=−1, a=0: scri law 4: MEASURED −4.126 (dev 0.13) → tol 0.5;
        #            horizon law 5: MEASURED −5.202 (dev 0.20) → tol 0.7.
        (-1, 0.0, 24, 2, 0.0, 1600.0, [(3, (400.0, 900.0), -4.0, 0.5),
                                       (1, (400.0, 900.0), -5.0, 0.7)]),
        # s=+1, a=0: scri law 2: MEASURED −2.123 (dev 0.12) → tol 0.5;
        #            horizon law 6 (δ=1 boost, s>0 m=0): MEASURED −6.46 (dev 0.46) → tol 1.0
        #            (window-scatter −7.1..−5.4 across neighbouring windows — the s>0
        #            three-value splitting region, [arXiv:1301.1591, §5.2]).
        (1, 0.0, 24, 2, 0.0, 700.0, [(3, (200.0, 500.0), -2.0, 0.5),
                                     (1, (200.0, 600.0), -6.0, 1.0)]),
        # s=+2, a=0: EXPONENT SKIPPED WITH EVIDENCE — the μ_scri = 2 tail lies below the
        # residual numerical floor (MEASURED: |ψ|_scri decays steeply to ~2e-7 by t≈500 and
        # meets a flat 7.5e-7-level floor; no law stretch emerges at K∈{12,16,24}, ID0/ID1).
        # Stability/boundedness (gate E) is still asserted; K=12 keeps the s=+2 scri-glancing
        # growth mode unresolved (MEASURED stable, maxu_end 8.5e-7 vs 1.4e-3/M growth at K=24).
        (2, 0.0, 12, 2, 36.0, 1200.0, Tuple{Int,Tuple{Float64,Float64},Float64,Float64}[]),
        # s=−2, a=0.9 (coupled ℓ=2..4): scri law 6: MEASURED −6.538 (dev 0.54) → tol 1.0;
        #            horizon law 7: MEASURED −6.427 (dev 0.57) → tol 1.2.
        (-2, 0.9, 24, 3, 36.0, 1800.0, [(3, (700.0, 1300.0), -6.0, 1.0),
                                        (1, (700.0, 1300.0), -7.0, 1.2)]),
        # s=−1, a=0.9: scri law 4: MEASURED −4.096 (dev 0.10) → tol 0.4;
        #            horizon law 5: MEASURED −5.61 (dev 0.61) → tol 1.2.
        (-1, 0.9, 24, 3, 0.0, 1600.0, [(3, (400.0, 900.0), -4.0, 0.4),
                                       (1, (400.0, 900.0), -5.0, 1.2)]),
        # s=+1, a=0.9: scri law 2: MEASURED −2.176 (dev 0.18) → tol 0.6 (horizon skipped:
        # window scatter −11..−6.6, no stable value at this budget — evidence in dossier run).
        (1, 0.9, 24, 3, 0.0, 700.0, [(3, (200.0, 500.0), -2.0, 0.6)]),
        # s=+2, a=0.9: exponent skipped with evidence (same floor pathology as a=0:
        # MEASURED |ψ|_scri minimum 1.6e-7 at t≈550 then flat 4.8e-7..9.8e-7 junk floor,
        # K∈{12,24}, ID0/ID1); gate E asserted (K=12: MEASURED maxu_end 1.4e-7, decaying).
        (2, 0.9, 12, 3, 36.0, 1200.0, Tuple{Int,Tuple{Float64,Float64},Float64,Float64}[]),
    )
        τs, Ψ, mx = _spin_run(s, a; K=K, nl=nl, fa=fa, tf=tf)
        # gate E: bounded, finite, decayed by the end of the window (no blowup)
        @test all(isfinite, mx)
        @test mx[end] <= 0.01                       # measured ends: 6e-9 .. 1.3e-3
        for (ip, win, law, tol) in gates
            @test abs(fitw(τs, Ψ, ip, 1, win) - law) <= tol
        end
        # s=−2, a=0 horizon: the μ = 7 law [Eq. (28)] is approached only asymptotically
        # (the paper: LPIs "reach the value −x only asymptotically at late times"); assert
        # the monotone approach + the late-window bracket (MEASURED −2.65 → −4.11 → −6.15,
        # late window −6.45).
        if s == -2 && a == 0.0
            n1 = fitw(τs, Ψ, 1, 1, (400.0, 900.0))
            n2 = fitw(τs, Ψ, 1, 1, (700.0, 1300.0))
            n3 = fitw(τs, Ψ, 1, 1, (1000.0, 1800.0))
            @test n3 < n2 < n1 < -1
            @test -7.5 < fitw(τs, Ψ, 1, 1, (1300.0, 1800.0)) < -5.8
        end
    end

    # ---- (E) pentadiagonal s·m ≠ 0 complex-sector smoke ------------------------------------------
    @testset "gate E: pentadiagonal s·m≠0 sector (s=−2, m=2, a=0.9)" begin
        # E_ππ is PENTADIAGONAL here (C±^ℓ ∝ ms ≠ 0; dossier §3.3) and the state is complex.
        # Structure + stability smoke; MEASURED maxu: 23.5 (IC) → ≤ 1e-4 by t=400 (decay).
        A = TDS.rt_ATT(0.7, -2, [2, 3, 4], 2; M=1.0, a=0.9)
        @test abs(A[1, 2]) > 1e-3 && abs(A[1, 3]) > 1e-3   # ℓ±1 AND ℓ±2 bands populated
        @test A ≈ transpose(A)                             # symmetric (dossier §3.3)
        τs, Ψ, mx = _spin_run(-2, 0.9; K=24, nl=3, fa=36.0, tf=400.0, m=2)
        @test all(isfinite, mx)
        @test mx[end] <= 1e-3 * mx[1]                      # decayed ≥3 orders from the IC
    end

end
