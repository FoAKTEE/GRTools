# DG-S2: SOURCED s = ±2 / s = 0 circular-equatorial point-particle sector at KERR a ≠ 0 — gates.
#
# Physics: closes the ONE [BLOCKED] item of the committed assessment (dossier §8.3): the
# spheroidal→spherical source re-projection. Route A (spinsource_kerr.jl provenance header):
# per-spheroidal jumps from the general-a GHP triplet + the DG-P5 jump algebra, b-projected
# onto the spherical TD band with Harmonics' SphericalExpansion coefficients, pushed through
# the Kerr RT chart (DERIVED-CK), evolved by the UNCHANGED DG-P5 machinery, and extracted per
# spheroidal mode by b back-projection.
#
# GATES (tolerances MEASURED on the exact configs below — scratchpad/dgs2_measure1.out +
# scratchpad/dgs2_flux_campaign.dat; headroom stated inline):
#   (K1) Kerr chart push-forward   : analytic H_a′(R) vs 4th-order central differences at
#        a ∈ {0.5, ±0.9}; H_a ≡ rt_F_of_R at a = 0 (closed-form identity).
#   (K2) b-representation sanity   : SphericalExpansion equator values/derivatives vs the
#        INDEPENDENT Leaver representation at the gate γ = aω; band orthonormality of b;
#        the ℓ ↔ ℓ±1 mixing is present and perturbative.
#   (K3) a → 0 reduction IDENTITY  : the Kerr jump chain at a = 0 equals the committed
#        `spin_source_jumps` and the committed `SpinSourcedProblem` jump state (measured
#        ≤ 2e-16 — same arithmetic up to λ-from-eigen-solve; the b-matrix is EXACTLY I).
#   (K4) jumps == FD Green's fn    : [X], [X′] == Z^∞R̃_up − Z^H R̃_in at r₀ for a ∈
#        {0.5, ±0.9} incl. s = +2 (unit-transmission Radial + oracle amplitudes — the
#        DECISIVE per-spheroidal-mode check of the general-a source + jump conventions).
#   (K5) DECISIVE TD flux vs FD oracle at a ≠ 0: s = −2 circular r₀ = 10M; (2,2) at
#        a = 0.9 through the (6,32) → (8,48) resolution ladder with ladder-convergence
#        asserted; (3,3) at a = 0.9; (2,2) at a = 0.5 with two-level convergence.
#   (K6) SKH literature anchor     : the FD-oracle |m|-mode sums and the TD band totals vs
#        the PUBLISHED fluxes of Sundararajan-Khanna-Hughes [arXiv:gr-qc/0703028, Table I
#        (label FDm2, |m| = 2) row (a=0.9, r0=10) Ė_FD = 4.456e-5 and row (a=0.5, r0=10)
#        Ė_FD = 4.812e-5; Table II (FDm3, |m| = 3) row (a=0.9, r0=10) Ė_FD = 5.044e-6] —
#        an ~0.1%-digit table (anchor, not a tight gate).
#
# The committed a = 0 sourced path and all vacuum paths are BIT-UNTOUCHED (additive
# dispatch only; asserted implicitly by the untouched DG-P5 suite running in this session).
# Honest boundaries: |s| = 1 sourced [FUTURE]; m = 0 fluxless; eccentric/inclined Kerr TD
# out of scope (monochromatic jump ansatz); (2,1)-class small-|Z| modes need (10,64)+
# resolution as at a = 0 (campaign evidence, not suite-gated).

using Test
using LinearAlgebra: I as IdK5
using Perturbation
using Perturbation: TeukolskyPointParticleModeCircular, TeukolskyRadial, teukolsky_radial_deriv
using Harmonics: SpinWeightedSpheroidalHarmonicS, swsh_derivative, expansion_coefficients
const TDS6 = Perturbation.TimeDomain

@testset "TimeDomain Kerr spin source (DG-S2)" begin

    # ---- (K1) Kerr chart push-forward -------------------------------------------------------------
    @testset "gate K1: Kerr RT chart maps (DERIVED-CK)" begin
        # MEASURED (dgs2_measure1.out): dH/dR vs 4th-order FD rel ≤ 6.9e-11 (FD-limited,
        # same class as the a=0 J1 gate) → gate 1e-7; a=0 identity ≤ 1.3e-16 → gate 1e-13.
        for a in (0.5, 0.9, -0.9), r in (8.0, 10.0)
            R = TDS6.rt_R_of_r(r)
            h = 1e-6
            fd = (8 * (TDS6.rt_H_of_R_kerr(R + h; a=a) - TDS6.rt_H_of_R_kerr(R - h; a=a)) -
                  (TDS6.rt_H_of_R_kerr(R + 2h; a=a) - TDS6.rt_H_of_R_kerr(R - 2h; a=a))) / (12h)
            @test abs(TDS6._rt_dH_dR_kerr(R; a=a) - fd) <= 1e-7 * abs(fd)
        end
        for r in (6.0, 10.0, 20.0)
            R = TDS6.rt_R_of_r(r)
            @test abs(TDS6.rt_H_of_R_kerr(R; a=0.0) - TDS6.rt_F_of_R(R)) <=
                  1e-13 * abs(TDS6.rt_F_of_R(R))
        end
        @test TDS6._kerr_phi_shift(10.0, 0.0) == 0.0        # φ̃ shift vanishes at a = 0
    end

    # ---- (K2) b-representation sanity -------------------------------------------------------------
    @testset "gate K2: spheroidal→spherical b-coefficients at gate γ = aω" begin
        # MEASURED: SphExp vs Leaver equator values ≤ 1.3e-15 / ∂θ ≤ 5.9e-14 → gate 1e-11
        # (≥170x); band-orthonormality residual ≤ 8.9e-16 → gate 1e-12; the ℓ̂ = 2 → ℓ = 3
        # mixing coefficient is 1.04e-2 (a=0.9) / 5.9e-3 (a=0.5): present AND perturbative.
        for (s, m, a, r0, jlo) in ((-2, 2, 0.9, 10.0, 2), (-2, 3, 0.9, 10.0, 3),
                                   (-2, 2, 0.5, 10.0, 2))
            J = TDS6.spin_source_jumps_kerr(s, jlo, m, a, r0)
            γ = a * J.omega
            B = Dict{Int,Dict{Int,Float64}}()
            for j in jlo:(jlo+5)
                Sc = SpinWeightedSpheroidalHarmonicS(s, j, m, γ)
                B[j] = Dict(k => Float64(v) for (k, v) in expansion_coefficients(Sc))
                if j <= jlo + 1     # cross-representation check (independent Leaver series)
                    ScL = SpinWeightedSpheroidalHarmonicS(s, j, m, γ; method=:Leaver)
                    @test abs(Sc(pi / 2, 0.0) - ScL(pi / 2, 0.0)) <=
                          1e-11 * abs(ScL(pi / 2, 0.0))
                    @test abs(swsh_derivative(Sc, 1, 0, pi / 2, 0.0) -
                              swsh_derivative(ScL, 1, 0, pi / 2, 0.0)) <=
                          1e-11 * abs(swsh_derivative(ScL, 1, 0, pi / 2, 0.0))
                end
            end
            worst = 0.0
            for j1 in jlo:(jlo+5), j2 in jlo:(jlo+5)
                acc = sum(get(B[j1], l, 0.0) * get(B[j2], l, 0.0)
                          for l in max(abs(s), abs(m)):(jlo + 11))
                worst = max(worst, abs(acc - (j1 == j2 ? 1.0 : 0.0)))
            end
            @test worst <= 1e-12
            bmix = abs(get(B[jlo], jlo + 1, 0.0))
            @test 1e-3 < bmix < 0.1                        # mixing present, perturbative
        end
    end

    # ---- (K3) a → 0 reduction identity ------------------------------------------------------------
    @testset "gate K3: Kerr jump chain at a=0 == committed spin_source_jumps" begin
        # MEASURED: ≤ 6.9e-17 rel across all fields/configs (identical arithmetic; the
        # only path differences are λ-from-eigen-solve and Ω-convention) → gate 1e-13
        # (>1000x, guards cross-machine libm variation).
        rel(x, y) = abs(x - y) / max(abs(y), floatmin(Float64))
        for (s, l, m, r0) in ((-2, 2, 2, 10.0), (-2, 3, 3, 8.0), (2, 2, 2, 10.0),
                              (0, 2, 2, 10.0), (-2, 2, 1, 10.0))
            J0 = TDS6.spin_source_jumps(s, l, m, r0)
            JK = TDS6.spin_source_jumps_kerr(s, l, m, 0.0, r0)
            @test rel(JK.omega, J0.omega) <= 1e-13
            @test rel(JK.Jpsi, J0.Jpsi) <= 1e-13 || abs(J0.Jpsi) == 0 && abs(JK.Jpsi) == 0
            @test rel(JK.Jpsip, J0.Jpsip) <= 1e-13
            @test rel(JK.JXp, J0.JXp) <= 1e-13
            if abs(J0.psid) == 0
                @test abs(JK.psid) == 0                    # s = 0: no δ-layer, exactly
            else
                @test rel(JK.psid, J0.psid) <= 1e-13
            end
        end
        # constructor-level: the assembled RT jump state and mesh addressing coincide
        # (MEASURED ≤ 1.7e-16 → gate 1e-13; b-matrix is EXACTLY the identity at a = 0).
        p0 = TDS6.SpinSourcedProblem(-2, [2], 2, 10.0; N=4, K=16)
        pk = TDS6.SpinSourcedKerrProblem(-2, [2], 2, 0.0, 10.0; N=4, K=16)
        @test rel(pk.prob.Jψ[1], p0.Jψ[1]) <= 1e-13
        @test rel(pk.prob.Jπ[1], p0.Jπ[1]) <= 1e-13
        @test rel(pk.prob.Jχ[1], p0.Jχ[1]) <= 1e-13
        @test pk.prob.omega == p0.omega
        @test pk.prob.kp == p0.kp
        pkb = TDS6.SpinSourcedKerrProblem(-2, [2, 3, 4, 5], 2, 0.0, 10.0; N=4, K=16)
        @test maximum(abs.(pkb.bmat - Matrix(1.0 * IdK5, 4, 4))) == 0.0
    end

    # ---- (K4) jumps == FD Green's-function representation at a ≠ 0 --------------------------------
    @testset "gate K4: Kerr jumps vs Z^∞R_up − Z^H R_in (FD construction, a ≠ 0)" begin
        # THE decisive per-spheroidal-mode check of the general-a source triplet + jump
        # algebra + λ/Ω/Υt conventions. MEASURED rel: 3.4e-13 (−2,2,2,a=0.9), 2.7e-13
        # (−2,3,3,0.9), 1.3e-13 (−2,2,2,0.5), 4.5e-14 (+2,2,2,0.9), 6.7e-14 (−2,2,2,−0.9)
        # — MST-series-limited → gate 1e-10 (≥290x headroom, J3-precedent class).
        for (s, l, m, a, r0) in ((-2, 2, 2, 0.9, 10.0), (-2, 3, 3, 0.9, 10.0),
                                 (-2, 2, 2, 0.5, 10.0), (2, 2, 2, 0.9, 10.0),
                                 (-2, 2, 2, -0.9, 10.0))
            J = TDS6.spin_source_jumps_kerr(s, l, m, a, r0)
            o = TeukolskyPointParticleModeCircular(s, l, m, a, r0)
            @test abs(J.lambda - o.lambda) <= 1e-10 * (1 + abs(o.lambda))
            R = TeukolskyRadial(s, l, m, a, J.omega; lambda=o.lambda)
            JX_FD = o.Zinf * R.Up(r0) - o.ZH * R.In(r0)
            JXp_FD = o.Zinf * teukolsky_radial_deriv(R.Up, r0) -
                     o.ZH * teukolsky_radial_deriv(R.In, r0)
            @test abs(J.JX - JX_FD) <= 1e-10 * abs(JX_FD)
            @test abs(J.JXp - JXp_FD) <= 1e-10 * abs(JXp_FD)
        end
    end

    # ---- constructor guards (honest boundaries) ---------------------------------------------------
    @testset "Kerr sourced constructor guards" begin
        @test_throws ArgumentError TDS6.SpinSourcedKerrProblem(-1, [2], 2, 0.9, 10.0)  # |s|=1 [FUTURE]
        @test_throws ArgumentError TDS6.SpinSourcedKerrProblem(-2, [2], 0, 0.9, 10.0)  # m=0 fluxless
        @test_throws ArgumentError TDS6.SpinSourcedKerrProblem(-2, [2], 2, 1.0, 10.0)  # subextremal
        @test_throws ArgumentError TDS6.spin_source_jumps_kerr(-2, 2, 2, 0.9, 1.3)     # no timelike circular geodesic
        @test_throws ArgumentError TDS6.SpinSourcedKerrProblem(-2, [2, 4], 2, 0.9, 10.0) # band must be consecutive (s·a≠0)
        p = TDS6.SpinSourcedKerrProblem(-2, [2, 3], 2, 0.9, 10.0; N=4, K=16)
        R0 = TDS6.rt_R_of_r(10.0)
        @test any(abs.(p.prob.base.mesh.edges .- R0) .<= 1e-14)                        # forced edge
        @test p.prob.base.mesh.edges[p.prob.kp+1] ≈ R0 atol = 1e-14
        @test p.rplus ≈ 1 + sqrt(1 - 0.81) atol = 1e-15
    end

    # ---- shared Kerr flux driver -------------------------------------------------------------------
    function _kerr_flux_case(s, m, a, r0, lmax; N, K, tf=1000.0, twin=(0.7, 1.0))
        band = collect(max(abs(s), abs(m)):lmax)
        kp = TDS6.SpinSourcedKerrProblem(s, band, m, a, r0; N=N, K=K)
        dt = TDS6.spin_cfl_dt(kp.prob.base; C=0.6)
        τs, Ψ, Π, maxu = TDS6.solve_spin_sourced(kp; tfinal=tf, dt=dt,
                                                 record_every=max(1, round(Int, 5.0 / dt)))
        @test all(isfinite, maxu)                          # boundedness on every run
        return TDS6.spin_sourced_flux_kerr(kp, τs, Π; twin=twin)
    end

    # ---- (K5) DECISIVE: TD flux vs FD oracle at a ≠ 0 ---------------------------------------------
    @testset "gate K5a: s=-2 (2,2) a=0.9 r0=10 ladder (6,32)→(8,48) vs FD oracle" begin
        # MEASURED (scratchpad/dgs2_runs/c01,c02.dat; band ℓ ∈ [2,4], tf=1000, twin=(0.7,1)):
        #   (N,K)=(6,32): EI rel +1.101e+0, EH rel −2.43e-2, devI 3.8e-2 — the a=0.9 scri
        #     signal is UNRESOLVED at this level (resolution cliff vs a=0's −1.1e-1 and
        #     a=0.5's +5.1e-2: honest ladder floor, gated loose as a convergence point);
        #   (N,K)=(8,48): EI rel −3.31e-2, EH rel +2.49e-4 — resolution-limited, NOT
        #     band-limited (lmax=4 vs 5 shifts EI rel by only 2.6e-5, c08.dat), and the
        #     (10,64) campaign level continues the convergence (figure data). Gates at
        #     measured + ≥2x headroom; ladder convergence asserted on Ė^∞.
        o = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.9, 10.0)
        errs = Float64[]
        for (N, K, tolI, tolH, dtol) in ((6, 32, 2.5, 6e-2, 0.08), (8, 48, 7e-2, 1e-3, 0.05))
            fx = _kerr_flux_case(-2, 2, 0.9, 10.0, 4; N=N, K=K)
            f2 = fx[1]                                     # ℓ̂ = 2 back-projection
            push!(errs, abs(f2.EI / o.EnergyFlux.I - 1))
            @test abs(f2.EI / o.EnergyFlux.I - 1) <= tolI
            @test abs(f2.EH / o.EnergyFlux.H - 1) <= tolH
            @test f2.devI <= dtol                          # steady state reached
        end
        @test errs[2] < errs[1]                            # ladder converges onto the oracle
    end

    @testset "gate K5b: s=-2 (3,3) a=0.9 r0=10 at (8,48)" begin
        # MEASURED (c03.dat; band ℓ ∈ [3,5], tf=1000): ĵ=3 EI rel −2.40e-2, EH rel
        # +4.93e-5, devI 5.1e-4 (and even the subdominant ĵ=4 sits at −2.43e-2 / +4.9e-4)
        # → gates 6e-2 (2.5x) / 5e-4 (10x).
        o = TeukolskyPointParticleModeCircular(-2, 3, 3, 0.9, 10.0)
        fx = _kerr_flux_case(-2, 3, 0.9, 10.0, 5; N=8, K=48)
        f3 = fx[1]
        @test abs(f3.EI / o.EnergyFlux.I - 1) <= 6e-2
        @test abs(f3.EH / o.EnergyFlux.H - 1) <= 5e-4
        @test f3.devI <= 0.05
    end

    # a = 0.5 runs shared by gates K5c and K6 (computed once in the enclosing testset body)
    fx6_a05 = _kerr_flux_case(-2, 2, 0.5, 10.0, 4; N=6, K=32)
    fx8_a05 = _kerr_flux_case(-2, 2, 0.5, 10.0, 4; N=8, K=48)

    @testset "gate K5c: s=-2 (2,2) a=0.5 r0=10 two-level convergence" begin
        # MEASURED (c05,c04.dat; band ℓ ∈ [2,4], tf=1000): (6,32) EI rel +5.06e-2;
        # (8,48) EI rel −3.47e-3, EH rel +5.53e-4, devI 6.1e-4 → gates 0.12 (2.4x) /
        # 1e-2 (2.9x) / 2.5e-3 (4.5x); two-level convergence asserted.
        o = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.5, 10.0)
        @test abs(fx6_a05[1].EI / o.EnergyFlux.I - 1) <= 0.12
        @test abs(fx8_a05[1].EI / o.EnergyFlux.I - 1) <= 1e-2
        @test abs(fx8_a05[1].EH / o.EnergyFlux.H - 1) <= 2.5e-3
        @test abs(fx8_a05[1].EI / o.EnergyFlux.I - 1) < abs(fx6_a05[1].EI / o.EnergyFlux.I - 1)
    end

    # ---- (K6) SKH gr-qc/0703028 published-table anchor --------------------------------------------
    @testset "gate K6: SKH [gr-qc/0703028] Table FDm2/FDm3 anchors" begin
        # FD-oracle |m|-mode sums vs the published FD column (Hughes h2000 code; 4
        # published digits). MEASURED rel (dgs2_skh_anchor.out): 1.36e-5 (a=0.9, |m|=2),
        # 6.44e-5 (a=0.5, |m|=2), 2.36e-5 (|m|=3) — i.e. agreement to ALL published digits
        # → anchor gates at 0.5% (published-accuracy class, NOT a tight gate). Cited rows:
        # Table I (label FDm2) (a=0.9, r0=10) Ė_FD=4.456e-5, (a=0.5, r0=10) 4.812e-5;
        # Table II (FDm3) (a=0.9, r0=10) 5.044e-6.
        sum09 = 2 * sum(TeukolskyPointParticleModeCircular(-2, l, 2, 0.9, 10.0).EnergyFlux.I
                        for l in 2:8)
        sum05 = 2 * sum(TeukolskyPointParticleModeCircular(-2, l, 2, 0.5, 10.0).EnergyFlux.I
                        for l in 2:8)
        sum09m3 = 2 * sum(TeukolskyPointParticleModeCircular(-2, l, 3, 0.9, 10.0).EnergyFlux.I
                          for l in 3:9)
        @test abs(sum09 / 4.456e-5 - 1) <= 5e-3
        @test abs(sum05 / 4.812e-5 - 1) <= 5e-3
        @test abs(sum09m3 / 5.044e-6 - 1) <= 5e-3
        # TD side (chain anchor): the a=0.5 (8,48) band total from the K5c run vs the SKH
        # row — MEASURED rel −3.43e-3 → gate 1.5e-2 (4.4x; TD discretization + band
        # truncation; the ℓ ≥ 5 tail is 2e-9 of the |m|=2 sum, dgs2_skh_anchor.out).
        td05 = 2 * sum(f.EI for f in fx8_a05)
        @test abs(td05 / 4.812e-5 - 1) <= 1.5e-2
    end
end
