# DG-P5: SOURCED s = ±2 / s = 0 circular-equatorial point-particle sector at a = 0 — gates.
#
# Physics: the s = −2 (and +2) point-particle δ/δ′/δ″ source enters the RT-chart DG
# evolution through analytically-known jump conditions at the particle interface
# (spinsource.jl provenance header (1)–(6)); the steady-state scri/horizon π amplitudes
# give Ė^∞/Ė^H, gated against the in-package FREQUENCY-DOMAIN oracle
# `TeukolskyPointParticleModeCircular` (AmplitudesFluxes: MST radial + GHP source
# convolution — an independent solution pipeline sharing only the source-projection
# coefficients).
#
# GATES (tolerances MEASURED on the exact configs below; headroom stated inline):
#   (J1) chart-transform derivatives : analytic F′(R), r(R)↔R closed forms.
#   (J2) jump algebra == Barack–Giudice: generic distributional jump algebra + GHP source
#        reproduces the INDEPENDENT closed forms [arXiv:1702.04204, Eqs. (A15)–(A17)]
#        (s = +2; different papers, different construction — validates the algebra AND
#        every Δ-weight/normalization convention).
#   (J3) jumps == FD Green's function : [X], [X′] == Z^∞R̃_up − Z^H R̃_in at r₀ (unit-
#        transmission Radial solutions + oracle amplitudes; validates the Green's-fn
#        representation the steady state must reproduce).
#   (S0) s = 0 reduction IDENTITY     : the generic jump formulas at s = 0 equal the
#        locked Flux1D δ-injection data ([ψ̂] = 0, f₀[ψ̂′_r] = Gpre) in closed form.
#   (F−2) s = −2 flux vs FD oracle    : resolution ladder (N,K) = (6,32)→(7,40)→(8,48)
#        →(10,64) for (2,2), two-level convergence for (3,3), converged-level (2,1);
#        all at r₀ = 10M; Ė^∞ and Ė^H gated at measured values with ≥2x headroom.
#   (F0) s = 0 sourced regression     : the SAME sourced machinery at s = 0 reproduces the
#        locked layer-chart Flux1D flux (baked fixture td_flux1d.jl) and the FD oracle.
#   (F+2) s = +2 flux vs FD oracle    : the numerically fragile spin sign (vacuum evidence
#        in test_timedomain_spin.jl); sourced run is signal-dominated and gates green.
#
# All runs a = 0, M = 1 (Kerr a ≠ 0 sourced = documented [FUTURE], dossier §8.3).

using Test
using Perturbation
using Perturbation: TeukolskyPointParticleModeCircular, TeukolskyRadial, teukolsky_radial_deriv
using Harmonics: SpinWeightedSpheroidalEigenvalue
const TDS5 = Perturbation.TimeDomain

@testset "TimeDomain spin source (DG-P5)" begin

    # ---- (J1) chart maps and derivatives ---------------------------------------------------------
    @testset "gate J1: RT chart maps for the jump push-forward" begin
        for r in (6.0, 8.0, 10.0, 20.0)
            @test abs(TDS5.rt_r_of_R(TDS5.rt_R_of_r(r)) - r) <= 1e-12 * r     # inverse pair
        end
        # analytic F′(R) vs 4th-order central difference (MEASURED rel 1.5e-10, FD-limited)
        for r in (8.0, 10.0)
            R = TDS5.rt_R_of_r(r)
            h = 1e-6
            fd = (8 * (TDS5.rt_F_of_R(R + h) - TDS5.rt_F_of_R(R - h)) -
                  (TDS5.rt_F_of_R(R + 2h) - TDS5.rt_F_of_R(R - 2h))) / (12h)
            @test abs(TDS5._rt_dF_dR(R) - fd) <= 1e-7 * abs(fd)
        end
        # SWSH eigenvalue identity used by the ODE coefficients and the flux layer
        for (s, l, m) in ((-2, 2, 2), (-2, 3, 3), (2, 2, 2), (0, 2, 2))
            @test SpinWeightedSpheroidalEigenvalue(s, l, m, 0.0) ≈ l * (l + 1) - s * (s + 1) atol = 1e-12
        end
    end

    # ---- (J2) generic jump algebra == Barack–Giudice closed forms (s = +2) -----------------------
    @testset "gate J2: jump algebra vs BG closed forms [1702.04204 (A15)-(A17)]" begin
        # MEASURED max rel across the four configs: 5e-16 (pure closed-form arithmetic on
        # both sides; SWSH values shared) → gate 5e-14 (100x headroom).
        for (l, m, r0) in ((2, 2, 10.0), (3, 3, 8.0), (2, 1, 10.0), (4, 2, 7.0))
            J = TDS5.spin_source_jumps(2, l, m, r0)
            B = TDS5.bg_jumps_s2(l, m, r0)
            @test abs(J.psid - B.psid) <= 5e-14 * abs(B.psid)
            @test abs(J.Jpsi - B.Jpsi) <= 5e-14 * abs(B.Jpsi)
            @test abs(J.Jpsip - B.Jpsip) <= 5e-14 * abs(B.Jpsip)
        end
    end

    # ---- (J3) jumps == FD Green's-function representation ----------------------------------------
    @testset "gate J3: jumps vs Z^∞R_up − Z^H R_in (FD construction)" begin
        # X(r>r₀) = Z^∞R̃_up, X(r<r₀) = Z^H R̃_in (unit transmission; spinsource.jl header
        # (6)). MEASURED rel: 1.9e-14 (s=−2, 2,2), 4.8e-13 (s=−2, 3,3), 6.6e-15 (s=+2)
        # — MST-series-limited → gate 1e-10 (≥200x headroom).
        for (s, l, m, r0) in ((-2, 2, 2, 10.0), (-2, 3, 3, 8.0), (2, 2, 2, 10.0))
            J = TDS5.spin_source_jumps(s, l, m, r0)
            o = TeukolskyPointParticleModeCircular(s, l, m, 0.0, r0)
            R = TeukolskyRadial(s, l, m, 0.0, J.omega; lambda=o.lambda)
            JX_FD = o.Zinf * R.Up(r0) - o.ZH * R.In(r0)
            JXp_FD = o.Zinf * teukolsky_radial_deriv(R.Up, r0) -
                     o.ZH * teukolsky_radial_deriv(R.In, r0)
            @test abs(J.JX - JX_FD) <= 1e-10 * abs(JX_FD)
            @test abs(J.JXp - JXp_FD) <= 1e-10 * abs(JXp_FD)
        end
    end

    # ---- (S0) s = 0 reduction identity vs the locked Flux1D injection ----------------------------
    @testset "gate S0: s=0 jump formulas == Flux1D δ-injection (closed form)" begin
        # Flux1D continuum source (flux1d.jl): π̇ = −∂_{r*}χ − Vψ + Gπδ(r*−r*_p) ⇒
        # [∂_{r*}ψ] = Gpre·e^{−imΩt}, [ψ] = 0, with Gpre = −4πq/(u^t r_p)·Y_{ℓm}(π/2,0).
        # The generic machinery must give Jpsi ≡ 0 and f₀·Jpsip ≡ Gpre EXACTLY
        # (MEASURED rel ≤ 3e-15 — same closed rationals) → gate 1e-12.
        for (l, m, r0) in ((2, 2, 10.0), (3, 3, 8.0), (4, 4, 10.0))
            J = TDS5.spin_source_jumps(0, l, m, r0)
            op = TDS5.orbital_params(r0)
            Gpre = -4pi / (op.ut * r0) * TDS5.Ylm_equator(l, m)
            f0 = 1 - 2 / r0
            @test abs(J.Jpsi) == 0.0                        # scalar field continuous
            @test abs(J.psid) == 0.0                        # no δ-layer at s = 0
            @test abs(f0 * J.Jpsip - Gpre) <= 1e-12 * abs(Gpre)
        end
    end

    # ---- constructor guards -----------------------------------------------------------------------
    @testset "sourced constructor guards (a=0 scope, honest [FUTURE] boundary)" begin
        @test_throws ArgumentError TDS5.SpinSourcedProblem(-1, [2], 2, 10.0)   # |s|=1 [FUTURE]
        @test_throws ArgumentError TDS5.spin_source_jumps(-2, 2, 0, 10.0)      # m=0 fluxless
        p = TDS5.SpinSourcedProblem(-2, [2], 2, 10.0; N=4, K=16)
        R0 = TDS5.rt_R_of_r(10.0)
        @test any(abs.(p.base.mesh.edges .- R0) .<= 1e-14)                     # forced edge
        @test p.base.mesh.edges[p.kp+1] ≈ R0 atol = 1e-14                      # kp addresses it
    end

    # ---- shared flux driver ------------------------------------------------------------------------
    function _flux_case(s, l, m, r0; N, K, fa, tf=1000.0, twin=(0.7, 1.0))
        p = TDS5.SpinSourcedProblem(s, [l], m, r0; N=N, K=K, filter_alpha=fa)
        dt = TDS5.spin_cfl_dt(p.base; C=0.6)
        τs, Ψ, Π, maxu = TDS5.solve_spin_sourced(p; tfinal=tf, dt=dt,
                                                 record_every=max(1, round(Int, 5.0 / dt)))
        @test all(isfinite, maxu)                    # boundedness on every run
        return TDS5.spin_sourced_flux(p, τs, Π; twin=twin)[1], maxu
    end

    # ---- (F−2) s = −2 flux ladder + modes vs FD oracle --------------------------------------------
    @testset "gate F-2: s=-2 flux ladder (2,2) r0=10 vs FD oracle" begin
        # MEASURED (tf=1000, twin=(0.7,1); scratchpad/dgp5_measure.jl log):
        #   (N,K)=(6,32):  EI rel −1.131e-1, EH rel −8.86e-3, devI 9.0e-3
        #   (N,K)=(7,40):  EI rel −1.093e-1, EH rel +1.58e-3, devI 3.2e-3
        #   (N,K)=(8,48):  EI rel +6.92e-3, EH rel +6.53e-4, devI 5.0e-4
        #   (N,K)=(10,64): EI rel +1.34e-3, EH rel −4.90e-5, devI 1.2e-4  (tf=1000);
        #                  tf=800 (suite config): EI rel +1.32e-3, EH rel −5.35e-5
        # → the flux CONVERGES through the ladder onto the FD oracle; gates at each level
        #   with ≥2x headroom; converged-level gate 5e-3 on Ė^∞ (3.8x) / 5e-4 on Ė^H (9x).
        o = TeukolskyPointParticleModeCircular(-2, 2, 2, 0.0, 10.0)
        errs = Float64[]
        for (N, K, tolI, tolH, tf, twin) in (
            (6, 32, 0.25, 0.02, 1000.0, (0.7, 1.0)),
            (7, 40, 0.25, 0.01, 1000.0, (0.7, 1.0)),
            (8, 48, 0.02, 5e-3, 1000.0, (0.7, 1.0)),
            (10, 64, 5e-3, 5e-4, 800.0, (0.68, 1.0)),
        )
            fx, _ = _flux_case(-2, 2, 2, 10.0; N=N, K=K, fa=36.0, tf=tf, twin=twin)
            push!(errs, abs(fx.EI / o.EnergyFlux.I - 1))
            @test abs(fx.EI / o.EnergyFlux.I - 1) <= tolI
            @test abs(fx.EH / o.EnergyFlux.H - 1) <= tolH
            @test fx.devI <= 0.05                       # steady state reached
        end
        @test errs[end] < errs[1]                        # ladder converges toward the oracle
        @test errs[end] < errs[3]                        # ... and beyond the (8,48) level
    end

    @testset "gate F-2: s=-2 mode (3,3) r0=10 two-level convergence" begin
        # MEASURED: (8,48) tf=1000: EI rel +4.54e-2, EH rel +6.64e-5 (small-|Z^∞| mode:
        # the scri signal |π| ≈ 5.7e-6 sits closer to the filter-debris floor → positive
        # power bias at the coarse level); (10,64) tf=1000: EI rel −2.80e-3, EH +6.88e-6;
        # (10,64) tf=800 (suite config): EI rel −2.77e-3, EH +6.24e-6, devI 4.7e-4.
        # Gates: (8,48) 0.09/1e-3 (2x/15x), (10,64) 1e-2/1e-4 (3.6x/16x), and the top
        # level must beat the coarse level (convergence onto the oracle).
        o = TeukolskyPointParticleModeCircular(-2, 3, 3, 0.0, 10.0)
        fx8, _ = _flux_case(-2, 3, 3, 10.0; N=8, K=48, fa=36.0)
        @test abs(fx8.EI / o.EnergyFlux.I - 1) <= 0.09
        @test abs(fx8.EH / o.EnergyFlux.H - 1) <= 1e-3
        fx10, _ = _flux_case(-2, 3, 3, 10.0; N=10, K=64, fa=36.0, tf=800.0, twin=(0.68, 1.0))
        @test abs(fx10.EI / o.EnergyFlux.I - 1) <= 1e-2
        @test abs(fx10.EH / o.EnergyFlux.H - 1) <= 1e-4
        @test abs(fx10.EI / o.EnergyFlux.I - 1) < abs(fx8.EI / o.EnergyFlux.I - 1)
    end

    @testset "gate F-2: s=-2 mode (2,1) r0=10 at (10,64) (small-|Z| mode)" begin
        # The (2,1) flux is 9.66e-8 (|Z^∞| ≈ 6e-5, scri |π| ≈ 1.9e-6) — at (8,48) the TD
        # scri signal is floor-contaminated (MEASURED: EI rel +5.70e-1 with devI 6.3e-2
        # at tf=1000; +1.12 at tf=1600; +3.34e-1 with a permanent ψ-penalty — a numerical
        # floor, NOT a settling issue). At (10,64), tf=1400 the floor drops below the
        # signal and the mode CONVERGES: MEASURED EI rel +1.52e-2 (devI 1.5e-3),
        # EH rel −4.36e-6 → gates 4e-2 (2.6x) / 1e-4 (23x).
        o = TeukolskyPointParticleModeCircular(-2, 2, 1, 0.0, 10.0)
        fx, _ = _flux_case(-2, 2, 1, 10.0; N=10, K=64, fa=36.0, tf=1400.0)
        @test abs(fx.EI / o.EnergyFlux.I - 1) <= 4e-2
        @test abs(fx.EH / o.EnergyFlux.H - 1) <= 1e-4
        @test fx.devI <= 2e-2
    end

    # ---- (F0) s = 0 sourced regression vs locked Flux1D -------------------------------------------
    @testset "gate F0: s=0 sourced machinery vs locked Flux1D flux (2,2) r0=10" begin
        # The SAME jump-enforced machinery at s = 0 (filter off — scalar sector is stable
        # unfiltered, cf. test_timedomain_spin.jl measured evidence) must reproduce the
        # LOCKED layer-chart Flux1D flux (baked fixture; solve_flux gated ≤ 1.4e-12 vs the
        # DG-Teulolsky repo) — a chart-independent physical-flux identity. MEASURED
        # (scratchpad/dgp5_measure.jl): (6,32): EI rel(F1D) −1.00e-8, EH rel +1.98e-7;
        # (8,48): EI rel −6.4e-10, EH rel +1.5e-10; vs FD oracle −7.0e-10.
        # Gates at 10–30x headroom.
        isdefined(Main, :TD_FLUX1D) || include(joinpath(@__DIR__, "fixtures", "td_flux1d.jl"))
        o0 = TeukolskyPointParticleModeCircular(0, 2, 2, 0.0, 10.0)
        fx6, _ = _flux_case(0, 2, 2, 10.0; N=6, K=32, fa=0.0)
        @test abs(fx6.EI / TD_FLUX1D.Escri_end - 1) <= 1e-7
        @test abs(fx6.EH / TD_FLUX1D.Ehor_end - 1) <= 5e-6
        fx8, _ = _flux_case(0, 2, 2, 10.0; N=8, K=48, fa=0.0)
        @test abs(fx8.EI / TD_FLUX1D.Escri_end - 1) <= 1e-8
        @test abs(fx8.EH / TD_FLUX1D.Ehor_end - 1) <= 5e-9
        @test abs(fx8.EI / o0.EnergyFlux.I - 1) <= 1e-8      # and the FD oracle agrees
    end

    # ---- (F+2) s = +2 flux vs FD oracle ------------------------------------------------------------
    @testset "gate F+2: s=+2 sourced flux (2,2) r0=10 (BG-jump-sourced)" begin
        # The numerically fragile spin sign (vacuum: scri-glancing growth, gate C+E
        # evidence) — but the sourced steady state is SIGNAL-dominated (jump scale ~4e2,
        # scri amplitude |Z^∞_{+2}| ~ 1.5e3) and the finite-window run is clean.
        # MEASURED (8,48), tf=1000: EI rel −3.83e-4 (devI 1.4e-7), EH rel −2.50e-2,
        # maxu 2.26e4 finite/steady → gates 2e-3 (5x) and 6e-2 (2.4x).
        o = TeukolskyPointParticleModeCircular(2, 2, 2, 0.0, 10.0)
        fx, maxu = _flux_case(2, 2, 2, 10.0; N=8, K=48, fa=36.0)
        @test abs(fx.EI / o.EnergyFlux.I - 1) <= 2e-3
        @test abs(fx.EH / o.EnergyFlux.H - 1) <= 6e-2
        @test fx.devI <= 1e-4                              # steady (measured 1.4e-7)
    end
end
