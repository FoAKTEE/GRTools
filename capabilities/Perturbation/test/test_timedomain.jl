# DG-P3: Perturbation.TimeDomain characterization vs the ORIGINAL DG-Teulolsky repo.
# Grounding: baked one-shot fixtures test/fixtures/td_*.jl (captured by
# scratchpad/capture_td_dg.jl running /data/haiyangw/claude/DG-Teulolsky in ITS OWN env at the
# exact configs recorded in each fixture header) + the repo's own unit testsets
# (Hyperboloidal / Backgrounds / sine-Gaussian / Mode1D short check / FlatWave convergence)
# ported verbatim. Tolerances are MEASURED port-vs-original deviations with ~10x headroom
# (dominant divergence source: NumericsCore's self-contained LGL rule vs FastGaussQuadrature at
# ~1e-15 node level + Spacetime's tighter inverse-tortoise stopping (8·eps vs 1e-13), amplified
# through 3e4–7e4 RK4 steps; measured numbers quoted at each testset).
using Test
using Perturbation
using Perturbation.TimeDomain
const TD = Perturbation.TimeDomain

include(joinpath(@__DIR__, "fixtures", "td_substitution.jl"))
include(joinpath(@__DIR__, "fixtures", "td_mode1d.jl"))
include(joinpath(@__DIR__, "fixtures", "td_flux1d.jl"))
include(joinpath(@__DIR__, "fixtures", "td_selfforce.jl"))
include(joinpath(@__DIR__, "fixtures", "td_kerrtail.jl"))
include(joinpath(@__DIR__, "fixtures", "td_kerrflux.jl"))
include(joinpath(@__DIR__, "flatwave_support.jl"))
using .TDFlatWave

# rel-vs-baked helper: |a-b|/|b| with the b=0 columns compared absolutely
_rel(a, b) = b == 0 ? abs(a - b) : abs(a - b) / abs(b)

# MEASURED port-vs-original deviations (scratchpad/measure_td_port.jl, 2026-07-02) and the
# tolerances derived with ~10-100x headroom:
#   substitutions: Ylm 3.9e-15, orbital_params ≤2.2e-15, tortoise 0 (bit-exact),
#                  inverse_tortoise 0 (bit-exact), potential/coupling/friction/hyperboloidal 0.
#   constructions: ALL probed C/H/A/S⁻¹/V/Minj/μ/Gpre/probe values 0 to ≤2.7e-16 rel.
#   evolutions   : flux1d ≤1.4e-12 rel, selfforce ≤4.8e-14 rel, kerrflux ≤7.1e-14 rel;
#                  decaying tails (mode1d/kerrtail) ≤7.2e-14 / ≤1.5e-14 of the run's PEAK
#                  amplitude (pure relative error is meaningless on samples decayed 10-70
#                  orders below peak; the ulp-level LGL difference NumericsCore-vs-
#                  FastGaussQuadrature is the irreducible noise floor there).
#   tail exponents: mode1d ≤1.6e-8; kerrtail n_l0 1.7e-8, n_l2 5.1e-3, n_l4 1.6e-3 (the l2/l4
#                  fit windows are dominated by sub-peak·1e-12 samples — same noise floor).
const TOL_MODE1D_PEAK  = 5e-12   # |Δψ| / max|ψ(run)|      (measured ≤ 7.2e-14)
const TOL_KERRTAIL_PEAK = 1e-12  # |Δψ| / max|ψ_mode(run)| (measured ≤ 1.5e-14)
const TOL_FLUX1D       = 2e-11   # rel                     (measured ≤ 1.4e-12)
const TOL_SELFFORCE    = 5e-13   # rel                     (measured ≤ 4.8e-14)
const TOL_KERRFLUX     = 1e-12   # rel                     (measured ≤ 7.1e-14)
const TOL_KERR_CONS    = 1e-13   # rel, construction probes (measured 0, bit-exact)
const TOL_TAIL_EXP     = 1e-6    # |Δn|, clean-window fits  (measured ≤ 1.6e-8)
const TOL_TAIL_EXP_NOISY = 5e-2  # |Δn|, noise-floor-dominated l2/l4 windows (measured ≤ 5.1e-3)

@testset "TimeDomain (DG-P3 s=0 time-domain DG Teukolsky port)" begin

    # ---- repo's own unit testsets, ported verbatim ---------------------------------------------
    @testset "Hyperboloidal layer (repo testset)" begin
        L = TD.Layer(30.0, 50.0, 4)
        @test TD.Hfun(L, 10.0) == 0.0          # below the layer
        @test TD.Ωfun(L, 10.0) == 1.0
        @test TD.Hfun(L, 50.0) ≈ 1.0 atol = 1e-12   # at scri
        @test TD.Ωfun(L, 50.0) ≈ 0.0 atol = 1e-12
        @test TD.denomΩ(L, 50.0) ≈ 50 * 4 / 20 atol = 1e-10  # s·P/(s-R)
        # H continuous & = 0 at the layer start (P=4 ⇒ Ω'(R)=0)
        @test TD.Hfun(L, 30.0) ≈ 0.0 atol = 1e-12
        # baked repo values (verbatim arithmetic ⇒ bit-exact; measured 0.0)
        F = TD_HYP_TEST_REF
        @test TD.Hfun(L, 10.0) == F.H10 && TD.Ωfun(L, 10.0) == F.Om10
        @test TD.Hfun(L, 50.0) == F.H50 && TD.Ωfun(L, 50.0) == F.Om50
        @test TD.denomΩ(L, 50.0) == F.dOm50 && TD.Hfun(L, 30.0) == F.H30
        for (R, s, P, ρ, Om, dOm, den, H) in TD_HYP_REF
            LL = TD.Layer(R, s, P)
            @test TD.Ωfun(LL, ρ) == Om
            @test TD.Ωprime(LL, ρ) == dOm
            @test TD.denomΩ(LL, ρ) == den
            @test TD.Hfun(LL, ρ) == H
            @test TD.omega_factor(LL, ρ) == 1 - H
        end
    end

    @testset "Backgrounds (repo testset, Spacetime-delegated tortoise)" begin
        # a=0 potential == standard Schwarzschild scalar potential
        for (r, l) in ((6.0, 2), (10.0, 0), (20.0, 3))
            @test TD.potential_V(r, l, 0; M=1.0, a=0.0) ≈
                  -(1 - 2 / r) * (l * (l + 1) / r^2 + 2 / r^3) rtol = 1e-12
        end
        # tortoise/inverse round-trip (Schwarzschild + Kerr) through the Spacetime shims
        for r in (2.5, 6.0, 100.0, 1000.0)
            @test TD._inverse_tortoise(TD._tortoise(r; M=1.0, a=0.0); M=1.0, a=0.0) ≈ r rtol = 1e-10
        end
        for r in (3.0, 50.0), a in (0.5, 0.9)
            @test TD._inverse_tortoise(TD._tortoise(r; M=1.0, a=a); M=1.0, a=a) ≈ r rtol = 1e-9
        end
        # coupling/friction vanish at a=0
        @test TD.coupling_f(10.0; a=0.0) == 0.0
        @test TD.friction_mu(10.0, 2; a=0.0) == 0.0
        # baked repo values (local potential_V/coupling_f/friction_mu are verbatim ⇒ bit-exact)
        for (a, r, l, m, V, f, mui) in TD_BG_REF
            @test TD.potential_V(r, l, m; M=1.0, a=a) == V
            @test TD.coupling_f(r; M=1.0, a=a) == f
            @test imag(TD.friction_mu(r, m; M=1.0, a=a)) == mui
        end
    end

    # ---- EXTRA GATE: reuse-substitution equivalence vs baked repo values ------------------------
    @testset "substitution gate: tortoise/inverse ≡ Spacetime (baked repo)" begin
        # tortoise: identical closed form through Spacetime.tortoise_kerr ⇒ measured 0.0
        # (bit-exact); assert ≤ 1e-15 rel (a few ulps of robustness)
        for (a, r, v) in TD_TORTOISE_REF
            @test _rel(TD._tortoise(r; M=1.0, a=a), v) <= 1e-15
        end
        # inverse: local verbatim repo iteration (deep-horizon floor is load-bearing — see
        # TimeDomain/background.jl) evaluating r* through Spacetime ⇒ measured 0.0 (bit-exact)
        for (a, rs, v) in TD_INVTORT_REF
            @test _rel(TD._inverse_tortoise(rs; M=1.0, a=a), v) <= 1e-15
        end
    end

    @testset "substitution gate: Ylm_equator ≡ Harmonics (baked repo, 18 (l,m))" begin
        # measured ≤ 3.9e-15 rel (the repo normalizes via exp(-Σlog k), Harmonics via the exact
        # reciprocal falling factorial — agreement well inside the ≤1e-12-class gate)
        for (l, m, v) in TD_YLM_REF
            if v == 0.0
                @test TD.Ylm_equator(l, m) == 0.0        # exact parity zero (odd l+|m|)
            else
                @test _rel(TD.Ylm_equator(l, m), v) <= 1e-13
            end
        end
    end

    @testset "substitution gate: orbital_params ≡ Spacetime (baked repo, a∈{0,0.5,0.9})" begin
        # measured ≤ 2.2e-15 rel across E/Lz/u^t/Ω, 12 (a, rp) points
        for (a, rp, E, Lz, ut, Om) in TD_ORBIT_REF
            op = TD.orbital_params(rp; M=1.0, a=a)
            @test _rel(op.E, E) <= 1e-13
            @test _rel(op.Lz, Lz) <= 1e-13
            @test _rel(op.ut, ut) <= 1e-13
            @test _rel(op.Ω, Om) <= 1e-13
        end
    end

    @testset "sine-Gaussian profile derivatives vs finite differences (repo testset)" begin
        h = 1e-5
        for u in (-12.0, -10.0, -8.0, -3.0, 0.0)
            f, f1, f2, f3 = profile_derivs(u)
            fp = profile_derivs(u + h); fm = profile_derivs(u - h)
            @test (fp[1] - fm[1]) / 2h ≈ f1 atol = 1e-5
            @test (fp[2] - fm[2]) / 2h ≈ f2 atol = 1e-5
            @test (fp[3] - fm[3]) / 2h ≈ f3 atol = 1e-4
        end
    end

    # ---- (a) Mode1D Schwarzschild Price tail (baked characterization + repo short check) --------
    @testset "characterization (a): Mode1D Schwarzschild ℓ=0 tail" begin
        F = TD_MODE1D
        p = TD.Mode1DProblem(0; N=8, K=120, ρL=-200.0, R=150.0, s=1200.0, P=4,
                             extract_rstar=[500.0, Inf])
        # construction probes (inverse_tortoise→Spacetime + verbatim C build); measured ≤ 2e-16
        @test p.probe_k == F.probe_k
        @test maximum(_rel.(p.probe_ρ, F.probe_rho)) <= 1e-13
        @test _rel(p.C[3, 17], F.C_3_17) <= 1e-13
        @test _rel(p.C[5, 60], F.C_5_60) <= 1e-13
        @test _rel(p.C[end, end], F.C_scri) <= 1e-14
        @test _rel(p.H[5, 100], F.H_5_100) <= 1e-14
        τs, Ψ = TD.solve_mode1d(p; tfinal=6000.0, dt=0.0866, record_every=100)
        @test length(τs) == F.nsamp
        pk1 = maximum(abs.(Ψ[1, :])); pk2 = maximum(abs.(Ψ[2, :]))
        for (i, j) in enumerate(F.idx)
            @test τs[j] == F.tau[i]
            @test abs(Ψ[1, j] - F.psi_r500[i]) <= TOL_MODE1D_PEAK * pk1
            @test abs(Ψ[2, j] - F.psi_scri[i]) <= TOL_MODE1D_PEAK * pk2
        end
        n1 = TD.fit_tail(τs, Ψ[1, :]; window=(2500.0, 6000.0))
        n2 = TD.fit_tail(τs, Ψ[2, :]; window=(2500.0, 6000.0))
        @test abs(n1 - F.n_r500) <= TOL_TAIL_EXP
        @test abs(n2 - F.n_scri) <= TOL_TAIL_EXP
        # repo's own physics assertion (runtests "Mode1D Schwarzschild tail"): n ≈ -2ℓ-3 = -3
        @test isfinite(n1)
        @test abs(n1 - (-3.0)) < 0.4
    end

    # ---- (b) Flux1D circular-orbit flux + self-force ---------------------------------------------
    @testset "characterization (b): Flux1D (2,2) rp=10 flux + selfforce_ft" begin
        F = TD_FLUX1D
        p = TD.Flux1DProblem(2, 2, 10.0; N=8, K=100, ρL=-150.0, R=100.0, s=200.0, P=4,
                             ramp_τ=400.0, ramp_δ=0.00025)
        # construction probes (tortoise/orbital_params→Spacetime, Ylm→Harmonics); measured ≤ 3e-16
        @test p.kp == F.kp
        @test _rel(p.rsp, F.rsp) <= 1e-14
        @test _rel(real(p.Gpre), real(F.Gpre)) <= 1e-13
        @test abs(imag(p.Gpre) - imag(F.Gpre)) <= 1e-15
        @test _rel(p.Ω, F.Omega) <= 1e-13
        @test _rel(p.C[4, 37], F.C_4_37) <= 1e-13
        @test _rel(p.C[end, end], F.C_scri) <= 1e-14
        @test _rel(TD.orbital_params(10.0; M=1.0, a=0.0).ut, F.ut) <= 1e-13
        τs, Es, Eh = TD.solve_flux(p; tfinal=1500.0, dt=0.01, record_every=200)
        @test length(τs) == F.nsamp
        for (i, j) in enumerate(F.idx)
            @test τs[j] == F.tau[i]
            @test _rel(Es[j], F.Escri[i]) <= TOL_FLUX1D
            @test _rel(Eh[j], F.Ehor[i]) <= TOL_FLUX1D
        end
        @test _rel(Es[end], F.Escri_end) <= TOL_FLUX1D
        @test _rel(Eh[end], F.Ehor_end) <= TOL_FLUX1D
        # physics sanity vs the frequency-domain BHPToolkit gold of the repo's T2.1 table (the
        # baked shrunk-config value itself sits ~2e-9 rel from gold — assert the same ballpark)
        @test abs(Es[end] - 3.3699774706034e-6) / 3.3699774706034e-6 < 1e-4
        ftL, ftR = TD.selfforce_ft(p; tfinal=1200.0, dt=0.01, twin=800.0)
        @test _rel(ftL, TD_SELFFORCE.ftL) <= TOL_SELFFORCE
        @test _rel(ftR, TD_SELFFORCE.ftR) <= TOL_SELFFORCE
    end

    # ---- (c) KerrModes m=0 coupled tail -----------------------------------------------------------
    @testset "characterization (c): KerrModes a=0.9 m=0 ℓs=[0,2,4] tail" begin
        F = TD_KERRTAIL
        p = TD.KerrTailProblem([0, 2, 4], 0; N=8, K=120, ρL=-400.0, R=149.35, s=400.0, P=4,
                               a=0.9, extract_rstar=[100.0])
        # mode-coupling closed forms (verbatim ⇒ bit-exact)
        @test TD.cminus(2, 0) == F.cminus_20
        @test collect(TD.mode_coupling(2, 0)) == F.coup_l2
        @test collect(TD.mode_coupling(4, 0)) == F.coup_l4
        @test TD.mode_list(0, 4; parity=:even) == [0, 2, 4]
        # construction probes (coupled A/S⁻¹/V via Spacetime tortoise + unified _Amat)
        @test maximum(_rel.(p.probe_ρ, F.probe_rho)) <= 1e-13
        @test maximum(_rel.(vec(p.Aflux[:, :, 5, 40]), F.Aflux_5_40)) <= TOL_KERR_CONS
        @test maximum(abs.(vec(p.Aflux[:, :, end, end]) .- F.Aflux_scri)) == 0  # exact ½-blocks
        @test maximum(_rel.(vec(p.Sinv[:, :, 5, 40]), F.Sinv_5_40)) <= TOL_KERR_CONS
        @test maximum(_rel.(p.Vmode[:, 5, 40], F.Vmode_5_40)) <= TOL_KERR_CONS
        τs, Ψ = TD.solve_kerr_tail(p; init_mode=2, tfinal=3000.0, dt=0.04, record_every=100,
                                   ic=(amp=1 / sqrt(2π * 6^2), ρ0=25.0, w=6.0))
        @test length(τs) == F.nsamp
        pk = [maximum(abs.(Ψ[1, q, :])) for q in 1:3]
        for (i, j) in enumerate(F.idx)
            @test τs[j] == F.tau[i]
            @test abs(Ψ[1, 1, j] - F.psi_l0[i]) <= TOL_KERRTAIL_PEAK * pk[1]
            @test abs(Ψ[1, 2, j] - F.psi_l2[i]) <= TOL_KERRTAIL_PEAK * pk[2]
            @test abs(Ψ[1, 3, j] - F.psi_l4[i]) <= TOL_KERRTAIL_PEAK * pk[3]
        end
        ns = [TD.fit_tail_window(τs, Ψ[1, q, :], (1200.0, 3000.0)) for q in 1:3]
        @test abs(ns[1] - F.n_l0) <= TOL_TAIL_EXP
        @test abs(ns[2] - F.n_l2) <= TOL_TAIL_EXP_NOISY   # window ψ_l2 ~ peak·1e-12 (noise floor)
        @test abs(ns[3] - F.n_l4) <= TOL_TAIL_EXP_NOISY   # window ψ_l4 ~ peak·1e-13 (noise floor)
        # physics: the (2,0)-seeded l=0 tail along r=const decays as τ^{-(l+l'+1)} = τ^{-3}
        @test abs(ns[1] - (-3.0)) < 0.4
    end

    # ---- (d) KerrFlux coupled circular-orbit flux --------------------------------------------------
    @testset "characterization (d): KerrFlux a=0.9 ℓs=[1,3] m=1 Ė¹¹" begin
        F = TD_KERRFLUX
        rp = 16.0914363989845
        p = TD.KerrFluxProblem([1, 3], 1, rp; N=8, K=120, ρL=-300.0, R=100.0, s=300.0, P=4,
                               a=0.9, ramp_τ=200.0, ramp_δ=0.001)
        op = TD.orbital_params(rp; M=1.0, a=0.9)
        @test p.kp == F.kp
        @test _rel(p.Ω, F.Omega) <= 1e-13
        @test _rel(op.E, F.E) <= 1e-13
        @test _rel(op.Lz, F.Lz) <= 1e-13
        @test _rel(op.ut, F.ut) <= 1e-13
        @test maximum(_rel.(real.(p.Gpre), F.Gpre_re)) <= 1e-13
        @test maximum(abs.(imag.(p.Gpre) .- F.Gpre_im)) == 0
        @test maximum(abs.(vec(p.Minjp) .- F.Minjp)) <= TOL_KERR_CONS
        @test maximum(abs.(vec(p.Minjm) .- F.Minjm)) <= TOL_KERR_CONS
        @test _rel(imag(p.μmode[1, 5, 40]), F.mu_5_40_im) <= TOL_KERR_CONS
        τs, Es = TD.solve_kerr_flux(p; tfinal=1400.0, dt=0.02, record_every=200, target=1)
        @test length(τs) == F.nsamp
        for (i, j) in enumerate(F.idx)
            @test τs[j] == F.tau[i]
            @test _rel(Es[j], F.Escri[i]) <= TOL_KERRFLUX
        end
        @test _rel(Es[end], F.Escri_end) <= TOL_KERRFLUX
        # physics sanity: the repo's T1.3 [1,3] gold is 1.835753e-6 (full config); the shrunk
        # config sits ~0.5% below it — assert the same ballpark
        @test abs(Es[end] - 1.835753e-6) / 1.835753e-6 < 0.02
    end

    # ---- (e) FlatWave verification harness (repo testset, verbatim) -------------------------------
    @testset "flat-wave hyperboloidal solver converges to exact (repo testset)" begin
        scri(K) = let
            p = FlatWaveProblem(; N=4, K=K, R=30.0, s=50.0, P=4)
            τs, Ψ, _ = solve_flatwave(p; tfinal=50.0, dt=1e-3, record_every=40)
            rel_L2_in_time(p, τs, Ψ, 3)       # paper's E_{N,K} at 𝓘⁺ (ρ=50)
        end
        e32, e64 = scri(32), scri(64)
        @test isfinite(e32) && isfinite(e64)
        @test e64 < e32                       # convergence under refinement
        @test e64 < 1e-3                      # N=4,K=64 already well-resolved at 𝓘⁺
    end
end
