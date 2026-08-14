# SpinFluxes: the 10 TeukolskySpinMode* entrypoints ported in S5-2 (design §4.3 — the full
# spinning-secondary surface: inclined analytical, numeric-orbit(-correction), spherical family,
# equatorial Darwin). GROUNDING (Phase-7B fixture rule):
# (1) PRIMARY characterization vs the BAKED BHPToolkit MATHEMATICA-oracle fixtures
#     (ref_incl_spin_mode / ref_incl_spin_mode_k10 / ref_p5_comparson_dh / ref_incl_correction /
#      ref_from_correction / ref_spin_mode_num / ref_spin_mode_correction / ref_spin_mode_spherical /
#      ref_spin_mode_spherical_correction{_num,,_analytic}) at the FROZEN legacy tolerances copied
#     VERBATIM from packages/TeukolskySpinFluxes/test/runtests.jl (per-testset comments below).
# (2) Julia-to-Julia method identities (FromCorrection s=0 == Analytical spar=0; Mode s=0 ==
#     Analytical spar=0; Correction monopole == FromCorrection s=0; SphericalCorrectionNum monopole
#     == Spherical s=0; SphericalCorrection == …Num FluxesCorrection; Darwin == Phases) — analytic
#     cross-checks independent of any oracle.
# (3) Exported & callable surface check (the 12 TeukolskySpinMode* symbols).
# No live test-time dependency on packages/* or ref-code/* (fixtures baked, provenance in-file).
using Perturbation
using Perturbation: KerrSpinOrbitCorrection, KerrSpinOrbitCorrectionSpherical,
                    _orbit_derivatives_spherical

include("ref_incl_spin_mode.jl")
include("ref_incl_spin_mode_k10.jl")
include("ref_p5_comparson_dh.jl")
include("ref_incl_correction.jl")
include("ref_from_correction.jl")
include("ref_spin_mode_num.jl")
include("ref_spin_mode_correction.jl")
include("ref_spin_mode_spherical.jl")
include("ref_spin_mode_spherical_correction_num.jl")
include("ref_spin_mode_spherical_correction.jl")
include("ref_spin_mode_spherical_correction_analytic.jl")

@testset "SpinFluxes modes :: 12 TeukolskySpinMode* symbols exported" begin
    # the full design-§4.3 SpinFluxes surface (2 pre-existing + 10 ported here)
    for sym in (:TeukolskySpinMode, :TeukolskySpinModeFromCorrection, :TeukolskySpinModeAnalytical,
                :TeukolskySpinModeCorrection, :TeukolskySpinModeCorrectionAnalytical,
                :TeukolskySpinModeSpherical, :TeukolskySpinModeSphericalCorrectionNum,
                :TeukolskySpinModeSphericalCorrection, :TeukolskySpinModeSphericalCorrectionAnalytical,
                :TeukolskySpinModeCircularCorrectionAnalytical,
                :TeukolskySpinModeEquatorialCorrectionAnalyticalPhases,
                :TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin)
        @test isdefined(Perturbation, sym)
        @test getfield(Perturbation, sym) isa Function
    end
end

@testset "SpinFluxes modes :: INCLINED TeukolskySpinModeAnalytical vs MMA oracle (5 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "INCLINED full mode ..." testset): ω 1e-9, CP 2e-4,
    # ℐ-flux 2e-4 (atol 1e-10/1e-11 = .wl Chop), ℋ-flux 2e-2 (R.Up Ctrans artifact + tiny values).
    for r in REF_INCL_SPIN_MODE
        md = TeukolskySpinModeAnalytical(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x, r.spar)
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test isapprox(md.Amplitudes.I, r.CP; rtol=2e-4)
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=2e-4, atol=1e-10)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=2e-4, atol=1e-11)
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=2e-2, atol=1e-10)
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=2e-2, atol=1e-10)
    end
end

@testset "SpinFluxes modes :: DH P5 polar harmonic k=-10 TeukolskySpinModeAnalytical vs oracle" begin
    # FROZEN legacy tolerances (runtests.jl "DH P5 mode k=−10" testset): ω 1e-9, high-|k| 5e-3.
    for r in REF_INCL_SPIN_MODE_K10
        md = TeukolskySpinModeAnalytical(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x, r.spar)
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test isapprox(md.Amplitudes.I, r.CP; rtol=5e-3)
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=5e-3)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=5e-3)
    end
end

@testset "SpinFluxes modes :: P5 published golden (comparson_DH.nb Out[50], gauge pex)" begin
    # FROZEN legacy tolerances (runtests.jl "P5 headline reproduction" testset): steps exact,
    # ω 1e-12, amplitudes/fluxes 1e-8 (near machine — same analytic radial as the cached .nb).
    a, p, e, x = 0.9, 12.0, 0.2, 0.5   # the DH P5 orbit (comparson_DH.nb In[25])
    for r in REF_P5_COMPARSON_DH
        md = TeukolskySpinModeAnalytical(r.l, r.m, r.n, r.k, a, p, e, x, r.spar; gauge=r.gauge)
        @test md.stepsr == r.stepsr && md.stepsz == r.stepsz
        @test isapprox(md.omega, r.omega; rtol=1e-12)
        @test isapprox(md.Amplitudes.I, r.AmpI; rtol=1e-8)
        @test isapprox(md.Fluxes.Energy.I, r.EnergyI; rtol=1e-8)
        @test isapprox(md.Fluxes.Energy.H, r.EnergyH; rtol=1e-8)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.AngMomI; rtol=1e-8)
    end
end

@testset "SpinFluxes modes :: INCLINED TeukolskySpinModeCorrectionAnalytical vs oracle" begin
    # FROZEN legacy tolerances (runtests.jl "INCLINED linear-in-spin flux correction" testset):
    # ω 1e-9; ℐ 2e-3 (k=0) / 3e-2 (k=−10, high-ω radial); ℋ 5e-2 atol 1e-24.
    for r in REF_INCL_CORRECTION
        md = TeukolskySpinModeCorrectionAnalytical(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x)
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        rt = r.k == 0 ? 2e-3 : 3e-2
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=rt)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=rt)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=5e-2, atol=1e-24)
    end
end

@testset "SpinFluxes modes :: TeukolskySpinModeFromCorrection vs MMA oracle (3 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "TeukolskySpinModeFromCorrection" testset).
    for r in REF_FROM_CORRECTION
        md = TeukolskySpinModeFromCorrection(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x, r.nmax, r.kmax, r.spar)
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsr == r.stepsr && md.stepsθ == r.stepsth   # FD step-count heuristic agrees
        @test isapprox(md.Amplitudes.I, r.CP; rtol=2e-4)        # CPlus
        @test isapprox(md.Amplitudes.H, r.CM; rtol=2e-4)        # CMinus
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=2e-4, atol=1e-12)            # ℐ Ė
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=2e-4, atol=1e-12)   # ℐ L̇z
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=2e-3, atol=1e-13)            # ℋ Ė
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=2e-3, atol=1e-13)   # ℋ L̇z
        @test isapprox(md.Fluxes.CarterConstantK.I, r.cfKI; rtol=2e-4, atol=1e-12)  # ℐ K̇
        @test isapprox(md.Fluxes.CarterConstantK.H, r.cfKH; rtol=2e-3, atol=1e-13)  # ℋ K̇
    end
end

@testset "SpinFluxes modes :: FromCorrection s=0 == TeukolskySpinModeAnalytical spar=0 (Julia-to-Julia)" begin
    # FROZEN legacy identity testset: geodesic-monopole path isolation (rtol 1e-9, ω 1e-12).
    a, p, e, x = 0.9, 12.0, 0.2, 0.5
    for (l, m, n, k) in ((2,2,0,0), (2,1,0,0))
        fc = TeukolskySpinModeFromCorrection(l, m, n, k, a, p, e, x, 4, 4, 0.0)
        an = TeukolskySpinModeAnalytical(l, m, n, k, a, p, e, x, 0.0)
        @test isapprox(fc.omega, an.omega; rtol=1e-12)
        @test isapprox(fc.Amplitudes.I, an.Amplitudes.I; rtol=1e-9)
        @test isapprox(fc.Fluxes.Energy.I, an.Fluxes.Energy.I; rtol=1e-9)
        @test isapprox(fc.Fluxes.AngularMomentum.I, an.Fluxes.AngularMomentum.I; rtol=1e-9)
    end
end

@testset "SpinFluxes modes :: TeukolskySpinMode (full finite-s) vs MMA oracle (3 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "TeukolskySpinMode: full finite-s" testset).
    for r in REF_SPIN_MODE_NUM
        md = TeukolskySpinMode(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x, r.nmax, r.kmax, r.spar)
        @test isapprox(md.omega, r.omega; rtol=1e-9)           # full spin-corrected BL ω
        @test md.stepsr == r.stepsr && md.stepsθ == r.stepsth  # FD step-count heuristic agrees
        @test isapprox(md.Amplitudes.I, r.CP; rtol=2e-4)       # CPlus
        @test isapprox(md.Amplitudes.H, r.CM; rtol=2e-4)       # CMinus
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=2e-4, atol=1e-12)            # ℐ Ė
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=2e-4, atol=1e-12)   # ℐ L̇z
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=2e-3, atol=1e-13)            # ℋ Ė
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=2e-3, atol=1e-13)   # ℋ L̇z
    end
end

@testset "SpinFluxes modes :: TeukolskySpinMode s=0 == TeukolskySpinModeAnalytical spar=0 (Julia-to-Julia)" begin
    a, p, e, x = 0.9, 12.0, 0.2, 0.5
    for (l, m, n, k) in ((2,2,0,0), (2,1,0,0))
        sm = TeukolskySpinMode(l, m, n, k, a, p, e, x, 4, 4, 0.0)
        an = TeukolskySpinModeAnalytical(l, m, n, k, a, p, e, x, 0.0)
        @test isapprox(sm.omega, an.omega; rtol=1e-12)
        @test isapprox(sm.Amplitudes.I, an.Amplitudes.I; rtol=1e-9)
        @test isapprox(sm.Fluxes.Energy.I, an.Fluxes.Energy.I; rtol=1e-9)
        @test isapprox(sm.Fluxes.AngularMomentum.I, an.Fluxes.AngularMomentum.I; rtol=1e-9)
    end
end

@testset "SpinFluxes modes :: TeukolskySpinModeCorrection (monopole + dE¹) vs MMA oracle (3 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "TeukolskySpinModeCorrection" testset): geodesic Fluxes
    # 1e-6 ℐ / 1e-5 ℋ (same MST backend both sides), dE¹/dLz¹ 8e-3 (FD-limited h=3e-4 vs hh=1e-13).
    for r in REF_SPIN_MODE_CORRECTION
        md = TeukolskySpinModeCorrection(r.l, r.m, r.n, r.k, r.a, r.p, r.e, r.x, r.nmax, r.kmax)
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsr == r.stepsr && md.stepsθ == r.stepsth
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=1e-6, atol=1e-14)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=1e-6, atol=1e-13)
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=1e-5, atol=1e-15)
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=1e-5, atol=1e-14)
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=8e-3, atol=1e-14)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=8e-3, atol=1e-15)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=8e-3, atol=1e-13)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=8e-3, atol=1e-14)
    end
end

@testset "SpinFluxes modes :: Correction monopole == FromCorrection s=0 (Julia-to-Julia)" begin
    a, p, e, x = 0.9, 12.0, 0.2, 0.5
    for (l, m, n, k) in ((2,2,0,0), (2,1,0,0))
        mc = TeukolskySpinModeCorrection(l, m, n, k, a, p, e, x, 4, 4)
        fc = TeukolskySpinModeFromCorrection(l, m, n, k, a, p, e, x, 4, 4, 0.0)
        @test isapprox(mc.omega, fc.omega; rtol=1e-12)
        @test isapprox(mc.Amplitudes.I, fc.Amplitudes.I; rtol=1e-9)        # CPlus0 == FromCorrection s=0
        @test isapprox(mc.Fluxes.Energy.I, fc.Fluxes.Energy.I; rtol=1e-9)
    end
end

@testset "SpinFluxes modes :: TeukolskySpinModeSpherical vs MMA oracle (3 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "TeukolskySpinModeSpherical" testset).
    for r in REF_SPIN_MODE_SPHERICAL
        md = TeukolskySpinModeSpherical(r.l, r.m, r.k, r.a, r.p, r.e, r.x, r.kmax, r.spar)
        @test md.n == 0
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsθ == r.stepsth                            # FD step-count heuristic agrees
        @test isapprox(md.Amplitudes.I, r.CP; rtol=2e-4)        # CPlus
        @test isapprox(md.Amplitudes.H, r.CM; rtol=2e-4)        # CMinus
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=2e-4, atol=1e-12)            # ℐ Ė
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=2e-4, atol=1e-12)   # ℐ L̇z
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=2e-3, atol=1e-13)            # ℋ Ė
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=2e-3, atol=1e-13)   # ℋ L̇z
    end
end

@testset "SpinFluxes modes :: TeukolskySpinModeSpherical s=0 finite, dipole active, pre-built oc agrees" begin
    # FROZEN legacy testsets: s=0 flux finite & positive (ℐ) and s=0.1 dipole changes CPlus; both
    # method signatures agree at rtol 1e-12; e≠0 guarded at the orbit builder.
    oc = KerrSpinOrbitCorrectionSpherical(0.9, 12.0, 0.0, 0.5, 4)
    for (l, m, k) in ((2,2,0), (2,2,2), (2,1,0))
        g  = TeukolskySpinModeSpherical(l, m, k, oc, 0.0)
        sp = TeukolskySpinModeSpherical(l, m, k, oc, 0.1)
        @test isfinite(g.Amplitudes.I) && isfinite(g.Amplitudes.H)
        @test isfinite(g.Fluxes.Energy.I) && g.Fluxes.Energy.I > 0           # ℐ energy flux positive
        @test isfinite(g.Fluxes.Energy.H) && isfinite(g.Fluxes.AngularMomentum.I)
        @test !isapprox(g.Amplitudes.I, sp.Amplitudes.I; rtol=1e-6)          # s=0.1 dipole changes CPlus
    end
    m1 = TeukolskySpinModeSpherical(2, 2, 0, oc, 0.1)
    m2 = TeukolskySpinModeSpherical(2, 2, 0, 0.9, 12.0, 0.0, 0.5, 4, 0.1)
    @test isapprox(m1.Amplitudes.I, m2.Amplitudes.I; rtol=1e-12)
    @test isapprox(m1.Fluxes.Energy.I, m2.Fluxes.Energy.I; rtol=1e-12)
    @test_throws ArgumentError TeukolskySpinModeSpherical(2, 2, 0, 0.9, 12.0, 0.2, 0.5, 4, 0.1)  # e≠0 guarded
end

@testset "SpinFluxes modes :: TeukolskySpinModeSphericalCorrectionNum vs MMA oracle (3 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "...SphericalCorrectionNum" testset): δω MUST equal the
    # oracle's 1e-6 (FD matched both sides); geodesic ℐ 2e-4 / ℋ 2e-3; dE¹/dLz¹ 1e-3.
    δω = 1 / 10^6
    for r in REF_SPIN_MODE_SPHERICAL_CORR_NUM
        md = TeukolskySpinModeSphericalCorrectionNum(r.l, r.m, r.k, r.a, r.p, r.e, r.x, r.kmax, δω)
        @test md.n == 0
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsθ == r.stepsth
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=2e-4, atol=1e-14)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=2e-4, atol=1e-13)
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=2e-3, atol=1e-15)
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=2e-3, atol=1e-14)
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=1e-3, atol=1e-14)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=1e-3, atol=1e-16)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=1e-3, atol=1e-13)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=1e-3, atol=1e-15)
    end
end

@testset "SpinFluxes modes :: SphericalCorrectionNum monopole == Spherical s=0 + pre-built oc (Julia-to-Julia)" begin
    oc = KerrSpinOrbitCorrectionSpherical(0.9, 12.0, 0.0, 0.5, 4)
    for (l, m, k) in ((2,2,0), (2,1,0))
        mc = TeukolskySpinModeSphericalCorrectionNum(l, m, k, oc, 1/10^6)
        sp = TeukolskySpinModeSpherical(l, m, k, oc, 0.0)
        @test isapprox(mc.omega, sp.omega; rtol=1e-12)
        @test isapprox(mc.Amplitudes.I, sp.Amplitudes.I; rtol=1e-9)        # CPlus0 == Spherical s=0
        @test isapprox(mc.Fluxes.Energy.I, sp.Fluxes.Energy.I; rtol=1e-9)
    end
    m1 = TeukolskySpinModeSphericalCorrectionNum(2, 2, 0, oc, 1/10^6)
    m2 = TeukolskySpinModeSphericalCorrectionNum(2, 2, 0, 0.9, 12.0, 0.0, 0.5, 4, 1/10^6)
    @test isapprox(m1.FluxesCorrection.Energy.I, m2.FluxesCorrection.Energy.I; rtol=1e-12)
    @test isapprox(m1.Amplitudes.I, m2.Amplitudes.I; rtol=1e-12)
    @test_throws ArgumentError TeukolskySpinModeSphericalCorrectionNum(2, 2, 0, 0.9, 12.0, 0.2, 0.5, 4, 1/10^6)
end

@testset "SpinFluxes modes :: TeukolskySpinModeSphericalCorrection (+ r/x derivatives) vs MMA oracle (2 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "...SphericalCorrection" testset): geodesic Fluxes 1e-9
    # (same FD steps both sides ⇒ near-machine), FluxesCorrection 1e-7, FluxesDerivatives 1e-5.
    for r in REF_SPIN_MODE_SPHERICAL_CORR
        md = TeukolskySpinModeSphericalCorrection(r.l, r.m, r.k, r.a, r.p, r.e, r.x, r.kmax)
        @test md.n == 0
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsθ == r.stepsth
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=1e-9, atol=1e-18)
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=1e-9, atol=1e-20)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=1e-9, atol=1e-16)
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=1e-9, atol=1e-18)
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=1e-7, atol=1e-16)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=1e-7, atol=1e-18)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=1e-7, atol=1e-15)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=1e-7, atol=1e-17)
        @test isapprox(md.FluxesDerivatives.r.Energy.I, r.dEdrI; rtol=1e-5, atol=1e-16)
        @test isapprox(md.FluxesDerivatives.r.Energy.H, r.dEdrH; rtol=1e-5, atol=1e-18)
        @test isapprox(md.FluxesDerivatives.r.AngularMomentum.I, r.dLdrI; rtol=1e-5, atol=1e-15)
        @test isapprox(md.FluxesDerivatives.r.AngularMomentum.H, r.dLdrH; rtol=1e-5, atol=1e-17)
        @test isapprox(md.FluxesDerivatives.x.Energy.I, r.dEdxI; rtol=1e-5, atol=1e-16)
        @test isapprox(md.FluxesDerivatives.x.Energy.H, r.dEdxH; rtol=1e-5, atol=1e-18)
        @test isapprox(md.FluxesDerivatives.x.AngularMomentum.I, r.dLdxI; rtol=1e-5, atol=1e-15)
        @test isapprox(md.FluxesDerivatives.x.AngularMomentum.H, r.dLdxH; rtol=1e-5, atol=1e-17)
    end
end

@testset "SpinFluxes modes :: SphericalCorrection == …CorrectionNum + pre-built oc/odv (Julia-to-Julia)" begin
    oc = KerrSpinOrbitCorrectionSpherical(0.9, 12.0, 0.0, 0.5, 4)
    for (l, m, k) in ((2,2,0), (2,1,0))
        ma = TeukolskySpinModeSphericalCorrection(l, m, k, 0.9, 12.0, 0.0, 0.5, 4)
        mn = TeukolskySpinModeSphericalCorrectionNum(l, m, k, oc, 1/10^6)
        @test isapprox(ma.omega, mn.omega; rtol=1e-12)
        @test isapprox(ma.Amplitudes.I, mn.Amplitudes.I; rtol=1e-10)              # CPlus0 == Num geodesic
        @test isapprox(ma.Fluxes.Energy.I, mn.Fluxes.Energy.I; rtol=1e-10)
        @test isapprox(ma.FluxesCorrection.Energy.I, mn.FluxesCorrection.Energy.I; rtol=1e-10)   # dE¹ == Num
        @test isapprox(ma.FluxesCorrection.AngularMomentum.I, mn.FluxesCorrection.AngularMomentum.I; rtol=1e-10)
    end
    odv = _orbit_derivatives_spherical(0.9, 12.0, 0.0, 0.5, 4; h=1e-5)
    m1 = TeukolskySpinModeSphericalCorrection(2, 2, 0, oc, odv)
    m2 = TeukolskySpinModeSphericalCorrection(2, 2, 0, 0.9, 12.0, 0.0, 0.5, 4)
    @test isapprox(m1.FluxesCorrection.Energy.I, m2.FluxesCorrection.Energy.I; rtol=1e-12)
    @test isapprox(m1.FluxesDerivatives.r.Energy.I, m2.FluxesDerivatives.r.Energy.I; rtol=1e-12)
    @test isapprox(m1.FluxesDerivatives.x.Energy.I, m2.FluxesDerivatives.x.Energy.I; rtol=1e-12)
    @test_throws Exception TeukolskySpinModeSphericalCorrection(2, 2, 0, 0.9, 12.0, 0.2, 0.5, 4)  # e≠0 guarded
end

@testset "SpinFluxes modes :: TeukolskySpinModeSphericalCorrectionAnalytical (+ CarterK) vs oracle (2 modes)" begin
    # FROZEN legacy tolerances (runtests.jl "...SphericalCorrectionAnalytical" testset): geodesic
    # Fluxes 1e-9 (FD-independent, near-machine), FluxesCorrection 2e-3 (h vs hh FD mismatch).
    for r in REF_SPIN_MODE_SPHERICAL_CORR_ANALYTIC
        md = TeukolskySpinModeSphericalCorrectionAnalytical(r.l, r.m, r.k, r.a, r.p, r.x)
        @test md.n == 0
        @test isapprox(md.omega, r.omega; rtol=1e-9)
        @test md.stepsθ == r.stepsth
        @test isapprox(md.Fluxes.Energy.I, r.efI; rtol=1e-9, atol=1e-18)
        @test isapprox(md.Fluxes.Energy.H, r.efH; rtol=1e-9, atol=1e-20)
        @test isapprox(md.Fluxes.AngularMomentum.I, r.jfI; rtol=1e-9, atol=1e-16)
        @test isapprox(md.Fluxes.AngularMomentum.H, r.jfH; rtol=1e-9, atol=1e-18)
        @test isapprox(real(md.Fluxes.CarterConstantK.I), r.kfI; rtol=1e-9, atol=1e-16)
        @test isapprox(real(md.Fluxes.CarterConstantK.H), r.kfH; rtol=1e-9, atol=1e-18)
        @test imag(md.Fluxes.CarterConstantK.I) == 0 && imag(md.Fluxes.CarterConstantK.H) == 0
        @test isapprox(md.FluxesCorrection.Energy.I, r.dE1I; rtol=2e-3, atol=1e-16)
        @test isapprox(md.FluxesCorrection.Energy.H, r.dE1H; rtol=2e-3, atol=1e-18)
        @test isapprox(md.FluxesCorrection.AngularMomentum.I, r.dL1I; rtol=2e-3, atol=1e-15)
        @test isapprox(md.FluxesCorrection.AngularMomentum.H, r.dL1H; rtol=2e-3, atol=1e-17)
        @test isapprox(real(md.FluxesCorrection.CarterConstantK.I), r.dK1I; rtol=2e-3, atol=1e-15)
        @test isapprox(real(md.FluxesCorrection.CarterConstantK.H), r.dK1H; rtol=2e-3, atol=1e-17)
    end
end

@testset "SpinFluxes modes :: DARWIN == Phases (Julia-to-Julia, transitively MMA-verified)" begin
    # FROZEN legacy testset (runtests.jl "DARWIN equatorial flux correction == Phases"): the Darwin
    # χr-grid sum must reproduce the (oracle-verified, ref_eq_spin_flux.jl) Mino-phase equatorial
    # correction — two parametrizations of the SAME flux with the SAME central-FD h=3e-4 (FD error
    # cancels ⇒ near-machine). Legacy tolerances: ω 1e-12, geodesic Fluxes 1e-9, dE¹/dL¹ 1e-7.
    a, p, e = 0.7, 10.0, 0.2
    for (l, m, n) in ((2,2,0), (2,2,1), (2,1,0))
        mdD = TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin(l, m, n, a, p, e)
        mdP = TeukolskySpinModeEquatorialCorrectionAnalyticalPhases(l, m, n, a, p, e)
        @test mdD.k == 0 && mdD.n == n                                    # Darwin NamedTuple shape
        @test isapprox(mdD.omega, mdP.omega; rtol=1e-12)
        @test isapprox(mdD.Fluxes.Energy.I, mdP.Fluxes.Energy.I; rtol=1e-9)            # geodesic Ė ℐ
        @test isapprox(mdD.Fluxes.AngularMomentum.I, mdP.Fluxes.AngularMomentum.I; rtol=1e-9)
        @test isapprox(mdD.FluxesCorrection.Energy.I, mdP.FluxesCorrection.Energy.I; rtol=1e-7)  # dE¹ ℐ
        @test isapprox(mdD.FluxesCorrection.AngularMomentum.I, mdP.FluxesCorrection.AngularMomentum.I; rtol=1e-7)
    end
end
