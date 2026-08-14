# Force.KerrFrontend — the a≠0 metric front-end THROUGH MetricReconstruction.Lorenz (the design
# §4.10 consumption edge; §7.3: AccretionDisc reaches the Lorenz machinery ONLY via
# MetricReconstruction, never SelfForce/Kernels). GROUNDING — the legacy ref_kerr_* gates at the
# legacy FROZEN tolerances, on a REDUCED mode set / grid (runtime budget; explicitly sanctioned:
# "small but physics-meaningful grids/mode sets"):
#   • G-K-F8 (a=0.5): the F8 torque balance ∂_r J̇^r_Adv == ∂T^φ_Mat closes mean_rel < 1e-6 — an
#     ANALYTIC OPERATOR IDENTITY of the forced master ODE, valid for ANY smooth real source, hence
#     unaffected by the reduced lmax/grid (legacy gdkf8 gate, frozen 1e-6).
#   • G-K-F8 a→0 regression: the SAME gate through the Kerr path at a=0 (legacy gdkf8_a0, 1e-6) —
#     ties the Kerr front-end to the corpus-validated a=0 closure.
#   • G-K-RES: |𝔉| peaks within 1.5M of the spin-shifted Kerr OLR (legacy gres, frozen 1.5), and
#     prograde spin pushes the OLR outward (legacy gres_spin_trend).
#   • G-K-PAR (§8): kerr_mode threaded == serial BITWISE (race-free per-slot + blas_single_thread).
#   • parity/packing: even (l+m) ⇒ 7 BLS fields, odd ⇒ 3 (the corpus convention).
# The ILR-side driver (f8_torque_balance_ilr) carries the legacy honest-scope disclosure: the
# off-equatorial source assembly inside the corotation band is slow/fragile at a≠0 — a MANUAL
# driver, not asserted here (verbatim the legacy suite's disclosed non-gate).
using Test
using AccretionDisc
using AccretionDisc.Force
using Printf

# 22-pt grid bracketing the OLR + the forced-solve span [10.4,16] with ≥2 guard points each end.
const KERR_RG = collect(8.6:0.4:17.0)

@testset "Force.KerrFrontend (a≠0 via MetricReconstruction.Lorenz)" begin
    m = 2; r_p = 8.0

    @testset "G-K-F8 + G-K-RES + G-K-PAR at a=0.5 (one shared closure build)" begin
        a = 0.5
        t0 = time()
        modes, build = kerr_modes(a, r_p, 2:3, m; rgrid = KERR_RG, lmax_proj = 8,
                                  return_build = true)
        @printf("  kerr_modes(a=0.5, l=2:3, lmax_proj=8, %d pts): %.1f s\n", length(KERR_RG), time() - t0)
        # parity/packing convention (l+m even ⇒ 7 even BLS fields; odd ⇒ 3)
        @test modes[2].parity === :even && modes[2].nfields == 7
        @test modes[3].parity === :odd  && modes[3].nfields == 3
        @test all(all(isfinite, modes[l].val_right[j]) for l in 2:3 for j in 1:modes[l].nfields)

        # §8 G-K-PAR: threaded == serial BITWISE on a subgrid, reusing the shared closure build
        sub = collect(8.6:1.6:16.6)
        mser = kerr_mode(a, r_p, 2, m; rgrid = sub, recs = build.recs, S = build.S, ls = build.ls,
                         threaded = false)
        mpar = kerr_mode(a, r_p, 2, m; rgrid = sub, recs = build.recs, S = build.S, ls = build.ls,
                         threaded = true)
        @test all(mser.val_right[j] == mpar.val_right[j] for j in 1:mser.nfields)   # BITWISE (§8)
        @test all(mser.dr_right[j] == mpar.dr_right[j] for j in 1:mser.nfields)

        # G-K-F8: the operator identity closes through the Kerr 𝔉 (legacy frozen 1e-6)
        d = disc(a)
        res = f8_torque_balance(d, m, r_p, modes, 3; rg = 10.2:0.8:16.0, N = 16000,
                                r_start = 16.0, r_end = 10.4, rmin = 10.6, rmax = 15.4,
                                verbose = false)
        @printf("  G-K-F8 a=0.5: mean_rel=%.3e max_rel=%.3e\n", res.mean_rel, res.max_rel)
        @test res.mean_rel < 1e-6

        # G-K-RES: |𝔉| peaks at the spin-shifted Kerr OLR (legacy frozen 1.5M)
        rI, rO = lindblad_locations(d, r_p, m)
        f5 = f5_source(d, m, r_p, modes, 3; rg = 8.8:0.1:13.0, verbose = false)
        @printf("  G-K-RES a=0.5: peak r=%.2f, Kerr OLR=%.3f, |Δ|=%.2f\n",
                Float64(f5.r_peak), Float64(rO), abs(Float64(f5.r_peak) - Float64(rO)))
        @test abs(Float64(f5.r_peak) - Float64(rO)) < 1.5
    end

    @testset "G-K-F8 a→0 regression (Kerr path reproduces the a=0 closure)" begin
        modes0 = kerr_modes(0.0, r_p, 2:3, m; rgrid = KERR_RG, lmax_proj = 6)
        d0 = disc(0.0)
        res = f8_torque_balance(d0, m, r_p, modes0, 3; rg = 10.2:0.8:16.0, N = 16000,
                                r_start = 16.0, r_end = 10.4, rmin = 10.6, rmax = 15.4,
                                verbose = false)
        @printf("  G-K-F8 a=0 (Kerr path): mean_rel=%.3e max_rel=%.3e\n", res.mean_rel, res.max_rel)
        @test res.mean_rel < 1e-6
    end

    @testset "Kerr OLR spin trend (prograde spin pushes the OLR outward)" begin
        o = [Float64(lindblad_locations(disc(a), r_p, m)[2]) for a in (0.0, 0.5, 0.9)]
        @printf("  OLR(a=0,0.5,0.9) = %.4f, %.4f, %.4f\n", o[1], o[2], o[3])
        @test o[2] > o[1] && o[3] > o[2]                              # legacy gres_spin_trend gate
    end
end
