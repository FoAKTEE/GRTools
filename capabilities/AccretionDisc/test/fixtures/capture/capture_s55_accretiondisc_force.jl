# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_s55_accretiondisc_force.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = 06f8d46b4ce73f9117a8a9d4bd05beb156508c9ec7923e356d86ef8526a2cf80
# fixture(s)  : h1-l{2..6}m2.h5
# oracle      : legacy packages/DiscEMRIForce @ fee0036b + packages/DiscModes @ 6ca6457a; corpus ref-code/h1Lorenz/run_r8
# run         : julia --project=$GRTOOLS_LEGACY_PACKAGES_ROOT/DiscEMRIForce -t4 test/fixtures/capture/capture_s55_accretiondisc_force.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const EMRIDISK_ROOT        = get(ENV, "GRTOOLS_EMRIDISK_ROOT",        "/data/haiyangw/claude/emri-disk")
const REFCODE_ROOT         = get(ENV, "GRTOOLS_REFCODE_ROOT",         "/data/haiyangw/claude/emri-disk/ref-code")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# =============================================================================================================
# ONE-SHOT CAPTURE (task S5-5, 2026-07-02) — legacy DiscEMRIForce + DiscModes oracle → baked fixtures/literals
# for the grtools AccretionDisc.Force pipeline port (design §4.10). Runs the LEGACY packages ONCE; the committed
# grtools tests then depend ONLY on the fixtures written here + the printed literals baked into the test files
# (Phase-7B fixture rule: no live packages/* or ref-code/* reach from committed tests).
#
# Run: julia --project=packages/DiscEMRIForce -t4 scratchpad/capture_s55_accretiondisc_force.jl
#
# PART A  — truncated a=0 h1 BLS corpus fixtures: ref-code/h1Lorenz/run_r8/fields/h1-l{2..6}m2.h5, right grid
#           windowed to r ∈ [r0=8, 17.6] (276 pts) + the last 3 left-grid points (r0index=3). Values are BYTE
#           COPIES of the corpus (no processing) — per-point pipeline quantities on the window are bit-identical
#           to the full-corpus run (verified in PART B).
# PART B  — truncation-invariance verification + oracle gate numbers on the fixtures:
#           G-D1 loader anchor, force_mode f_r / Vvec V_r at r≈11 (the legacy cov_source_fluid refs),
#           F8 torque balance (G-D2 ~1e-7), F8 flux cross-check (G-D3), F5 source peak-vs-OLR.
# PART C  — DiscModes off-equator oracle literals (OmegaK_z / alpha_z / Uvec_z / pressure_z / eden_z /
#           build_Qinv_z / build_Q_z / heightavg) for the grtools Modes.OffEquator characterization tests.
# =============================================================================================================
using DiscEMRIForce, Printf
import DiscModes
import HDF5

const OUT = "$(FIXTURE_DIR)"
mkpath(OUT)
const H1DIR = DiscEMRIForce.H1DIR
const RMAX = 17.6
const NLEFT = 3

println("=== PART A: truncated corpus fixtures (l=2..6, m=2, right window r0..$(RMAX), $(NLEFT) left pts) ===")
for l in 2:6
    HDF5.h5open(joinpath(H1DIR, "h1-l$(l)m2.h5"), "r") do f
        g = HDF5.read(f["grid"]); r0i = Int(HDF5.read_attribute(f["grid"], "r0index")[1])  # 0-based
        R = HDF5.read(f["inhom_right"]); L = HDF5.read(f["inhom_left"])
        Cinf = haskey(f, "C_inf") ? HDF5.read(f["C_inf"]) : Float64[]
        Chor = haskey(f, "C_horiz") ? HDF5.read(f["C_horiz"]) : zeros(Float64, length(Cinf))
        nkeep = findlast(x -> x <= RMAX, g[r0i+1:end])
        # new grid: NLEFT pts strictly left of r0, then r0, then the right window
        newgrid = vcat(g[r0i+1-NLEFT:r0i], g[r0i+1:r0i+nkeep])
        newL = L[:, r0i+1-NLEFT:r0i+1]        # old left cols ↔ grid[1:r0i+1]; keep the last NLEFT+1 (incl r0)
        newR = R[:, 1:nkeep]                  # old right col k ↔ grid[r0i+k]; col 1 = r0
        HDF5.h5open(joinpath(OUT, "h1-l$(l)m2.h5"), "w") do o
            o["grid"] = newgrid
            HDF5.attributes(o["grid"])["raindex"] = Int32(0)
            HDF5.attributes(o["grid"])["r0index"] = Int32(NLEFT)
            HDF5.attributes(o["grid"])["rbindex"] = Int32(length(newgrid) - 1)
            o["inhom_left"] = newL
            o["inhom_right"] = newR
            o["C_inf"] = Cinf
            o["C_horiz"] = Chor
            HDF5.attributes(o["grid"])["provenance"] = "one-shot S5-5 capture 2026-07-02: byte-copied window " *
                "[r0,$(RMAX)] of ref-code/h1Lorenz/run_r8/fields/h1-l$(l)m2.h5 (iter94 corpus) via " *
                "scratchpad/capture_s55_accretiondisc_force.jl"
        end
        @printf("  l=%d: grid %d→%d pts (r0=%.6f, last=%.6f), nf=%d\n",
                l, length(g), length(newgrid), newgrid[NLEFT+1], newgrid[end], size(R, 1) ÷ 4)
    end
end

println("\n=== PART B: truncation invariance + oracle gate numbers ===")
d = disc_a0(); m = 2; r_p = 8.0
ω0 = m * DiscModes.OmegaK(d, r_p)
modes_fix  = load_modes(2:6, m; dir = OUT)
modes_full = load_modes(2:6, m)

# raw-value byte identity on the window
for l in 2:6
    mf = modes_fix[l]; mo = modes_full[l]
    n = length(mf.r_right)
    @assert mf.r_right == mo.r_right[1:n] "grid mismatch l=$l"
    @assert all(mf.val_right[j] == mo.val_right[j][1:n] for j in 1:mf.nfields) "val mismatch l=$l"
    @assert all(mf.dr_right[j] == mo.dr_right[j][1:n] for j in 1:mf.nfields) "dr mismatch l=$l"
    @assert mf.r0 == mo.r0 && mf.parity == mo.parity
end
println("  fixture right-grid values BYTE-IDENTICAL to full corpus on the window: OK")

# G-D1 anchor
v = modes_fix[2].val_right[1][1]
@printf("  G-D1 anchor field1(r0) fixture = %.17g %+.17gim   |Δ vs 3.16711-0.23710im| = %.3e\n",
        real(v), imag(v), abs(v - (3.16711 - 0.23710im)))

# force_mode / Vvec at r≈11 (legacy cov_source_fluid refs)
gf = modes_fix[2].r_right; k = argmin(abs.(gf .- 11.0))
gg = modes_full[2].r_right; kfull = argmin(abs.(gg .- 11.0))
f_fix  = force_mode(d, modes_fix[2], k, 2, 2, ω0, 0.0)
f_full = force_mode(d, modes_full[2], kfull, 2, 2, ω0, 0.0)
@printf("  force_mode f_r (r=%.6f): fixture = %.17g %+.17gim (bitmatch full: %s)\n",
        gf[k], real(f_fix[2]), imag(f_fix[2]), f_fix == f_full)
V_fix  = Vvec(d, modes_fix, k, m, r_p, ω0, 0.0; lmax = 6)
V_full = Vvec(d, modes_full, kfull, m, r_p, ω0, 0.0; lmax = 6)
@printf("  Vvec V_r (r=%.6f): fixture = %.17g %+.17gim (bitmatch full: %s)\n",
        gf[k], real(V_fix[2]), imag(V_fix[2]), V_fix == V_full)
sF_fix  = source_F(d, modes_fix, k, m, r_p, ω0; lmax = 6)
sF_full = source_F(d, modes_full, kfull, m, r_p, ω0; lmax = 6)
@printf("  source_F (r=%.6f): fixture = %.17g %+.17gim (bitmatch full: %s)\n",
        gf[k], real(sF_fix), imag(sF_fix), sF_fix == sF_full)

# G-D2 F8 torque balance
res_fix  = f8_torque_balance(d, m, r_p, modes_fix, 6;  N = 16000, verbose = false)
res_full = f8_torque_balance(d, m, r_p, modes_full, 6; N = 16000, verbose = false)
@printf("  G-D2 F8: fixture mean_rel=%.3e max_rel=%.3e mean_int=%.3e | full mean_rel=%.3e\n",
        res_fix.mean_rel, res_fix.max_rel, res_fix.mean_rel_integral, res_full.mean_rel)

# G-D3 flux cross-check
fx_fix = f8_flux_crosscheck(d, m, r_p, modes_fix, 6; N = 12000, verbose = false)
@printf("  G-D3 flux: fixture cv_full=%.3e cv_wave=%.3e (wave<full: %s)\n",
        fx_fix.cv_full, fx_fix.cv_wave, fx_fix.cv_wave < fx_fix.cv_full)

# F5 source structure
f5_fix = f5_source(d, m, r_p, modes_fix, 6; rg = 8.2:0.25:15.0, verbose = false)
@printf("  F5: fixture peak r=%.3f OLR=%.6f |Δ|=%.3f\n",
        f5_fix.r_peak, f5_fix.OLR, abs(f5_fix.r_peak - f5_fix.OLR))

# source_F_grid probes (cov_kerr_grid analogues on the fixture window)
Fg = source_F_grid(d, modes_fix, m, r_p, ω0; lmax = 6, npts = 9, spectral = false)
rI, rO = DiscModes.lindblad_locations(d, r_p, m)
kprot = argmin(abs.(gf .- rO))
Fp_point = source_F(d, modes_fix, kprot, m, r_p, ω0; lmax = 6)
ksm = argmin(abs.(gf .- (r_p + 6.0)))
Fsm_point = source_F(d, modes_fix, ksm, m, r_p, ω0; lmax = 6)
@printf("  source_F_grid: byte@OLR=%s  rel_wing@r≈14=%.3e  finite=%s\n",
        Fg[kprot] == Fp_point, abs(Fg[ksm] - Fsm_point)/(abs(Fsm_point) + eps()), all(isfinite, Fg))

println("\n=== PART C: DiscModes off-equator oracle literals (grtools Modes.OffEquator characterization) ===")
dd = DiscModes.RelativisticDisc(a = 0.6, h0 = 0.1, r_I = 1.0, eden_I = 1.0)
r = 12.0; z = 0.05; rp = 15.0; mm = 2
@printf("  probe (a=0.6, h0=0.1, r=12, z=0.05, r_p=15, m=2):\n")
@printf("  OmegaK_z   = %.17g\n", DiscModes.OmegaK_z(dd, r, z))
@printf("  alpha_z    = %.17g\n", DiscModes.alpha_z(dd, r, z))
U = DiscModes.Uvec_z(dd, r, z)
@printf("  Uvec_z     = [%.17g, %.17g, %.17g, %.17g]\n", U[1], U[2], U[3], U[4])
@printf("  pressure_z = %.17g\n", DiscModes.pressure_z(dd, r, z))
@printf("  eden_z     = %.17g\n", DiscModes.eden_z(dd, r, z))
Qi = DiscModes.build_Qinv_z(dd, r, z, rp, mm)
Q  = DiscModes.build_Q_z(dd, r, z, rp, mm)
for (i, j) in ((1,1), (2,2), (1,4), (4,2), (3,3), (2,1))
    @printf("  Qinv[%d,%d] = %.17g %+.17gim    Q[%d,%d] = %.17g %+.17gim\n",
            i, j, real(Qi[i,j]), imag(Qi[i,j]), i, j, real(Q[i,j]), imag(Q[i,j]))
end
# equatorial (z=0) probe for the z→0 reduction cross-check vs grtools qinv_operator
Qi0 = DiscModes.build_Qinv_z(dd, r, 0.0, rp, mm)
@printf("  z=0: Qinv[1,1]=%.17g %+.17gim  Qinv[2,2]=%.17g %+.17gim  Qinv[4,2]=%.17g %+.17gim\n",
        real(Qi0[1,1]), imag(Qi0[1,1]), real(Qi0[2,2]), imag(Qi0[2,2]), real(Qi0[4,2]), imag(Qi0[4,2]))
# a=0 probe in the Force-test region
d0 = DiscModes.RelativisticDisc(a = 0.0, h0 = 0.1, r_I = 1.0, eden_I = 1.0)
Q0 = DiscModes.build_Q_z(d0, 11.0, 0.02, 8.0, 2)
@printf("  a=0 (r=11, z=0.02, r_p=8, m=2): Q[1,1]=%.17g %+.17gim  Q[2,2]=%.17g %+.17gim\n",
        real(Q0[1,1]), imag(Q0[1,1]), real(Q0[2,2]), imag(Q0[2,2]))
# heightavg oracle samples (Gaussian moments)
@printf("  heightavg(z->z^2, 0.1)          = %.17g  (exact h^2 = 0.01)\n", DiscModes.heightavg(zz -> zz^2, 0.1))
@printf("  heightavg(z->cos(z), 0.1)        = %.17g\n", DiscModes.heightavg(zz -> cos(zz), 0.1))
@printf("  heightavg(z->z^4, 0.1; order=4)  = %.17g  (exact 3h^4 = 3e-4)\n",
        DiscModes.heightavg(zz -> zz^4, 0.1; order = 4))
@printf("  heightavg(z->cos(z), 0.1; hz=0.02) = %.17g\n", DiscModes.heightavg(zz -> cos(zz), 0.1; hz = 0.02))

println("\ncapture done.")
