# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_h1.jl (untracked one-shot, mtime 2026-06-30)
#               sha256(original) = 58f7dbfebdf045c9a4a4f54f67ed41dd74f920a0842c72f88661dd0c42fa0d05
# fixture(s)  : h1lorenz/h1-l{0..3}m*.h5  (the H1Lorenz 4-part writer layout)
#               NOTE: the original one-shot wrote into .../test/fixtures/h1; that directory was
#               later RENAMED to h1lorenz when the 6-column SecondOrderRicci corpus took the
#               name `h1` (see test/lorenz_solver.jl:4-8).  The output path here follows the
#               current tree.  Verified 2026-08-19: the 6 modes that finished inside the re-run
#               budget (l0m0, l1m0, l1m1, l2m0, l2m1, l2m2) come back BIT-IDENTICAL to
#               fixtures/h1lorenz/ (all datasets + attributes, max abs difference 0.0).
# oracle      : legacy packages/H1Lorenz @ 394b7844 (compute_one_mode/write_mode_h5)
# run         : julia test/fixtures/capture/capture_h1.jl   # self-activates a temp env
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const LEGACY_PACKAGES_ROOT = get(ENV, "GRTOOLS_LEGACY_PACKAGES_ROOT", "/data/haiyangw/claude/emri-disk/packages")
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# One-shot offline capture: regenerate GENUINE h1 fixtures from the legacy H1Lorenz oracle on the
# small 7-point grid the SelfForce.Solver test uses (grid=[4,5,6,8,10,12,14], r0=8 → r0index=3,
# ra=0, rb=6). Writes the legacy (4nf,nL)/(4nf,nR) layout + grid attrs via the legacy write_mode_h5.
# NOT part of the test (H1Lorenz is a HEAVY oracle, never a live dep) — run by hand to bake fixtures.
import Pkg
Pkg.activate(mktempdir())
Pkg.develop(path="$(LEGACY_PACKAGES_ROOT)/H1Lorenz")
Pkg.add("HDF5")
using H1Lorenz

grid = [4.0,5.0,6.0,8.0,10.0,12.0,14.0]
ra, r0i, rb = 0, 3, 6   # r0 = grid[r0i+1] = 8
p = 8.0
outdir = "$(FIXTURE_DIR)/h1lorenz"

modes = [(0,0),(1,0),(1,1),(2,0),(2,1),(2,2),(3,0),(3,1),(3,2),(3,3)]
for (l,m) in modes
    nmode, orbit = compute_one_mode(p, l, m; grid=grid, indices=(ra,r0i,rb))
    path = joinpath(outdir, "h1-l$(l)m$(m).h5")
    write_mode_h5(path, nmode, orbit)
    println("wrote $path  nf=$(nmode.nf)  kind=$(nmode.cset_kind)  inhom[1,1,end]=", nmode.inhom[1,1,end])
end
println("DONE")
