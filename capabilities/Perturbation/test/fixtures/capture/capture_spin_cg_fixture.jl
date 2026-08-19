# ==========================================================================================================
# TRACKED FIXTURE CAPTURE SCRIPT  (GRR::g2-fixture-capture-scripts, repro-audit gap G2)
#
# origin      : emri-disk/scratchpad/capture_spin_cg_fixture.jl (untracked one-shot, mtime 2026-07-02)
#               sha256(original) = 10176231a33e352910446bf7daad8adff55aa31bf903342a8bb03dac844b7e16
# fixture(s)  : spin_cg_reference.jl
# oracle      : core/Wigner.clebsch_gordan (BigInt-exact Racah), GRTools in-tree
# run         : julia --project=core/Wigner capabilities/Perturbation/test/fixtures/capture/capture_spin_cg_fixture.jl
#
# Every absolute path below was rewritten into an ENV-overridable constant (see the
# preamble): GRTOOLS_LEGACY_PACKAGES_ROOT, GRTOOLS_REFCODE_ROOT, GRTOOLS_DGTEUK_ROOT,
# GRTOOLS_EMRIDISK_ROOT, GRTOOLS_ROOT, GRTOOLS_FIXTURE_OUT_DIR.  Defaults reproduce the
# capture host.  Set GRTOOLS_FIXTURE_OUT_DIR to a scratch dir to re-run WITHOUT
# overwriting the committed fixture, then diff.
# The script body is otherwise VERBATIM from the original one-shot capture.
# ==========================================================================================================

# --- GRTools tracked-capture preamble (ENV-overridable paths) ---------------------------------
const FIXTURE_DIR          = get(ENV, "GRTOOLS_FIXTURE_OUT_DIR",      normpath(joinpath(@__DIR__, "..")))
# ----------------------------------------------------------------------------------------------

# ONE-SHOT capture (DG-P4 gate B): bake Clebsch-Gordan-form reference values of the
# spin-weighted cosθ coupling matrix elements into a Perturbation test fixture.
#
#   ⟨ₛY_{ℓ′m}|cosθ|ₛYₗₘ⟩ = √((2ℓ+1)/(2ℓ′+1)) · ⟨ℓ,m;1,0|ℓ′,m⟩ · ⟨ℓ,−s;1,0|ℓ′,−s⟩
#
# (dossier DG_SPIN_FORMULATION.md §2.2 DERIVED/CHECK3 — CG equivalent of arXiv:1702.04204
# Eq. (34)). Runs in the grtools core/Wigner environment (BigInt-exact Racah CG); Wigner is
# NOT a Perturbation dep — the baked fixture keeps the firewall.
#
# Usage:  julia --project=grtools/core/Wigner scratchpad/capture_spin_cg_fixture.jl
#
# Grid: s ∈ {0,±1,±2}, ℓ = |s|..12, |m| ≤ ℓ, ℓ′ ∈ {ℓ−1, ℓ, ℓ+1} (ℓ′ ≥ 0; unphysical ℓ′
# rows are exact zeros via the CG selection rules — they gate the closed form's guards too).

using Wigner

const OUT = joinpath(FIXTURE_DIR, "spin_cg_reference.jl")

cgcos(s, l, m, lp) = sqrt((2l + 1) / (2lp + 1)) *
                     clebsch_gordan(l, m, 1, 0, lp, m) *
                     clebsch_gordan(l, -s, 1, 0, lp, -s)

open(OUT, "w") do io
    println(io, "# AUTO-GENERATED one-shot fixture — DO NOT EDIT.")
    println(io, "# Baked from grtools core/Wigner clebsch_gordan (BigInt-exact Racah) by")
    println(io, "# scratchpad/capture_spin_cg_fixture.jl (julia $(VERSION)).")
    println(io, "# Formula: <sY_l'm|cosθ|sY_lm> = sqrt((2l+1)/(2l'+1))·CG(l,m;1,0|l',m)·CG(l,-s;1,0|l',-s)")
    println(io, "# — dossier DG_SPIN_FORMULATION.md §2.2 (DERIVED, CHECK3), the CG-form gate of")
    println(io, "# the closed-form coefficients [arXiv:1702.04204, Eq. (34)].")
    println(io, "# Grid: s ∈ {0,±1,±2}, l = |s|..12, |m| ≤ l, l' ∈ {l-1,l,l+1} (l' ≥ 0).")
    println(io, "# Tuple layout: (s, l, m, lp, value).")
    println(io, "const SPIN_CG_REF = [")
    for s in (-2, -1, 0, 1, 2), l in abs(s):12, m in -l:l, lp in (l - 1, l, l + 1)
        lp < 0 && continue
        v = cgcos(s, l, m, lp)
        println(io, "    ($s, $l, $m, $lp, $(repr(v))),")
    end
    println(io, "]")
end
println("wrote $(OUT)")
