# arXiv:2506.19917 (Fairbairn & Dittmann 2025) — shared HOLE stub of the scaffold rows.
#
# Every MANIFEST.toml row of this paper points here until its own thin driver lands
# (ECC::fd25-repro-scaffold, emri-disk EMRI-ECC ledger).  It contains NO physics and writes NO
# golden: `--check` and `--regen` both exit 3 (documented HOLE), which run_all counts as SKIP,
# never as PASS.  When a driver lands, the row's `script` moves to its own `figNN.jl` and its
# `status` follows the ledger node.
#
#   julia --project=reproduce reproduce/arXiv-2506.19917/_future.jl --check

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.hole("arXiv:2506.19917 scaffold — no driver and no golden yet for this row " *
                  "(MANIFEST status future, mirrors the ECC::fd25-* ledger nodes); " *
                  "see reproduce/arXiv-2506.19917/README.md")
