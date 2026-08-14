# WaSABI fixture data (Phase-7B fixture invariant)

Byte-identical copies of the exact upstream WaSABI files the Hybrid test suite reads,
baked here so the committed tests carry NO live dependency on `ref-code/WaSABI`
(a fresh clone runs green without it).

## Provenance

- Source: `/data/haiyangw/claude/emri-disk/ref-code/WaSABI`
  (clone of https://github.com/BlackHolePerturbationToolkit/WaSABI.git)
- Git state at capture: commit `5a76f6f5d1033e8a315fdc161dab688d55c0b2e8`
  ("also adding - sign convention for the 22-mode", 2026-01-29), working tree clean.
- Captured: 2026-07-02, verified byte-identical (`cmp`) to the source at copy time.

## Contents (only what the tests read)

Model catalogue plugins — `InspiralModels/*.m`, `AmplitudeModels/*.m` (7 each):
used by `test_registry.jl` to characterize the upstream `FileNames["*.m", dir]`
plugin-discovery contract (basenames == the sorted 7-model catalogue). Contents are
carried verbatim for provenance; only the `sf_data` inspiral root (see below) is
actually parsed at test time via `get_inspiral("WaSABI-C", …)`.

Interpolation data read by `HybridC.wasabic_model(inspiral_datadir, amp_datadir)`
(see `src/HybridC/wasabic_deps.jl` and `src/Models/inspiral_0PAKerr.jl`), used by
`test_hybridc_deps.jl`, `test_hybridc_codegen.jl`, `test_hybridc_inspiral.jl`,
`test_registry.jl`, `test_introspection.jl`:

- `InspiralModels/sf_data/1SF_Flux/Kerr_Circ/SMRfluxdata2025_36x36.data`
- `InspiralModels/sf_data/1SF_Flux/Schwarz_Circ/1SFCircShwarzDotEInf.m`
- `InspiralModels/sf_data/2SF_Flux/Schwarz_Circ/2SFCircShwarzDotEInf.m`
- `InspiralModels/sf_data/Leading_S2_Flux/Kerr_Circ/SMRfluxdataWithSecondarySpin2025_36x36.m`
- `InspiralModels/sf_data/1SF_Local_Invariants/Kerr_Circ/KerrCircularEquatorialInvariants2025b_36x36.data`
- `InspiralModels/sf_data/stop_conditions/StopHybridFixedChi.m`
- `AmplitudeModels/sf_amp_data/1SF/Kerr_Circ/ReAmp1SFKerr22NewGrid.m`
- `AmplitudeModels/sf_amp_data/1SF/Kerr_Circ/ImAmp1SFKerr22NewGrid.m`

The `InspiralModels/sf_data` ↔ sibling `AmplitudeModels/sf_amp_data` layout is
load-bearing: `Registry._wc_amp_dir` derives the amplitude root from the inspiral
root (`joinpath(dirname(dirname(sf)), "AmplitudeModels", "sf_amp_data")`).

## Tolerance rationale

All oracle VALUES in the tests are pre-existing baked WL literals (oracle16a/16b,
`wolframscript` over the upstream paclet) — this fixture only replaces the DATA
INPUT path; every tolerance is unchanged from the original characterization.
