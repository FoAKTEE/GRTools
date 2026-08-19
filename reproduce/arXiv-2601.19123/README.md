# arXiv:2601.19123 — Dyson & D'Orazio, relativistic disc modes and the EMRI torque balance

| legacy script (`results/emri-disk/figures/2601.19123/…`) | driver | result id | tier | status | vs legacy |
|---|---|---|---|---|---|
| `repro_DampedRegions.jl` | `figDampedRegions.jl` | `figDampedRegions` | fast | solid | SAME (identical resonance condition) |
| `verify_isco_markers.jl` | `check_isco_markers.jl` | `check_isco_markers` | fast | solid | SAME |
| `repro_OperatorPlot.jl`, `repro_F2_operatorplot.jl` | `figOperatorPlot.jl` | `figOperatorPlot` | fast | preliminary | SAME leading-order potential; B^{(0)} still [HOLE] |

## Reproduced

* **F3 fig:DampedRegions** — the m = 1…15 Lindblad intervals at r_p = 15M for a = 0.01, 0.60, 0.99
  and the Newtonian limit. The bands narrow onto r_p with m and widen with spin
  (m=1: [8.689, 22.671] at a=0.01 → [6.560, 23.117] at a=0.99).
* **F7/F9 ISCO markers** — `Spacetime.kerr_isco`, prograde and retrograde, max relative deviation
  from the Bardeen textbook values 8.8e-5 (prograde, limited by the 5-digit textbook entries) and
  1.2e-5 (retrograde).
* **F2 fig:OperatorPlot** — the leading geometric potential B^{(−2)}_m/h² (promoted as
  `AccretionDisc.B_m_geometric`): its zeros sit on the Lindblad turning points to 4e-15 and its
  corotation pole coefficient matches `B_m_corotation` to 1.5e-7.

## Not reproduced here — corpus- or PMLA-gated (documented, never fabricated)

| published result | blocker |
|---|---|
| F1 fig:metricdata | needs the KerrLorenzCirc per-ℓ metric records; the ℓ-ramp is blocked by FD-tower roundoff + a ₂F₁ MST hang (a package fix, mission node `migrate-2601.19123` open obligation) |
| F5 fig:sourceplot, F8 fig:torquebalance, RT-A7/A8/A9 | need the `h1Lorenz run_r8` field corpus (139 MB, 136 HDF5 files, untracked) and/or the KerrLorenzCirc `.jls` metric; the drivers must locate them through the `Data` package root and call `ReproHarness.hole(...)` when absent |
| F4 fig:homogenoussolutons, F6 fig:densities2d, F7 fig:densitiesspindep | need the unpublished PMLA height-averaged operator B^{(0)} ([HOLE: H-DISC-FULLOP]); only the morphology/turning-point structure is currently grounded |
| F9 fig:AdvAngSpinVariation | regenerable from `AccretionDisc.lindblad_torque` but ~hours at the published spin/m grid — heavy tier, not yet driven |

The `derivations/` directory is reserved for the E2 CAS certificates (`E2cas/*.wls`, `E2_GATE.md`)
once the wolframscript gate is available in this environment.
