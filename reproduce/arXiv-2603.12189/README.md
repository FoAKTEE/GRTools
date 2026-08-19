# arXiv:2603.12189 — the shifted-geodesic approximation for a spinning secondary

Reproduced with GrTools `Perturbation.SpinFluxes` on top of `Spacetime.SpinningSecondary`.
Legacy sources: `emri-disk/results/gw-flux/figures/2603.12189/repro_*.jl` (read-only).

| published object | driver | result id | tier | status | wall (-t4) | vs legacy |
|---|---|---|---|---|---|---|
| fig. 3 `fig:fluxcomparen` | `fig03.jl` | `fig03` | fast | solid | 124 s | VLM: SAME |
| fig. 7 `fig:comparepvalues_leadingharmonics` | `fig07.jl` | `fig07` | fast | solid | 206 s | — |

`[N/A]`: fig. 1a–e are schematic cartoons with no data.
`[FUTURE]`: figs. 11 and 12 (inspiral dephasing) need the un-ported 0PA driver.
Not yet ported (legacy status recorded in
`emri-disk/results/gw-flux/reference/FIGURE_REPRODUCTION_MATRIX.md`): fig. 2 left panel and
figs. 8/9 are `(p, e)` grids the legacy analysis measured as infeasible (days of compute); fig. 2
right, fig. 4, fig. 5 and fig. 6 are heavy drivers still to be written.

## Promotion made for these drivers (one home, no duplicates)

* `Spacetime.SpinningSecondary.shifted_geodesic_correction` — commit `dd4965d`. The
  shifted-geodesic trajectory (the orbit correction with its libration zeroed) was defined inline
  in `repro_fig4_E1.jl` and `repro_fig7.jl`; it now has one package home with its own testset, and
  swapping `fig07.jl` onto it reproduces the golden bitwise.

## What the drivers reproduce

* **fig. 3** — the EXACT `(2,2,n,0)` linear-in-spin flux corrections versus the radial harmonic
  n = 1..11 in four channels (E and Lz, at infinity and at the horizon), against the
  SHIFTED-GEODESIC approximation. The exact curve is the full numeric orbit correction (with the
  oscillatory δχ/δr libration); SG is the central spin finite difference of the **pex-gauge** full
  flux (secular shift only). Measured: SG tracks the exact correction to 5.0e-5 (E∞), 7.7e-6 (EH),
  5.2e-5 (Lz∞), 7.7e-6 (LzH) at n = 1 and deviates to 5.1 % / 6.4 % / 5.1 % / 6.4 % at n = 11 — the
  published message that the shifted geodesic degrades with harmonic number.
  Note the gauge: the CTilde-locked correction entry point is a *different* gauge at this orbit and
  lands ~7.7× off with the opposite sign, so it is not the SG curve.
* **fig. 7** — the error of three approximate spinning trajectories (shifted geodesic, and the
  orbit correction truncated at n^s = k^s = 4 and 6) against the exact n^s = k^s = 8 trajectory,
  for the single (2,2,0,0) mode over p ∈ [6, 15]. The driver registers, per point, whether the
  published ordering SG > n^s=4 > n^s=6 holds: it does, at **19/19 points in all four channels**.
  At p = 6 the E-at-infinity errors are 1.92e-6 (SG), 3.24e-10 (n^s=4) and 4.79e-13 (n^s=6) — three
  to four decades per added harmonic pair.
