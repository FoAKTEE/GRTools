# arXiv:gr-qc/0703028 — Sundararajan, Khanna & Hughes, *Towards adiabatic waveforms for inspiral into Kerr black holes: I. A new model of the source for the time domain perturbation equation*

This directory holds the published objects of the **DG-Teukolsky s = −2 Kerr sourced-flux campaign**
migrated from emri-disk `results/dg-teukolsky-spin/` (`README-flux-figures.md`, `data-flux/`).

## Why this paper directory and not a "DG flux figure" one

The legacy campaign produced two artefacts, `F1.7-analogue_s-2_kerr_ell_convergence` and
`T1.3-analogue_s-2_kerr_flux_table`. Its own verifier records that they are **not reproductions of a
published object** — the DG-Teukolsky papers publish no s = −2 figure or table (a recorded P0
finding of that campaign); their reference is the in-package frequency-domain oracle plus the
published SKH table. The reproducible published object is therefore the SKH table, and it lives here.
Carrying the "analogues" as reproductions would be a false claim.

## Results

| published object | driver | result id | tier | status |
|---|---|---|---|---|
| Table I (column FD, \|m\| = 2), rows (a = 0.9, r₀ = 10M) and (a = 0.5, r₀ = 10M); Table II (FD, \|m\| = 3), row (a = 0.9, r₀ = 10M) | `check_skh_tables.jl` | `skh_tables` | fast | solid |

Measured relative differences against the published four-digit values: **1.363e-5** (a = 0.9, |m| = 2),
**6.435e-5** (a = 0.5, |m| = 2), **2.358e-5** (|m| = 3) — agreement to every published digit, and the
same numbers the legacy campaign recorded in `data-flux/dgs2_skh_anchor.out`.

## The time-domain side, and what is not migrated

The campaign's `data-flux/c01…c15.dat` (15 configurations, 44 rows) are runs of the Kerr sourced
time-domain sector `Perturbation.TimeDomain.spin_sourced_flux_kerr` against the frequency-domain
oracle. They are **not migrated as reproduce drivers**: a single (N, K) = (10, 64) configuration cost
1105 s of wall time per row-block and the campaign is ≈ 4 h; the physics is already gated inside the
package by `capabilities/Perturbation/test/test_timedomain_spinsource_kerr.jl` (gates K1–K6, 105/105,
7 m 53 s), which includes the SKH anchor gate this driver now also registers as a golden.
Re-running the campaign is `spin_sourced_flux_kerr` at the configurations listed in
`README-flux-figures.md`. [FUTURE — compute-bound, physics in the package and gated]
