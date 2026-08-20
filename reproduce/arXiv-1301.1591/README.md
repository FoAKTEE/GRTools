# arXiv:1301.1591 — Harms, Bernuzzi & Brügmann, *Numerical solution of the 2+1 Teukolsky equation on a hyperboloidal and horizon penetrating foliation of Kerr*

Reproduced with the generic-spin (s = ±1, ±2) VACUUM m-mode ℓ-coupled sector of
`Perturbation.TimeDomain` (the paper's own RT horizon-penetrating hyperboloidal chart:
`SpinModesProblem` / `solve_spin_modes`, `Spacetime.RTChart` maps).

Migrated from the legacy campaign `emri-disk results/dg-teukolsky-spin/` (scripts
`run_spin_tails.jl`, `run_convergence.jl`, `run_fig14.jl`, `run_scatter.jl`, `fit_rates.jl`,
`driver_common.jl`), which already ran on this solver. The migration moved the campaign's
figure-local ESTIMATORS into the package — they produce published numbers, so they are physics:

| legacy figure-local function | GrTools home |
|---|---|
| `driver_common.envelope_rate` | `Perturbation.TimeDomain.fit_tail_envelope` |
| `driver_common.adaptive_rate` | `Perturbation.TimeDomain.fit_tail_adaptive` |
| `fit_rates.sliding_fits` | `Perturbation.TimeDomain.fit_tail_sliding` |
| `fit_rates.fit_rate_K` | `Perturbation.TimeDomain.fit_convergence_order` |
| `driver_common.excise_policy` (+ its private `rt_R_of_r` copy) | `Perturbation.TimeDomain.spin_excise_policy` (chart maps reused from `Spacetime.RTChart`) |

gated by `capabilities/Perturbation/test/test_timedomain_tailfits.jl` (45/45), whose gate L
recomputes the legacy campaign's published self-convergence rates from its tracked ladders.

## Results

| legacy artefact | driver | result id | tier | status | vs legacy |
|---|---|---|---|---|---|
| `F-S1.4_Tails_Kerr_s-2_a0.9_finite` + `F-S1.4s_…_scri`, `data/tails_s-2_a0p9_K24.csv`, `verifier_tails.csv` (a=0.9 rows) | `tails_a0p9.jl` | `tails_a0p9` | fast | solid | **SAME** — `dt` and `maxu_end` bitwise, all nine fit rows to every printed digit |
| `F-S1.4n_Tails_spin_ladder`, `data/tails_s-2_a0p{99,995,9995}_K24.csv` | `tails_spin_ladder.jl` | `tails_spin_ladder` | heavy | solid | **SAME** — all nine ℓ=2 fits to every printed digit |
| `F-S1.1_Convergence_s-2_RT`, `data/conv_s-2_a0p9_N{2,4,6}_tf50.csv`, `verifier_conv.csv` | `fig_convergence.jl` | `convergence` | heavy | solid | **SAME** — every ladder entry to 9 digits, every rate to 2 decimals |
| `F-S1.14_LPI_radii_s-2_a0.9_lp4`, `data/fig14_s-2_a0p9_lp4.csv` | `fig14_lpi_radii.jl` | `fig14_lpi_radii` | heavy | solid | **SAME** — `dt`/`maxu_end` bitwise; the documented Float64 NO-GO reproduced (envelope LPIs stall at ≈ −5.5 against the expected −8/−7) |
| `T-S1.1_ScatteringMatrix_s-2_a0.8` (m=2 ID1 rows), `data/scatter_s-2_a0p8_m2_lp{2,3}_ID1.csv` | `scatter_m2_ID1.jl` | `scatter_m2_ID1` | heavy | solid | **SAME** — fits equal `verifier_scatter.csv` to every printed digit |
| `T-S1.1_ScatteringMatrix_s-2_a0.8` (m=0 ID0 rows), `data/scatter_s-2_a0p8_m0_lp{2,3,4,5}_ID0.csv` | `scatter_m0_ID0.jl` | `scatter_m0_ID0` | heavy | see MANIFEST | see MANIFEST |
| `figures/diag_kappa_sweep.png`, `data/diag*.csv` | — | — | — | NOT MIGRATED | diagnosis record, not a published object (see below) |

`common_tails.jl` / `common_scatter.jl` hold the shared RUN CONFIGURATION (the committed P4 gate
config) and the printed law tables; they contain no physics.

## What the paper's numbers are, and what is measured

The reproduced published objects are the empirical late-time tail laws Eqs. (28),(29) and the
printed decay-rate tables (Tab. `s2_LPI` for ID0/m = 0, Tab. `m_neq_0_rates_sm2` for ID1/m = 2),
plus Fig. 14's LPI-vs-extraction-radius structure.

Two boundaries are MEASURED and carried into the goldens rather than hidden:

* **the κ-mediated |s| = 2 oscillation** — the committed constraint damping `constraint_kappa = 0.3`
  drives a mesh-frequency oscillation that rides every s = −2 tail (absent at s = 0), so a raw
  log-log fit lands within ≈ 0.7 of the law instead of on it. That is the committed gate tolerance
  (±1.0), and it is why `fit_tail_envelope` exists. No knob is retuned per case.
* **the Float64 junk floor** — entries of the scattering table with printed rates μ ≳ 8 sit below the
  numerical floor and cannot be measured on this stack; the paper itself marks its μ ≥ 10 entries ×
  even at quadruple precision, and its Fig. 14 required quadruple precision. `fig14_lpi_radii` is
  therefore a reproduction of a documented NO-GO: the golden stores the measured (non-matching) LPIs.

## Figure checks (VLM, per the recurring reproduction procedure)

* `tails_a0p9.png` vs legacy `F-S1.4_Tails_Kerr_s-2_a0.9_finite.png` — **STRUCTURE MATCH**: same
  quantity and frame (log10|ψ_ℓ| vs log10 τ), same ψ₂/ψ₃/ψ₄ colour-and-dash scheme, same
  "(expected law, measured envelope LPI)" legend, same oscillation envelope, same ordering and the
  same floor onset near log10 τ ≈ 3.2. Two DIFFERENCES, both stated rather than smoothed over:
  (i) the legacy render carries an inset with the RUNNING envelope LPI, not reproduced here — the
  golden stores the envelope LPI per probe and ℓ instead; (ii) the legend's envelope LPIs differ
  (ψ₂ −6.558 here vs −6.325 there) because the legacy RENDERER fitted its own window while this
  driver reports the verifier's committed P4 window (700, 1300) — and −6.558 is exactly the value
  in the legacy verifier table `data/verifier_tails.csv`, which is the authoritative record.
* `convergence.png` vs legacy `F-S1.1_Convergence_s-2_RT.png` — **MATCH**: same 1×3 panel structure
  (R₊ / R_mid / scri), same three per-N ladders with the same colours and markers, same per-N slope
  legends carrying identical rates and classes, same ordering (scri steeper than the interior at
  every N), same data ranges.

## Not migrated

`figures/diag_kappa_sweep.png` and `data/diag*.csv` (legacy `run_diag.jl` / `run_diag2.jl`, 16 solver
knob variants) are the DIAGNOSIS record that identified the κ-oscillation, not a reproduction of a
published object; the finding they established is documented above and in the drivers. Re-running
the sweep costs ≈ 16 × 200 s at -t4 and is reproduced by
`run_diag2.jl`-equivalent variations of `tails_a0p9.jl`'s configuration (`constraint_kappa`,
`filter_alpha`, `N`, `K`, `excise`), which the driver exposes through `common_tails.run_tail`.
