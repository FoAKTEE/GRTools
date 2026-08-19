# arXiv:2608.17003 — Dittmann, Hegade K. R. & Fairbairn, "Extreme-Mass-Ratio Inspirals in Gaseous Disks"

Mirrors: tex `ref-paper/arXiv-2608.17003/src/` · supplementary notebook
`ref-code/arXiv-2608.17003-anc/` (emri-disk repo).  Physics lives in
`capabilities/AccretionDisc` (`SpiralWaves` = the paper's linear theory; `Modes` = the
pressureless analytic baselines it is compared against).  Drivers here are THIN.

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
| 1  | fig:spirals — δU maps                | `fig01.jl` | — | heavy | not started | — | |
| 2  | fig:spirals_conv — m_max convergence | `fig02.jl` | — | heavy | not started | — | |
| 3  | fig:angular_momentum_balance         | `fig03.jl` | — | fast  | not started | — | |
| 4  | fig:dtdr_ap — torque density vs a•    | `fig04.jl` | — | heavy | not started | — | |
| 5  | fig:TGR_vs_TNewt                     | `fig05.jl` | — | heavy | not started | — | |
| 6  | fig:torque_chi_ap                    | `fig06.jl` | — | heavy | not started | — | |
| 7  | fig:torque_ap_h0                     | `fig07.jl` | — | heavy | not started | — | |
| 8  | fig:m_95_plot — m95 violins          | `fig08.jl` | — | heavy | not started | — | |
| 9  | fig:num_vs_analytical                | `fig09.jl` | — | heavy | **blocked** | — | analytic half landed (`Modes.hkr_*`); the numerical T_L needs `SpiralWaves.solve_mode`/diagnostics (EGD::sw-mode-solver, EGD::sw-diagnostics) |
| 10 | fig:LR_torque_and_reversal           | `fig10.jl` | — | heavy | not started | — | |
| 11 | fig:convergence                      | `fig11.jl` | — | heavy | not started | — | |
| 12 | fig:pole_disp_errors                 | `fig12.jl` | — | heavy | not started | — | |
| 13 | fig:variation_with_respect_to_b      | `fig13.jl` | — | heavy | not started | — | |
| 14 | fig:m95_vs_b                         | `fig14.jl` | — | heavy | not started | — | |
| —  | printed-number gate                  | `check_printed_numbers.jl` | none (by design) | fast | **HOLE** (18/18 uncovered) | n/a | asserts every printed number of the paper that the figure goldens cover |

Style: `_style/emri_style.jl` + the VLM style manual
`progress/emri-gasdisk/FIGURE_STYLE_MANUAL.md` (emri-disk repo); every figure goes through the
3-step `vlm-figure-repro` procedure before its status may leave "preliminary".

## `check_printed_numbers.jl` — the printed-number gate (EGD::printed-numbers-gate)

Loads every `expected/<id>.json` in this directory and asserts the paper's printed numbers.
Exit codes: **0** all claims covered and passing · **1** a covered claim violated ·
**3** (HOLE) at least one claim's golden key does not exist yet.  It stores no number of its own,
so it is registered in `MANIFEST.toml` WITHOUT an `expected` key.

### Golden-key contract

A figure driver satisfies a claim by registering the named keys with `@repro_result`.  Missing
keys are reported by name — a claim is never silently dropped.

| claim | src | provenance | golden keys | paper value (tol) |
|---|---|---|---|---|
| `M919-peak-offset`        | M:919  | text | `fig04.peak_offset_over_h0a_newt_outer`, `fig04.peak_offset_over_h0a_newt_inner` | ±1.0 h0a• (±0.15) |
| `M919-peak-shift`         | M:919  | text | `fig04.peak_inward_shift_over_h0a_a10` | 0.2 h0a• (±0.1) |
| `M919-peak-growth`        | M:919  | text | `fig04.peak_growth_factor_a10` | 1.5 (±0.25) |
| `M925-no-reversal`        | M:925  | text | `fig07.total_torque_all_positive` | true |
| `M947-h0-turnover`        | M:947  | text | `fig07.torque_peak_h0`, `fig07.torque_monotone_up_to_peak` | h0 = 0.125, monotone |
| `M1025-m95h0-band`        | M:1025 | text | `fig08.m95_h0_min`, `fig08.m95_h0_max` | [2.0, 3.5] (±0.25) |
| `M992-m95-inverse-h0`     | M:992  | text | `fig08.m95_h0_spread_fraction` | ≤ 0.35 |
| `M999-m45-cutoff`         | M:999  | text | `fig09.m45_times_h0` | 0.35 (±0.1) |
| `M1002-ratio-pq00-a12`    | M:1002 | figure-read | `fig09.ratio_pq00_a12` | 1.5 (±0.15) |
| `M1002-ratio-pq00-a100`   | M:1002 | figure-read | `fig09.ratio_pq00_a100` | 0.92 (±0.06) |
| `M1002-ratio-pq01-band`   | M:1002 | figure-read | `fig09.ratio_pq01_min`, `fig09.ratio_pq01_max` | [0.88, 0.96] (±0.04) |
| `M1014-p5-negative-h005`  | M:1014 | text | `fig10.m_last_negative_p5_h005` | 10 (±3) |
| `M1014-p5-negative-h003`  | M:1014 | text | `fig10.m_last_negative_p5_h003` | 15 (±4) |
| `M1014-p5-total-positive` | M:1014 | text | `fig10.total_torque_p5_positive` | true |
| `A162-pole-scaling`       | A:162  | text | `fig12.relerr_ratio_1e3_over_1e4`, `fig12.relerr_ratio_1e4_over_1e5` | ≈10 each (factor 3) |
| `A164-torque-vs-b`        | A:164  | figure-read | `fig13.torque_h005_b03`, `fig13.torque_h005_b08` | 1.28 → 0.97 (±0.12) |
| `A164-m95-vs-b`           | A:164  | figure-read | `fig14.m95h0_b03`, `fig14.m95h0_b08` | 2.9 → 1.95 (±0.25) |
| `M893-ringing`            | M:893  | text | `fig02.ringing_residual_m10`, `fig02.ringing_residual_m30` | residual(30) < residual(10)/3 |

`provenance = text` means the value is printed verbatim at the quoted tex line;
`figure-read` means it was read off the reference figure in the wave-0 brief §5 — those carry
looser tolerances and are labelled `fig` in the gate table, never promoted to printed numbers.
