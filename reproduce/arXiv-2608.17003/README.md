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
| 3  | fig:angular_momentum_balance         | `fig03.jl` | `expected/fig03.json` | heavy | **solid** | **match** | a•=15M, χ=0, h0=0.05, (𝔭,𝔮)=(0,0), m=1:120. Balance (71) away from corotation 4.0e-5; the corotation spike is gated on the invariant ∫(dF_J/dr−d𝔗/dr)dr = 0.0266073 F_{J,0} = T_L−T of (75) to 2.7e-5 (its plotted height is a plot-sampling number). 605 s at -t4. |
| 4  | fig:dtdr_ap — torque density vs a•    | `fig04.jl` | — | heavy | not started | — | |
| 5  | fig:TGR_vs_TNewt                     | `fig05.jl` | `expected/fig05.json` | heavy | **solid** | **match** | T_GR/T_Newt vs χ (10 points, −0.9…0.9) for the five (𝔭,𝔮) at h₀=0.05, a•=12M, b=0.6, m=1:80, grid (30,30,3000). T_Newt is the FULLY Newtonian pipeline (`SpiralWaves` `model = :newtonian`, EGD::sw-newtonian-limit) and is χ-independent by construction. Caption claim M:930-932 (𝔭<0 enhances, 𝔭≥0 suppresses) PASSES on all 50 markers. Against the published markers (read at 200 dpi, ±0.015): 46/50 within 0.05, mean |Δ| ≤ 0.034 per series. FINDING: the (1.5,0) series drifts to −0.107 at χ=+0.9 (−12 %), not a T_Newt normalisation (the factor is not constant along the series) — see the evidence dir. 3470 s at -t8. |
| 6  | fig:torque_chi_ap                    | `fig06.jl` | — | heavy | not started | — | |
| 7  | fig:torque_ap_h0                     | `fig07.jl` | — | heavy | not started | — | |
| 8  | fig:m_95_plot — m95 violins          | `fig08.jl` | — | heavy | not started | — | |
| 9  | fig:num_vs_analytical                | `fig09.jl` | `expected/fig09.json` | heavy | **preliminary** | **partial** (reduced a• grid) | Analytic branch = `Modes.hkr_lindblad_torque` (arXiv:2510.03564 (18b)+(20) summed by (24)); numerical branch = `SpiralWaves.torque_survey` (Eq. 75); cutoff m₄₅ from our own mode sum. ALL FIVE printed numbers PASS: m₄₅h₀ = 0.35 (paper 0.35), (0,0) ratio 1.4844 @ 12M (1.5) and 0.9235 @ 100M (0.92), (0,1) band [0.8737, 0.9363] ([0.88, 0.96]). REDUCED to 7 a• points (paper ~19) for compute; m NOT reduced (top-10 of 100 modes ≤ 1.01e-3 of the torque). Bit-identical across two independent runs. 5350 s at -t8. |
| 10 | fig:LR_torque_and_reversal           | `fig10.jl` | `expected/fig10.json` | heavy | **preliminary** | **match** (structure) | T_{L,m} for (𝔭,h₀) ∈ {(5,0.03),(5,0.05),(0,0.03),(0,0.05)}, χ=0, 𝔮=0, a•=15M, b=0.6; plotted m=4:40, totals summed to ⌈4/h₀⌉ at the reference grid (240,240,24000). M:1014 PASSES (last m with T_{L,m}<0 = 10 at h₀=0.05 and 16 at h₀=0.03; total torque (74) = +1.539/+1.841 F_{J,0}h₀ > 0). FINDING: our T_{L,m} is ≈13 % below the published figure at low m (ratio 0.87 at m=4–10 → 1.02 at m=40); domain-extent and boundary-branch hypotheses refuted (evidence dir). ~28 min at -t8. |
| 11 | fig:convergence                      | `fig11.jl` | `expected/fig11.json` | heavy | **solid** | **match** | a•=25M, χ=0, h0=0.05, m=1:200, triples (30,3000)/(60,6000)/(120,12000) vs (240,240,24000). Median rel errors 1.8e-4 / 5.0e-6 / 9.9e-7. UNIT FINDING: the unlabelled absolute panel of the paper is in F_{J,0}, not F_{J,0}h0 — measured from a tick-calibrated pixel read (factor 1/h0 = 20 between the two panels' agreement). 798 s at -t4 for the reference triple. |
| 12 | fig:pole_disp_errors                 | `fig12.jl` | — | heavy | not started | — | |
| 13 | fig:variation_with_respect_to_b      | `fig13.jl` | — | heavy | not started | — | |
| 14 | fig:m95_vs_b                         | `fig14.jl` | — | heavy | not started | — | |
| —  | printed-number gate                  | `check_printed_numbers.jl` | none (by design) | fast | **HOLE** (14 PASS / 0 FAIL / 4 uncovered) | n/a | asserts every printed number of the paper that the figure goldens cover; still waiting on fig08 (m₉₅ band) and fig13/fig14 (softening) |

`_profile_assembly.jl` is a paper-local, PRESENTATION-side include shared by `fig03.jl` and
`fig04.jl` (so the assembly is written once): it runs the m-loop over `SpiralWaves.solve_mode` and
reduces each mode's `amf` / `cumulative_torque` to cell averages on a common plotting grid — the
mode sum of (73).  It contains no physics; its natural long-term home is
`AccretionDisc.SpiralWaves.diagnostics` (flagged [FUTURE] in its header).

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
