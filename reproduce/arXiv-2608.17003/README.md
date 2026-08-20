# arXiv:2608.17003 — Dittmann, Hegade K. R. & Fairbairn, "Extreme-Mass-Ratio Inspirals in Gaseous Disks"

Mirrors: tex `ref-paper/arXiv-2608.17003/src/` · supplementary notebook
`ref-code/arXiv-2608.17003-anc/` (emri-disk repo).  Physics lives in
`capabilities/AccretionDisc` (`SpiralWaves` = the paper's linear theory; `Modes` = the
pressureless analytic baselines it is compared against).  Drivers here are THIN.

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
| 1  | fig:spirals — δU maps                | `fig01.jl` | `expected/fig01.json` | heavy | **solid** | **match** | 2×3 δU M/m• maps, a•=15M, m=1:120, App.-C triple (120,120,12000), polar mesh r_ISCO…40M on 601² px. Winding compared with the published raster BY NUMBER (arm tracked in `spiral_grid.png` itself): turns(18M→35M) = −0.718/−1.084/−1.483 (χ=−0.9), −0.707/−1.067/−1.461 (χ=+0.9) vs published −0.721/−1.097/−1.504, −0.707/−1.081/−1.425 — worst \|Δ\| 0.036 turns. All three caption claims are driver gates and pass. Amplitude agrees to 4–6 % in the one unclipped column. 3936 s at -t8. |
| 2  | fig:spirals_conv — m_max convergence | `fig02.jl` | `expected/fig02.json` | heavy | **solid** | **match** | 𝔭=2, 𝔮=1, h0=0.1, a•=15M; m_max ∈ {10,30,120} as partial sums of ONE 120-mode solve per row. M:893 gate: ringing residual(10)=0.5095, residual(30)=0.1245 < residual(10)/3=0.1698, residual(40)=0.0615 < 0.10; ladder 0.746/0.510/0.250/0.125/0.061/0.015/0.002 monotone. Winding independent of the ceiling: turns(18M→35M) = −1.475/−1.482/−1.483 vs published −1.445/−1.478/−1.456. m_max=120 column reproduces fig01 panel 3 independently. 1801 s at -t8. |
| 3  | fig:angular_momentum_balance         | `fig03.jl` | `expected/fig03.json` | heavy | **solid** | **match** | a•=15M, χ=0, h0=0.05, (𝔭,𝔮)=(0,0), m=1:120. Balance (71) away from corotation 4.0e-5; the corotation spike is gated on the invariant ∫(dF_J/dr−d𝔗/dr)dr = 0.0266073 F_{J,0} = T_L−T of (75) to 2.7e-5 (its plotted height is a plot-sampling number). 605 s at -t4. |
| 4  | fig:dtdr_ap — torque density vs a•    | `fig04.jl` | `expected/fig04.json` | heavy | **solid** | **match** | χ=𝔭=𝔮=0, h0=0.05, a•∈{10,30,100,300,1000}M, m=1:120. All four M:919 numbers pass; pixel-calibrated reference agreement 0.00–0.25 % on the measurable extrema. 3779 s at -t4. |
| 5  | fig:TGR_vs_TNewt                     | `fig05.jl` | `expected/fig05.json` | heavy | **solid** | **match** | T_GR/T_Newt vs χ (10 points, −0.9…0.9) for the five (𝔭,𝔮) at h₀=0.05, a•=12M, b=0.6, m=1:80, grid (30,30,3000). T_Newt is the FULLY Newtonian pipeline (`SpiralWaves` `model = :newtonian`, EGD::sw-newtonian-limit) and is χ-independent by construction. Caption claim M:930-932 (𝔭<0 enhances, 𝔭≥0 suppresses) PASSES on all 50 markers. Against the published markers (read at 200 dpi, ±0.015): 46/50 within 0.05, mean |Δ| ≤ 0.034 per series. FINDING: the (1.5,0) series drifts to −0.107 at χ=+0.9 (−12 %), not a T_Newt normalisation (the factor is not constant along the series) — see the evidence dir. 3470 s at -t8. |
| 6  | fig:torque_chi_ap                    | `fig06.jl` | `expected/fig06.json` | heavy | **solid** | **match** | T (F_{J,0}h₀) vs a•/M ∈ [10,1000] (40 log-spaced points, counted from the reference marker pitch 13.565 px on a 264.5 px/decade axis) for χ ∈ {0.9, 0, −0.9}, panels (1.5,1) and (0,0), h₀=0.05, b=0.6, m=1:80, grid (30,30,3000). M:925 (no torque reversal) PASSES for all 240 configurations; M:936 spin ordering at a•=10M PASSES; the three spins converge at a•=1000M to 1.6 %. Plateau T(χ=0, a•=1000M) = 1.9664 / 1.1242 F_{J,0}h₀ against ≈1.97 / ≈1.13 read off the reference (0.2 % / 0.5 %). FINDING: the 1–2 smallest radii sit 5–15 % below the published markers (worst at a•=10M, χ ≤ 0) — same regime as the fig05/fig10 findings. 8616 s at -t8. NOTE: the y-limits were re-measured off the reference PDF after the run, so the golden's `script_sha256` predates that plotting-only edit; the numbers are untouched (redrawn with `--plot-only`). |
| 7  | fig:torque_ap_h0                     | `fig07.jl` | `expected/fig07.json` | heavy | **solid** | **match** | T (F_{J,0}h₀) vs a•/M ∈ [10,1000] (same 40 points as fig06) for h₀ ∈ {0.05,0.075,0.1,0.125,0.15} at χ=0, panels (1.5,1) and (0,0), b=0.6, m ≤ ⌈4/h₀⌉, grid (30,30,3000); 400 configurations, 6224 s at -t8. **Both printed claims PASS**: M:925 (T > 0 everywhere) and M:947 (the plateau torque of the (0,0) panel is 1.1196 → 1.1262 → 1.1722 → **1.2027** → 1.1335 over h₀ = 0.05…0.15, i.e. peak at 0.125, monotone below) — these are the `total_torque_all_positive` / `torque_peak_h0` / `torque_monotone_up_to_peak` keys the printed-number gate reads. Plateau ordering and values agree with the reference to ≈1 %. Cross-check: the h₀=0.05 series of the (0,0) panel is BITWISE the χ=0 series of fig06 (same 40 configurations, two drivers). FINDING (same as fig06): the smallest radii sit 5–15 % low. |
| 8  | fig:m_95_plot — m95 violins          | `fig08.jl` | `expected/fig08.json` | heavy | **preliminary** | **match** (structure) | m₉₅h₀ violins over h₀ ∈ {0.05,0.075,0.1,0.125,0.15} × χ ∈ {−0.9,0,0.9} for 5 (𝔭,𝔮) × a• ∈ {15,30,120}M, b=0.6; 225 surveys / 10 485 mode solves at (60,60,6000) with an in-driver reference-triple cross-check (same m₉₅, 3/3); 110 min at -t8. M:992 PASSES (worst spread 0.302 ≤ 0.35). M:1025 FAILS as the gate implements it: our band is [1.80, 3.00] and the gate wants max ∈ [3.25, 3.75] — but a pixel read-back of the PUBLISHED figure gives [2.208, 3.157], which fails the same criterion, so the gate band is the defect (see the evidence dir). FINDING: at the ⌈6/h₀⌉ ceiling our band is [1.95, 3.15], reproducing the published maximum 3.157 to 0.007 — with fig14 and fig13 this is three independent signs that the paper summed ≈6/h₀ modes. |
| 9  | fig:num_vs_analytical                | `fig09.jl` | `expected/fig09.json` | heavy | **preliminary** | **partial** (reduced a• grid) | Analytic branch = `Modes.hkr_lindblad_torque` (arXiv:2510.03564 (18b)+(20) summed by (24)); numerical branch = `SpiralWaves.torque_survey` (Eq. 75); cutoff m₄₅ from our own mode sum. ALL FIVE printed numbers PASS: m₄₅h₀ = 0.35 (paper 0.35), (0,0) ratio 1.4844 @ 12M (1.5) and 0.9235 @ 100M (0.92), (0,1) band [0.8737, 0.9363] ([0.88, 0.96]). REDUCED to 7 a• points (paper ~19) for compute; m NOT reduced (top-10 of 100 modes ≤ 1.01e-3 of the torque). Bit-identical across two independent runs. 5350 s at -t8. |
| 10 | fig:LR_torque_and_reversal           | `fig10.jl` | `expected/fig10.json` | heavy | **preliminary** | **match** (structure) | T_{L,m} for (𝔭,h₀) ∈ {(5,0.03),(5,0.05),(0,0.03),(0,0.05)}, χ=0, 𝔮=0, a•=15M, b=0.6; plotted m=4:40, totals summed to ⌈4/h₀⌉ at the reference grid (240,240,24000). M:1014 PASSES (last m with T_{L,m}<0 = 10 at h₀=0.05 and 16 at h₀=0.03; total torque (74) = +1.539/+1.841 F_{J,0}h₀ > 0). FINDING: our T_{L,m} is ≈13 % below the published figure at low m (ratio 0.87 at m=4–10 → 1.02 at m=40); domain-extent and boundary-branch hypotheses refuted (evidence dir). ~28 min at -t8. |
| 11 | fig:convergence                      | `fig11.jl` | `expected/fig11.json` | heavy | **solid** | **match** | a•=25M, χ=0, h0=0.05, m=1:200, triples (30,3000)/(60,6000)/(120,12000) vs (240,240,24000). Median rel errors 1.8e-4 / 5.0e-6 / 9.9e-7. UNIT FINDING: the unlabelled absolute panel of the paper is in F_{J,0}, not F_{J,0}h0 — measured from a tick-calibrated pixel read (factor 1/h0 = 20 between the two panels' agreement). 798 s at -t4 for the reference triple. |
| 12 | fig:pole_disp_errors                 | `fig12.jl` | — | heavy | not started | — | |
| 13 | fig:variation_with_respect_to_b      | `fig13.jl` | `expected/fig13.json` | heavy | **solid** | **match** | T (74) in F_{J,0}h₀ vs b ∈ {0.3…0.8} for h₀ ∈ {0.05,0.075,0.1,0.125,0.15}, 𝔭=𝔮=χ=0, a•=100M, m = 1:⌈4/h₀⌉ on the (60,60,6000) triple (A:139) with an in-driver cross-check at the reference triple (rel 2–5e-6). A:164 PASSES: T(h₀=0.05) = 1.2610 (b=0.3) → 0.9714 (b=0.8), paper 1.28 → 0.97; all five series decrease monotonically. Pixel read-back of the reference: values within 1.7 % (median 0.3–1.1 %) and identical series ordering at every b. ~33 min at -t8. |
| 14 | fig:m95_vs_b                         | `fig14.jl` | `expected/fig14.json` | heavy | **solid** | **match** | m₉₅h₀ vs b ∈ {0.3…0.8}, 𝔭=𝔮=χ=0, a•=30M, h₀=0.05, m = 1:⌈4/h₀⌉ = 1:80 on (60,60,6000) with an in-driver reference-triple cross-check (same m₉₅). A:164 PASSES: 2.70 (b=0.3) → 1.90 (b=0.8), paper 2.9 → 1.95, monotone. Pixel read-back: four of six points differ by exactly one mode (0.05). FINDING: m₉₅ depends on the mode ceiling — at ⌈6/h₀⌉ the b=0.3 point is 2.900, the published value (registered as `m95_h0_Ksens_b03`). ~13 min at -t8. |
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
