# arXiv:2608.17003 — Dittmann, Hegade K. R. & Fairbairn, "Extreme-Mass-Ratio Inspirals in Gaseous Disks"

Mirrors: tex `ref-paper/arXiv-2608.17003/src/` · supplementary notebook
`ref-code/arXiv-2608.17003-anc/` (emri-disk repo).  Physics lives in
`capabilities/AccretionDisc` (`SpiralWaves` = the paper's linear theory; `Modes` = the
pressureless analytic baselines it is compared against).  Drivers here are THIN.

> **Provenance of the master system (2026-08-21, `EGD::jd-swap`).**  The relativistic master system
> (47) and the diagnostics integrands (72a)/(72b) that every driver below evaluates are **derived in
> Julia** by `capabilities/AccretionDisc/derivation/` (a standalone `Symbolics.jl` project that is
> *not* a package dependency) and emitted into `src/SpiralWaves/master_derived.jl` and
> `r2_derived.jl`.  **No external CAS — and no `wolframscript` — is in the build path of anything
> this directory reproduces.**  The authors' supplementary notebook (route R1) and the independent
> Mathematica derivation (route R2) survive as ORACLES in
> `capabilities/AccretionDisc/test/SpiralWaves/fixtures/oracle/`, gated forever on-shell by
> `test_derived_equivalence.jl`.  The goldens below are UNCHANGED by the swap: the derived and the
> notebook entries agree to Float64 roundoff (worst 2.6e-12 relative on the 50-point physical
> fixture, where the notebook form itself measures 7.0e-10).

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
| 1  | fig:spirals — δU maps                | `fig01.jl` | `expected/fig01.json` | heavy | **solid** | **match** | 2×3 δU M/m• maps, a•=15M, m=1:120, App.-C triple (120,120,12000), polar mesh r_ISCO…40M on 601² px. Winding compared with the published raster BY NUMBER (arm tracked in `spiral_grid.png` itself): turns(18M→35M) = −0.718/−1.084/−1.483 (χ=−0.9), −0.707/−1.067/−1.461 (χ=+0.9) vs published −0.721/−1.097/−1.504, −0.707/−1.081/−1.425 — worst \|Δ\| 0.036 turns. All three caption claims are driver gates and pass. Amplitude agrees to 4–6 % in the one unclipped column. 3936 s at -t8. |
| 2  | fig:spirals_conv — m_max convergence | `fig02.jl` | `expected/fig02.json` | heavy | **solid** | **match** | 𝔭=2, 𝔮=1, h0=0.1, a•=15M; m_max ∈ {10,30,120} as partial sums of ONE 120-mode solve per row. M:893 gate: ringing residual(10)=0.5095, residual(30)=0.1245 < residual(10)/3=0.1698, residual(40)=0.0615 < 0.10; ladder 0.746/0.510/0.250/0.125/0.061/0.015/0.002 monotone. Winding independent of the ceiling: turns(18M→35M) = −1.475/−1.482/−1.483 vs published −1.445/−1.478/−1.456. m_max=120 column reproduces fig01 panel 3 independently. 1801 s at -t8. |
| 3  | fig:angular_momentum_balance         | `fig03.jl` | `expected/fig03.json` | heavy | **solid** | **match** | a•=15M, χ=0, h0=0.05, (𝔭,𝔮)=(0,0), m=1:120. Balance (71) away from corotation 4.0e-5; the corotation spike is gated on the invariant ∫(dF_J/dr−d𝔗/dr)dr = 0.0266073 F_{J,0} = T_L−T of (75) to 2.7e-5 (its plotted height is a plot-sampling number). 605 s at -t4. |
| 4  | fig:dtdr_ap — torque density vs a•    | `fig04.jl` | `expected/fig04.json` | heavy | **solid** | **match** | χ=𝔭=𝔮=0, h0=0.05, a•∈{10,30,100,300,1000}M, m=1:120. All four M:919 numbers pass; pixel-calibrated reference agreement 0.00–0.25 % on the measurable extrema. 3779 s at -t4. |
| 5  | fig:TGR_vs_TNewt                     | `fig05.jl` | `expected/fig05.json` | heavy | **solid** | **match** | T_GR/T_Newt vs χ (10 points, −0.9…0.9) for the five (𝔭,𝔮) at h₀=0.05, a•=12M, b=0.6, m=1:80, grid (30,30,3000). T_Newt is the FULLY Newtonian pipeline (`SpiralWaves` `model = :newtonian`, EGD::sw-newtonian-limit) and is χ-independent by construction. Caption claim M:930-932 (𝔭<0 enhances, 𝔭≥0 suppresses) PASSES on all 50 markers. Against the published markers (read at 200 dpi, ±0.015): 46/50 within 0.05, mean |Δ| ≤ 0.034 per series. FINDING [OPEN, under investigation]: the (1.5,0) series drifts to −0.107 at χ=+0.9 (−12 %), not a T_Newt normalisation (the factor is not constant along the series) — see the evidence dir. This is **not** covered by E-F1: fig05's abscissa is χ, read off tick-calibrated markers, and its a• is the single fixed value 12M, so no abscissa misalignment is possible. It remains an open residual alongside fig10's low-m offset. Wave 6 CLOSED the boundary candidate at fig05's own configuration ((1.5,0), χ=+0.9, a•=12M, r_ISCO = 2.32088M — the deepest prograde ISCO in the paper's range): over four inner-boundary variants the modal torque moves by ≤ 1.14 % (m=4) and \|R_m\| ≤ 7.5e-04, against a −12 % deviation, so a boundary treatment is REFUTED here too and fig05 and fig10 are consistent with ONE model-level difference in the driven wave amplitude (`evidence/fig05-tgr-tnewt/BC_TRANSFER.md`). 3470 s at -t8; `--check` re-ran all 55 configurations independently and reproduced the golden **BITWISE** (37/37 values, max relative difference 0). |
| 6  | fig:torque_chi_ap                    | `fig06.jl` | `expected/fig06.json` | heavy | **solid** | **match** | T (F_{J,0}h₀) vs a•/M on the **PUBLISHED abscissae** (`_published_grids.jl`: integer a• = 11, 12, 13 M then geometric to ≈1063 M, 40 points) for χ ∈ {0.9, 0, −0.9}, panels (𝔭,𝔮) = (1.5,1) and (0,0), h₀=0.05, b=0.6, m=1:80, grid (30,30,3000); 240 configurations, 5719 s at -t8. M:925 and M:936 PASS; spin spread 2.0 % at a•=1063M. The recorded 5–15 % strong-field deficit is WITHDRAWN — it was an abscissa misalignment; agreement over the 15 measured published markers is median 0.995, range 0.978–1.024, worst 2.4 %. |
| 7  | fig:torque_ap_h0                     | `fig07.jl` | `expected/fig07.json` | heavy | **solid** | **match** | The same 40 PUBLISHED abscissae as fig06 (shared `_published_grids.jl`) for h₀ ∈ {0.05…0.15} at χ=0; 400 configurations, 5151 s at -t8. M:925 and M:947 PASS (plateau 1.1201, 1.1265, 1.1723, 1.2031, 1.1337 → peak at h₀=0.125). Cross-driver identity with fig06 now EXACT: max\|Δ\| = 0.000e+00 over the 40 points. |
| 8  | fig:m_95_plot                        | `fig08.jl` | `expected/fig08.json` | heavy | **solid** | **match** | Violins of m₉₅h₀ over h₀ ∈ {0.05…0.15} × χ ∈ {−0.9,0,0.9} for five (𝔭,𝔮) × three a•, at the **PRIMARY ⌈6/h₀⌉ ceiling** (wave-6 D4): 225 surveys, 15 660 mode solves, 118 min at -t8. FULL K=6 band **[1.95, 3.15]** — exactly the wave-5 two-configuration prediction, so a full K=6 survey does not widen it. M:1025 passes all three parts (containment; coverage 76.7 % vs the published 63.3 %; pixel cross-check within ±0.30 of [2.208, 3.157] — 0.258 / 0.007); M:992 spread 0.268. Neglected tail 1.57e-2 → 1.07e-3. K=4 band [1.80, 3.00] registered. |
| 9  | fig:num_vs_analytical                | `fig09.jl` | `expected/fig09.json` | heavy | **solid** | **match** | Analytic branch = `Modes.hkr_lindblad_torque` (arXiv:2510.03564 (18b)+(20) summed by (24)); numerical branch = `SpiralWaves.torque_survey` (Eq. 75); cutoff m₄₅ from our own mode sum. **FULL 19-point a• grid** over [12,100] M, m = 1:100, 38 configurations, 8245 s at -t8. ALL FIVE printed numbers PASS: m₄₅h₀ = 0.35 (paper 0.35), (0,0) ratio 1.4844 @ 12M (1.5) and 0.9235 @ 100M (0.92), (0,1) band [0.8737, **0.9573**] ([0.88, 0.96]) — the (0,1) maximum moved 0.9363 → 0.9573 because the full grid samples a• = 16.889M where the published bump sits (pixel read-back 0.9576). VLM **match**: 18 of 19 markers agree with a 400-dpi read-back of the published PDF to ≤ 0.0101 (blue) / 0.0034 (orange). FINDING: a• = 12M is 2.6 %/2.9 % low — the ISCO-proximity deviation the paper predicts at M:1002. Determinism: the wave-3 7-point grid is an exact subset and all 56 shared values are BITWISE identical across the two runs. m₄₅ registered BOTH ways (Lindblad 0.35/h₀ = the paper's number; total-torque 0.65/h₀), settling the M:991 ambiguity. `--check` on the full grid PASSES **23/23** values at rel err 0.00e+00. |
| 10 | fig:LR_torque_and_reversal           | `fig10.jl` | `expected/fig10.json` | heavy | **preliminary** | **match** (structure) | T_{L,m} for (𝔭,h₀) ∈ {(5,0.03),(5,0.05),(0,0.03),(0,0.05)}, χ=0, 𝔮=0, a•=15M, b=0.6, plotted m=4:40 on the reference triple (240,240,24000); totals summed to ⌈6/h₀⌉ since wave-6 D4 (the plotted per-mode array is independent of that ceiling and unchanged bitwise). M:1014 all PASS. PRELIMINARY for a characterised residual: our T_{L,m} runs **5–18 % below** a 200-dpi tick-calibrated read-back (a **2–9 % amplitude** deficit), largest at low m and at the smaller h₀ and NOT converging by m=40. SIX hypotheses refuted — read-back, domain, E-M4 branch, E-M5 grid wavenumber, the ISCO boundary (\|R_m\| ≤ 1e-3 over 104 cells) and the softening b. |
| 11 | fig:convergence                      | `fig11.jl` | `expected/fig11.json` | heavy | **solid** | **match** | a•=25M, χ=0, h0=0.05, m=1:200, triples (30,3000)/(60,6000)/(120,12000) vs (240,240,24000). Median rel errors 1.8e-4 / 5.0e-6 / 9.9e-7. UNIT FINDING: the unlabelled absolute panel of the paper is in F_{J,0}, not F_{J,0}h0 — measured from a tick-calibrated pixel read (factor 1/h0 = 20 between the two panels' agreement). 798 s at -t4 for the reference triple. |
| 12 | fig:pole_disp_errors                 | `fig12.jl` | `expected/fig12.json` | heavy | **solid** | **match** | a•=100M, χ=𝔭=𝔮=0, ε_p∈{1e-3,1e-4,1e-5} vs 1e-6. O(ε_p) reproduced: median ratios 9.998 / 10.985 (reference pixel-read 9.94 / 10.95). FINDING: the n_C=2400 grid printed at A:162 resolves the ε_p=1e-6 pole only to m≈5, so the errors saturate there; the driver uses the paper's reference n_C=24000 (A:140) and its own cutoff m_max=⌈3/h0⌉ (M:105) — `GRTOOLS_FIG12_NC=2400` restores A:162. 3227 s at -t4. |
| 13 | fig:variation_with_respect_to_b      | `fig13.jl` | `expected/fig13.json` | heavy | **solid** | **match** | T (74) in F_{J,0}h₀ vs b ∈ {0.3…0.8} for h₀ ∈ {0.05…0.15}, at the **⌈6/h₀⌉ ceiling**: A:164 reads 1.2766 → 0.9724 against the published 1.28 → 0.97 (1.2610 → 0.9714 at ⌈4/h₀⌉, i.e. 1.5 % → 0.34 % from the published value); neglected tail 1.39e-2 → 1.3e-3; both K=4 values registered as `torque_h005_b0*_Ksens`. Reference-triple cross-check 2.4e-06 / 4.9e-06 relative. |
| 14 | fig:m95_vs_b                         | `fig14.jl` | `expected/fig14.json` | heavy | **solid** | **match** | m₉₅h₀ vs b ∈ {0.3…0.8}, 𝔭=𝔮=χ=0, a•=30M, h₀=0.05, at the **⌈6/h₀⌉ ceiling**: A:164 is hit EXACTLY — 2.900 (b=0.3) and 1.950 (b=0.8) against the published 2.9 and 1.95 (2.700 / 1.900 at ⌈4/h₀⌉); curve 2.900, 2.650, 2.450, 2.250, 2.100, 1.950, monotone. Neglected tail 1.47e-2 → 1.56e-3; reference triple gives the same m₉₅ = 58, 39. `--check` PASSES 22/22 bitwise. |
| —  | printed-number gate                  | `check_printed_numbers.jl` | none (by design) | fast | **preliminary** (18 PASS / 0 FAIL / 0 uncovered, exit 0) | n/a | asserts every printed number of the paper that the figure goldens cover; all 18 claims are covered and pass. `M1025-m95h0-band` was CORRECTED in wave 5 (the old criterion was falsified by the published figure itself) and its pixel cross-check RE-POINTED in wave 6 at the PRIMARY full-survey band, now the ⌈6/h₀⌉ one — [1.95, 3.15] against the published [2.208, 3.157], within ±0.30 on each edge, with [1.80, 3.00] (K = 4) printed beside it; that discharges the wave-5 scope caveat. `A164-m95-vs-b` now reads 2.9 → 1.95, the printed values exactly. Promotion to solid is refused by ONE preliminary predecessor, fig10. |

> **Table currency (2026-08-20, wave 6).** The `status` column mirrors the mission knowledge ledger
> (`results/ledgers/knowledge/paper_EMRI-GASDISK/`, emri-disk): **13 of the 14 figures are `solid`**
> (fig08 was promoted when the full ⌈6/h₀⌉ re-survey landed), fig10 is `preliminary` for a measured,
> bounded residual, and the printed-number gate is `preliminary` only because the ledger refuses to
> promote a node above its predecessors — fig10 is now its only preliminary predecessor. The
> fig06/fig07 rows and the `MANIFEST.toml` notes are BOTH current: the drivers were re-run on the
> published abscissae measured by erratum **E-F1**, which retired the strong-field "5–15 % deficit",
> and the manifest notes were rewritten with them.

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
| `M1025-m95h0-band`        | M:1025 | text | `fig08.m95_h0_min`, `fig08.m95_h0_max`, `fig08.m95_h0_min_Ksens`, `fig08.m95_h0_max_Ksens` | ⊆[2.0,3.5]±0.25, ≥55 % cover, the PRIMARY (⌈6/h₀⌉ full-survey) band within ±0.30 of the published pixels [2.208, 3.157]; the ⌈4/h₀⌉ band printed beside it |
| `M992-m95-inverse-h0`     | M:992  | text | `fig08.m95_h0_spread_fraction` | ≤ 0.35 |
| `M999-m45-cutoff`         | M:999  | text | `fig09.m45_times_h0` (Lindblad reading; the total-torque reading `fig09.m45_T_times_h0` = 0.65 is registered but NOT asserted — see the fig09 row) | 0.35 (±0.1) |
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

### `M1025-m95h0-band` — why the criterion changed (wave 5)

The first version of this claim required `m95_h0_min ∈ [1.75, 2.25]` **and**
`m95_h0_max ∈ [3.25, 3.75]`, i.e. that our extremes sit *at* the edges of the paper's rounded
prose interval "m₉₅h₀ ≃ 2.0–3.5" (M:1025).  A pixel read-back of the **published** Fig. 8 —
recorded while our own `fig08.jl` survey was still running, so it cannot be a post-hoc
accommodation (evidence: `results/emri-gasdisk/paper_EMRI-GASDISK/evidence/fig08-m95-violins/`
`REFERENCE_READBACK.md` + `ref_fig08_readout.json` in the emri-disk repo) — measures the published
band as **[2.208, 3.157]** across its fourteen resolved violins.  *The published figure fails the
old criterion* (3.157 < 3.25): the prose interval is a rounded envelope, not a pair of attained
extremes.  A criterion the paper's own figure violates cannot judge a reproduction of it.

The claim now asserts three things, all visible in the gate row:

1. **containment** — our band inside [2.0, 3.5] widened by ±0.25 (the rounding scale the read-back
   exhibits: the published band sits 0.208 above the lower edge and 0.343 below the upper one);
2. **coverage** — the overlap of our band with [2.0, 3.5] is ≥ 55 % of its width; the published
   band covers 63.3 %, so the threshold is one the reference passes with margin;
3. **pixel cross-check** — the PRIMARY band `m95_h0_min/max` within ±0.30 of [2.208, 3.157].

Part 3 compares like with like: the published read-back is the band of the *published full survey*,
so it is held against *our full survey* at the primary ceiling.  **Wave 6** (decision D4 of
`progress/emri-gasdisk/DECISIONS_W6.md`, resolving D2 of `DECISIONS_W5.md`) made ⌈6/h₀⌉ the primary
ceiling of `fig08`, `fig10`, `fig13` and `fig14`: fig08's own ceiling scan, fig14 (m₉₅h₀ at b = 0.3:
2.700 → 2.900) and fig13 (a 0.3–1.7 % torque deficit) independently indicate the paper summed
≈ 6/h₀ azimuthal modes, and the low-m deficit investigation
(`results/emri-gasdisk/paper_EMRI-GASDISK/evidence/fig10-per-mode-reversal/`, emri-disk) established
that the ceiling and that deficit are INDEPENDENT — a per-mode `T_{L,m}` does not depend on where the
mode SUM is truncated.  The pre-registered ⌈4/h₀⌉ band stays registered as `m95_h0_*_Ksens` and is
printed beside the primary band, so the ceiling difference is never absorbed silently.

The wave-5 **scope caveat is discharged**: part 3 is no longer a two-configuration comparison — it
now uses the extremes of the full ⌈6/h₀⌉ survey.
