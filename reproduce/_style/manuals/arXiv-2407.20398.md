# FIGURE STYLE MANUAL — arXiv:2407.20398 (Fairbairn & Rafikov 2024/25, "Eccentric planet-disc interactions: orbital migration and eccentricity evolution", arXiv v1 = submitted text; MNRAS 537, 1779) — FR24

VLM-inspected from the source figure PDFs in
`/data/haiyangw/claude/emri-disk/ref-paper/arXiv-2407.20398/src/Figures/`,
rendered at 200 dpi, and cross-checked against the compiled PDF
`2407.20398.pdf` (pages 8–25 at 150 dpi). The 17 published figures use these
files:

| Fig | file |
|---|---|
| 1 | `circular_morphology_torques.pdf` |
| 2 | `analytical_theory_torque.pdf` |
| 3 | `circular_torques.pdf` |
| 4 | `time_dependent_torque_q_0.0_p_0.0.pdf` |
| 5 | `torque_density_integrated_torque.pdf` |
| 6 | `non_coorbital_features.pdf` |
| 7 | `torque_eccentricity_h_0.06.pdf` |
| 8 | `eccentric_timescales_h_0.06.pdf` |
| 9 | `1d_migration.pdf` |
| 10 | `a_e_tracks.pdf` |
| 11 | `tau_comparison.pdf` |
| A1 | `convergence.pdf` |
| B1 | `circular_torque_q1.pdf` |
| C1 | `torque_eccentricity_h_0.08.pdf` |
| C2 | `torque_eccentricity_h_0.1.pdf` |
| C3 | `eccentric_timescales_h_0.08.pdf` |
| C4 | `eccentric_timescales_h_0.1.pdf` |

All other files in `src/Figures/` (`Figure_Archive/`, `tau_m_comparison.pdf`,
`torque_eccentricity.pdf`, …) are NOT used by `main.tex`. Captions are from
`src/Sections/*.tex`. Mission brief:
`emri-disk/results/emri-eccentric/paper_EMRI-ECC/evidence/acquire-2407.20398/BRIEF.md`.

Tags: [SOLID] = file metadata / pdffonts / measured pixels (colormaps matched
by sampling the colourbar against matplotlib colormaps; hex codes = dominant
saturated pixels). [PRELIMINARY] = estimated by eye.

READ THE RELEVANT SECTION BEFORE GENERATING ANY REPRODUCTION FIGURE; then
VLM-compare per the `vlm-figure-repro` skill.

---

## 0. Shared style conventions (whole paper)

**Toolchain [SOLID]:** every figure is **matplotlib 3.8.0** PDF with
**`text.usetex=True` and the DEFAULT sans-serif family**. `pdffonts` shows CMSS
(CM Sans) together with CMR/CMMI/CMSY. As a result:
- plain text in labels, legends and panel letters is **CM Sans**: "and" in
  `T and F_J`, `(a)`, `Numerical`, `Analytical`, `Papaloizou and Larwood`;
- math and tick labels are **CM serif** (CMR/CMMI).

Reproduce this with `text.usetex=True` and do NOT set `font.family='serif'`.
This mix of a sans word inside a serif label is the most recognisable FR24
tic.

**Axes [SOLID]:** matplotlib default ticks: **outward**, left and bottom only,
no minor ticks except on log axes, full black box, no grid. Exceptions:
- the δΣ maps in Fig 1 show only the left and bottom spines;
- Figs 9 and A1 add a secondary TOP x axis (t/yr; |l−m|).

**Palette [SOLID, measured]:** named matplotlib colours: black `'k'`, red
`'r'` #ff0000, blue `'b'` #0000ff, green `'g'` #008000, `'orange'` #ffa500.
The fixed SEMANTIC triple in Figs 3, 5c, 7, 8, C1–C4 is:
- **black = T_net** (or τ_L⁻¹);
- **red = T_L** (or τ_a⁻¹);
- **blue = T_C** (or τ_e⁻¹/100).

Fig 4 alone uses the default tab10 cycle (#1f77b4, #ff7f0e, #2ca02c, #d62728).
Reference lines (zero line, transonic e = h_p) are **dashed 'tab:gray'
#7f7f7f**. Fig 6's resonance markers are red at alpha 0.5 (#ff7f7f).

**Colormaps [SOLID]:**

| use | colormap | limits |
|---|---|---|
| δΣ maps (Fig 1) | viridis | ±2.0 |
| e-family of curves (Fig 5) | jet | 0–0.12 |
| track time (Fig 10) | inferno | 0–350 |
| error heat map (A1) | viridis | 0–0.87 |

Sampling distances: viridis 0.0044 / 0.0054, jet 0.0116, inferno 0.0080
(mean ΔRGB; the next best candidates are ≥ 0.16).

**Legends [SOLID]:** matplotlib default framed legend (rounded, white face,
grey edge, framealpha 0.8). Figs 7, 8 and C1–C4 use one figure-level legend
above the grid (`ncol=3`). Fig B1 has NO legend: its colours are defined only
in the caption.

**Panel letters:** `(a)`, `(b)`, … in CM Sans, inside the axes, with no box.
- Lower-left: Figs 5, 7, C1, C2, 9.
- Upper-left: Figs 8, C3, C4.
- Upper-right: Fig 11.
- Above the map: Fig 1.

**Units:**
- torques in `F_{J,0}`; torque density in `F_{J,0}/a`;
- timescales as `\tau^{-1}/\tau_0^{-1}`;
- τ_e is plotted divided by 100 in the 4×4 grids (Figs 8, C3, C4) but NOT in
  Fig 11.

**Figure sizes [SOLID, PDF page boxes; inches]:**

| figure(s) | size (in) | note |
|---|---|---|
| 1 | 9.97 × 7.17 | |
| 2 | 4.82 × 3.79 | |
| 3 | 10.7 × 6.4 | |
| 4 | 4.9 × 3.4 | |
| 5 | 12.3 × 3.3 | |
| 6 | 6.8 × 4.7 | |
| 7, 8, C1–C4 | 10.5–10.6 × 13.96 | included at `0.9\textwidth` |
| 9 | 5.0 × 6.5 | |
| 10 | 4.9 × 4.9 | |
| 11 | 5.1 × 6.2 | |
| A1 | 5.0 × 8.0 | exact `figsize=(5, 8)`, no tight bbox |
| B1 | 5.9 × 10.0 | |

Single-column figures are included at `\columnwidth`, so their fonts appear
smaller than nominal. Font sizes [PRELIMINARY]: ticks ≈ 12 pt, labels ≈ 14–17
pt (`CMR17`/`CMSS17` appear in the big grids).

Reference rcParams block:
```python
import matplotlib as mpl
mpl.rcParams.update({
    'text.usetex': True,                         # font.family left at default 'sans-serif'
    'font.size': 12, 'axes.labelsize': 16,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'axes.grid': False, 'legend.frameon': True, 'legend.fancybox': True,
})
TRIPLE = {'net': 'k', 'L': 'r', 'C': 'b'}        # semantic, never permute
REFLINE = dict(color='tab:gray', ls='--')
```

---

## 1. Plot type A — circular δΣ maps over radial T/F_J profiles (Fig 1)

**Layout [SOLID]:**
- 2 rows × 3 columns: (a) (q, p) = (0, 0), (b) (0, 1.5), (c) (0.5, 1.5); all
  at h_p = 0.06.
- **ONE horizontal viridis colourbar** across the top, ±2.0, ticks every 0.5,
  labelled ABOVE the bar `$(M_\mathrm{p}/M_\mathrm{th})^{-1}\delta\Sigma/\Sigma$`.
- Panel letters `(a)`, `(b)`, `(c)` in sans, centred above the maps.

**Top row (maps):**
- **Field:** a rasterised square-aspect δΣ/Σ map on X, Y ∈ [−2, 2] (ticks
  −2 … 2).
- **Labels:** serif italic `$X$` under every map; `$Y$` on the first map
  only.
- **Axes:** left and bottom spines only, ticks out.
- **Disc:** r = 0.05 … 2. The disc edge is visibly pixelated and the field
  outside r = 2 is white.
- **Planet:** a black filled dot at (1, 0).

**Bottom row (profiles):**
- **Curves:** `$T$` (one-sided torque, black solid, lw ≈ 2) and `$F_J$` (red
  dashed, lw ≈ 2) vs `$r/a$` ∈ [0, 3].
- **Axes:**
  - y: label `$T$ and $F_J$ [$F_{J,0}$]` (first column only); range 0–0.8
    in (a) and (b), 0–1.5 in (c);
  - full box.
- **Insets** in (a) and (b) only, lower right, thin frame:
  - x 0.95–1.05 (ticks 0.95, 1.00, 1.05), y 0–0.2;
  - small serif ticks.

**Dominant visual STRUCTURE:**
- **Maps:** ONE thin bright trailing arm per map, seen as a fine yellow line
  on a teal background.
  - The outer arm runs from the planet down and clockwise to r ≈ 1.8.
  - The inner arm wraps counter-clockwise inward to a ring at r ≈ 0.5.
  - Inside r ≈ 0.4 there are tight noisy concentric rings around a white
    central hole.
  - The three maps are visually almost identical.
- **Profiles:**
  - Both curves are V-shaped notches to 0 at r = a, with flat plateaus
    elsewhere.
  - (a): T plateaus 0.57 (inner) / 0.65 (outer); F_J plateaus 0.58 / 0.74.
    The inset shows **F_J jumping discontinuously at r = 1.00** (from 0 to
    ≈ 0.09: the corotation torque) while T is a smooth parabola.
  - (b): F_J lies above T everywhere. In the inset the F_J minimum is ≈ 0.045,
    above T's zero.
  - (c): F_J (red) diverges as r → 0 and decays outward from 0.75 to 0.47
    (the AMF is NOT conserved for q ≠ 0). T is flat at 0.51 / 0.73.
- **Consistency check:** outer minus inner plateau of T = T_net (0.08 in
  (a), matching Fig 3).

---

## 2. Plot type B — per-mode spectrum, log–log (Fig 2)

- **Axes:** x = `$m$`, log, m ∈ [2, 150] (labelled ticks 10¹, 10²); y =
  `$T_\mathrm{L}/F_{J,0}$`, log, 10⁻⁷ … 3×10⁻², with decade labels.
- **Series:** `Numerical` = black solid; `Analytical` = red dashed (PL00 / A93
  WKB formula, eq. 42). Framed legend, upper right.
- **Structure:**
  - **Black:** 1.15e-2 at m = 2, peak ≈ 1.4e-2 at m ≈ 3–4, 5.5e-3 at 10,
    ≈ 1e-4 at 70, ≈ 2.5e-7 at 150. It wiggles slightly at m > 100.
  - **Red:** peaks higher (≈ 2e-2 at m ≈ 5) and lies above black for
    3 ≲ m ≲ 85, then crosses BELOW black at m ≈ 85–90 (≈ 1.3e-7 at 150).
  - Both curves turn over steeply (the cut-off) beyond m ≈ 30.

---

## 3. Plot type C — 3 × 3 parameter-slice scan (Fig 3)

**Layout:**
- **Rows:** `$T_\mathrm{net}/F_{J,0}$`, `$T_\mathrm{L}/F_{J,0}$`,
  `$T_\mathrm{C}/F_{J,0}$` (y labels on the left column).
- **Columns:** x = `$q$` (0–0.75), `$p$` (0–1.5), `$h_\mathrm{p}$`
  (0.06–0.10), with x labels on the bottom row only. Column titles:
  `$p = 0.0,\ h_\mathrm{p} = 0.06$`, `$q = 0.0,\ h_\mathrm{p} = 0.06$`,
  `$q = 0.0,\ p = 0.0$`.
- Every panel has its own y range, and the x limits are tight to the data,
  so the end markers are half-clipped.

**Series [SOLID]:**
- **black** filled circles joined by a solid line (`'ko-'`, ms ≈ 8) = with
  m = 1;
- **red** crosses `'rx'` with no line = without m = 1;
- **blue dashed** (#0000ff) = the hyperplane fits, eqs (43)–(44).
- There is no legend; the encoding is given in the caption.

**Structure:**
- All rows are near-linear in each parameter.
- **T_net row:** red crosses lie ABOVE black (e.g. at q = 0: 0.087 vs 0.078).
- **T_L row:** red lies BELOW black.
- **T_C row:** values are negative; red is less negative than black. In the
  p column T_C rises linearly from −0.082 to 0.000 at p = 1.5 (zero
  vortensity gradient), and red reaches 0 there too.
- Blue fit lines hug the black points, except for slight offsets in the T_L
  vs p column (fit above data at p = 1.5) and the h column (fit below at
  0.10).

---

## 4. Plot type D — radial torque-density / AMF profiles (Figs 4, 5, 6, B1)

**Fig 4 — instantaneous dT/dr at 4 phases (e = 0.01, q = p = 0):**
- **Axes:** x `$r/a$` 0.8–1.2; y `$dT/dr\quad[F_{J,0}/a]$` (a wide space
  before the bracket), range ±43, ticks ±40, ±20, 0.
- **Series:** tab10 solid lines for `$0.0T$` (C0 blue), `$0.25T$` (C1
  orange), `$0.5T$` (C2 green), `$0.75T$` (C3 red), plus a black dashed
  `$\langle dT/dr\rangle$`.
- **Reference line:** a tab:gray dashed horizontal line at 0.
- **Planet markers:** coloured filled dots ON the zero line at r_p(t): blue
  0.99, orange/red ≈ 1.00 (overlapping), green 1.01.
- **Legend:** framed, upper left.
- **Structure:**
  - orange (0.25T): trough −39 at 0.99, crest +39 at 1.01;
  - red (0.75T): mirror image, +31 at 0.99 and −32 at 1.01;
  - blue and green: weak (±8);
  - the orbit average (black dashed) is the familiar small antisymmetric
    S-curve: −7 at 0.93, +7 at 1.05.

**Fig 5 — e-family (h_p = 0.06, q = p = 0), 1 × 3 + colourbar:**
- **Titles** above the panels (no y labels): (a) `$dT/dr$ [$F_{J,0}/a$]`, (b)
  `$T(r)$ [$F_{J,0}$]`, (c) `$T$ [$F_{J,0}$]`.
- **Panels (a) and (b):**
  - x `$r/a$` 0.5–1.5;
  - 13 lines coloured by **jet(e/0.12)** for e = 0, 0.01, …, 0.12;
  - vertical colourbar between (b) and (c), titled `$e$` on top, ticks
    0.000–0.100 step 0.025.
- **Panel (c):**
  - x `$e$` 0–0.12;
  - `'o-'` black `$T_\mathrm{net}$`, red `$T_\mathrm{LB}$` (sic), blue
    `$T_\mathrm{C}$`;
  - tab:gray dashed zero line;
  - legend at the left middle, small.
- Italic panel letters `$(a)$`, `$(b)$`, `$(c)$` at lower left.
- **Structure:**
  - **(a):** low e (blue) = standard pair, trough −7 at 0.95 and crest +7 at
    1.05. High e (dark red) = REVERSED sign: crest +13.5 at 0.90 and trough
    −24 at 1.10. Curves pivot through 0 at r = 1.
  - **(b):** low-e plateaus are positive on both sides (≈ 0.55 / 0.65). High-e
    plateaus are negative (≈ −0.7 inner, −1.55 outer). All curves pass
    through 0 at r = 1 (a V or Λ shape).
  - **(c):** black falls from 0.08 to 0 at e ≈ 0.06, then to −0.84. Red runs
    from 0.16 to −0.79. Blue is nearly flat (−0.08 to −0.05). Blue crosses
    red and black at e ≈ 0.065–0.07.

**Fig 6 — F_J(r) at e = 0.06 with resonance zooms:**
- **Main axes:** x `$r/a$` 0.3–1.7; y `$F_J/F_{J,0}$` 1.40–1.62; one black
  line.
- **Resonance markers:** three vertical dashed light-red (#ff7f7f) lines at
  r = 0.76314, 1.0, 1.31037.
- **Insets** (thin frames, small serif ticks):
  - lower left: x 0.75814–0.76814 (ticks at exactly those three values), y
    with offset text `+1.546` and ticks 0.0002–0.0006;
  - upper right: x 1.30537–1.31537, y 1.6140–1.6141.
- **Structure:**
  - plateau 1.545 inside, a double trough (1.41 at 0.93, 1.415 at 1.07);
  - a vertical **jump at r = 1** (≈ 1.505 → 1.563);
  - plateau ≈ 1.61 outside;
  - each inset shows a small upward step at its dashed line (≈ 1e-4).

**Fig B1 — m = 1 response (h_p = 0.06, b = 0.3, p = 0, e = 0):**
- **Layout:** two stacked panels, **hspace = 0**, shared x `$r/a$` 0–5
  (ticks 1 … 5). The top panel is `$F_J/F_{J,0}$` (−0.015 … 0.013); the
  bottom is `$T/F_{J,0}$` (−0.021 … 0.011).
- **Series:** red = q = 1.0, black = 1.005, green #008000 = 0.995 (black and
  green almost overlap). There is NO legend.
- **Reference line:** a grey dashed vertical line at the OLR, r = 2^{2/3} ≈
  1.587.
- **Structure:**
  - a sharp spike pair at r = 1 (+0.0068 / −0.014 in F_J);
  - red F_J has a DISCONTINUOUS jump at the OLR (−0.0007 → 0.0085);
  - all curves show decaying oscillations outward (≈ 13 wiggles to r = 5);
  - in the T panel red settles at −0.0165, while black and green oscillate
    about −0.003.

---

## 5. Plot type E — 4 × 4 (q, p) grids vs eccentricity (Figs 7, 8, C1–C4)

**Layout [SOLID]:**
- 16 panels ≈ square, wspace ≈ 0.3, hspace ≈ 0.2 [PRELIMINARY].
- Columns: `$p = 0.0$`, `0.5`, `1.0`, `1.5` as titles. Rows: `$q = 0.0$`,
  `0.25`, `0.5`, `0.75`, written rotated on the RIGHT edge.
- Panel letters (a)–(p) row-major, sans.
- y label `$T$` (Figs 7, C1, C2) or `$\tau^{-1}$` (Figs 8, C3, C4) on the left
  column only; x `$e$` on the bottom row only (0–0.12, ticks 0.00, 0.05,
  0.10).
- Every panel has an independent y range.
- One framed figure legend centred above the grid, 3 columns:
  - torques: `$T_\mathrm{net}$` black, `$T_\mathrm{L}$` red, `$T_\mathrm{C}$`
    blue;
  - timescales: `$\tau_L^{-1}/\tau_0^{-1}$` black, `$\tau_a^{-1}/\tau_0^{-1}$`
    red, `$\tau_e^{-1}/(100\tau_0^{-1})$` blue.
- Series are `'o-'` (ms ≈ 6, lw ≈ 1.5).
- In every panel:
  - a tab:gray dashed vertical line at the transonic **e = h_p** (0.06 in
    Figs 7, 8; 0.08 in C1, C3; 0.10 in C2, C4);
  - a tab:gray dashed horizontal line at 0.
- **Data:** 13 points, e = 0 … 0.12. The τ_e series (blue, timescale figures)
  starts at e = 0.01, because the e = 0 entry is undefined (it is stored as 0
  in the cube).

**Dominant visual STRUCTURE (Fig 7, h_p = 0.06):**
- T_net and T_L stay positive and flat-to-rising up to e ≈ h_p, then plunge
  steeply negative. They reach −0.5 … −0.9 by e = 0.12.
- The zero crossing moves to larger e with p: ≈ 0.06 (p = 0), 0.075 (0.5),
  0.08 (1.0), 0.085 (1.5).
- For p ≥ 1 there is a hump peaking at e ≈ 0.06–0.07.
- T_C (blue) is small and negative, with a slight upward trend in e. It is
  **≡ 0 at (q, p) = (0, 1.5)**, panel (d), where black is hidden under red.
- In C1/C2 the same patterns are shifted to the larger transonic e. In C2
  (h_p = 0.10) T_net crosses zero only at e ≈ 0.10–0.12.

**Dominant visual STRUCTURE (Fig 8):**
- Black τ_L⁻¹ repeats the T_net shape.
- Red τ_a⁻¹ rises to a sharp peak at e ≈ 0.07 (0.66 for p = 0, up to 1.73
  for p = 1.5 / q = 0.75), then falls (it goes negative at e ≈ 0.11–0.12 for
  p ≤ 0.5).
- Blue τ_e⁻¹/100 is **nearly identical in all 16 panels**: 0.26 at e = 0.01,
  rising to **1.00 at e = 0.08**, then 0.51 at 0.12. This is genuine: the
  values are in `src/anc/eccentric/tau_e_inv.out`, block (q, p) = (0, 0),
  h = 0.06, peak 100.2 at e = 0.08. It is NOT a plotting artefact.
- In C3 (h = 0.08) the blue peak is ≈ 0.73 at e ≈ 0.11.
- In C4 (h = 0.10) red and blue rise monotonically to the edge, with no peak
  inside 0.12. Black is flat, then turns down after e ≈ 0.08–0.10. In the
  p = 0 column it crosses zero at e ≈ 0.10–0.11 and reaches ≈ −0.24 at 0.12.

---

## 6. Plot type F — orbital evolution (Figs 9, 10)

**Fig 9 — migration time series:**
- **Layout:** two stacked panels, hspace = 0, shared x `$t/T_0$` 0–800
  (bottom). There is a secondary top axis `$t/\mathrm{yr}$` (0–9000, ticks
  0 … 8000) on the top panel.
- **Axes:** (a) `$a/a_0$` 0.96–1.00; (b) `$e/e_0$` 0–1.
- **Series:** solid lines, black `$q = 0.0,\ p = 0.0$`, red `$q = 0.0,\
  p = 1.5$`, blue `$q = 0.75,\ p = 0.0$`. Framed legend in panel (b), upper
  right.
- **Panel letters:** sans (a)/(b), lower left.
- **Structure:**
  - **(a):** black and blue first rise slightly above 1 (to ≈ 1.0003 at
    t ≈ 30), then decline. Red declines monotonically and fastest.
    Values at 800 T_0: black 0.985, blue 0.978, red 0.9615.
  - **(b):** all three e/e_0 curves are almost indistinguishable: a steep
    drop to 0.5 by t ≈ 110, then ≈ 0.1 at 800.

**Fig 10 — (a, e) tracks:**
- **Axes:** x `$a/a_0$` 0.995–1.0005; y `$e$` 0.008–0.125.
- **Tracks:** six tracks for e_0 = 0.02, 0.04, …, 0.12, drawn as dense
  overlapping **scatter markers** (large 'o', no edge). Their colour is
  **inferno(t/350)**: they start BLACK at a/a_0 = 1 and turn yellow late.
- **Colourbar:** horizontal, on top, label `$t/T_{a_0}$` above the bar, ticks
  0 … 350 step 50.
- **Reference line:** a grey dashed vertical line at a/a_0 = 1.
- **Inset** (upper left, framed): `$t_{a=0.995}/T_{a_0}$` vs `$e_0$`
  (0.02–0.12), black line.
- **Structure:**
  - all tracks move left (inward) while e decays;
  - the e_0 = 0.12 track first hooks RIGHT (to a/a_0 ≈ 1.0001), then joins
    the e_0 = 0.10 track at e ≈ 0.09;
  - end points at a = 0.995: e ≈ 0.047 (the merged e_0 = 0.10/0.12 pair),
    0.040, 0.030, 0.020, 0.010, so five visible ends for six starts;
  - end colour runs from pink (≈ 150, high e_0) to pale yellow (≈ 350,
    e_0 = 0.02);
  - the inset falls from ≈ 365 at e_0 = 0.02 to a minimum ≈ 135 at
    e_0 ≈ 0.095, then rises to ≈ 170.

---

## 7. Plot type G — comparison with published prescriptions (Fig 11)

- **Case:** (q, p, h_p, b) = (1.0, 1.5, 0.05, 0.5).
- **Layout:** two stacked panels, hspace = 0, shared x `$e$` 0.01–0.12.
  Large sans panel letters (a)/(b) at upper right.
- **(a)** `$\tau_L^{-1}/\tau_0^{-1}$` (−0.4 … 0.4):
  - linear theory = black `'ko-'` (ms ≈ 9), with NO legend entry;
  - PL00 eq (52) = green dashed `Papaloizou and Larwood`;
  - CN06 eq (54) = blue dashed `Cresswell and Nelson`;
  - framed legend, upper right.
- **(b)** `$\tau_e^{-1}/\tau_0^{-1}$` (0–160):
  - black `'ko-'`;
  - green dashed PL00 (53);
  - CN06 (55) dashed, red / blue / orange #ffa500 =
    `Cresswell and Nelson $Q_e = 0.1$` / `0.2` / `1.0`;
  - framed legend, upper right.
- **Structure:**
  - **(a):** black rises 0.145 → 0.168 (e = 0.04), crosses 0 at e ≈ 0.07,
    reaches −0.205 at 0.12. Green is much steeper (0.37 → −0.38 at 0.095,
    crossing black at e ≈ 0.045). Blue is shallower (0.13 → −0.135),
    crossing black at e ≈ 0.088.
  - **(b):** black is the only curve that RISES at low e (21 → 46 at
    e = 0.07). It then decays (18 at 0.12) and merges with green and blue
    for e ≥ 0.07. The prescriptions all decrease monotonically, in the order
    red (156) > green (85) > blue (78) > orange (15.5) at e = 0.01.

---

## 8. Plot type H — convergence heat map + slices (Fig A1)

- **Size and layout:** portrait 5 × 8 in.
  - Top: a horizontal viridis colourbar (0–0.87, ticks 0.0 … 0.8, label
    `$T_\mathrm{error}$` BELOW the bar).
  - Under it: the heat map, `pcolormesh` over x `$m$` (≈ 1–175) and y
    `$|l-m|$` (≈ 0–40), on a non-uniform grid.
- **Grid lines:** faint light cell-edge lines are visible (`edgecolors`
  ≈ white, low alpha).
- **Slice markers:**
  - a red dashed horizontal line at |l−m| = 40 (top edge);
  - a blue dashed vertical line at m = 175 (right edge).
- **Bottom axes:** `$T_\mathrm{error}$` (−0.05 … 0.92) vs `$m$` (bottom axis),
  plus a twin TOP axis `$|l-m|$` (5 … 40).
  - Red solid = slice in m at |l−m| = 40.
  - Blue solid = slice in |l−m| at m = 175.
- **Structure:**
  - **Heat map:** yellow (≈ 0.87) along the left edge (m ≲ 10) and the
    bottom rows (|l−m| ≲ 2); a diagonal green/teal transition; dark purple
    (→ 0) in the upper-right corner.
  - **Slices:** red has a small initial bump (0.85 → 0.88 at m ≈ 5), then a
    smooth decay to 0 at 175. Blue decays faster and is ≈ 0 by |l−m| ≈ 30
    (top axis).

---

## 9. Per-figure table

Provenance: **LIN** = linear theory (the target of reproduction). **FIT** =
closed-form fit or prescription evaluated from the paper's equations. **DATA**
= also tabulated in `src/anc/` (golden comparison only). **DER** = derived from
the LIN data cube (interpolation + ODE). FR24 contains NO hydro and NO
digitised overlays.

| Fig | panels | plot type | what is plotted | provenance per series | structure to check |
|---|---|---|---|---|---|
| 1 | 2×3 + top cbar | A | δΣ/Σ maps (q,p) = (0,0), (0,1.5), (0.5,1.5), h = .06; T(r), F_J(r) + coorbital insets | LIN | one thin trailing arm; F_J jump at r = 1 in inset (a); F_J divergence at small r in (c) |
| 2 | 1 | B | T_L per m vs m, (q,p,h,b) = (0,0,.06,.3) | black LIN; red FIT (eq 42; shifted LR positions not given) | red above black for 3 ≲ m ≲ 85, crossing near 85–90 |
| 3 | 3×3 | C | T_net, T_L, T_C vs q, p, h_p | black LIN (with m = 1), red × LIN (no m = 1), blue dashed FIT (43–44); DATA in anc/circular | linear trends; T_C(p = 1.5) = 0; red above/below black by row |
| 4 | 1 | D | instantaneous dT/dr at 0, .25, .5, .75 T + orbit average, e = .01 | LIN | ±39 antisymmetric orange/red pair; dots at r_p on zero line |
| 5 | 1×3 + cbar | D (+E-like c) | (a) dT/dr, (b) T(r) for e = 0 … .12 (jet); (c) T_net, T_LB, T_C vs e | LIN; (c) DATA | sign reversal of dT/dr pattern at high e; plateaus flip sign; (c) crossing ≈ .065 |
| 6 | 1 + 2 insets | D | F_J(r), e = .06, resonances (m,l) = (2,3), (3,2) | LIN | dashed at 0.76314, 1, 1.31037; ≈ 1e-4 steps in insets |
| 7 | 4×4 | E | T_net, T_L, T_C vs e, h_p = .06 | LIN; DATA | plunge after e ≈ h_p; T_C ≡ 0 in (d); transonic line at .06 |
| 8 | 4×4 | E | τ_L⁻¹, τ_a⁻¹, τ_e⁻¹/100 vs e, h_p = .06 | LIN; DATA | blue identical in all panels, peak 1.00 at e = .08; red peak at .07 |
| 9 | 2 stacked | F | a/a_0, e/e_0 vs t/T_0 (top axis t/yr), 3 discs | DER | red fastest a-decay; e/e_0 curves coincide |
| 10 | 1 + inset + cbar | F | (a, e) tracks coloured by t (inferno); inset t(a = .995) vs e_0 | DER | 6 tracks start black at a = 1; e_0 = .12 hooks right; inset minimum at e_0 ≈ .095 |
| 11 | 2 stacked | G | τ_L⁻¹ and τ_e⁻¹ vs e, (q,p,h,b) = (1,1.5,.05,.5) | black LIN; green/blue/red/orange dashed FIT (52–55) | black τ_e⁻¹ rises at low e and merges with PL00/CN06(Q = .2) for e ≥ .07 |
| A1 | heat map + slices | H | T_error over (m_max, \|l−m\|_max), e = .12 | LIN | yellow left/bottom edges, purple top-right; red/blue slices |
| B1 | 2 stacked | D | m = 1 F_J and T for q = 1.0, 1.005, 0.995 | LIN | red jump at OLR r = 1.587; no legend |
| C1 | 4×4 | E | as Fig 7, h_p = .08 | LIN; DATA | transonic line at .08 |
| C2 | 4×4 | E | as Fig 7, h_p = .10 | LIN; DATA | transonic line at .10; late zero crossings |
| C3 | 4×4 | E | as Fig 8, h_p = .08 | LIN; DATA | blue peak ≈ .73 near e ≈ .11 |
| C4 | 4×4 | E | as Fig 8, h_p = .10 | LIN; DATA | red and blue rise monotonically to e = .12; black turns down after ≈ .08–.10 |

**Per-figure VLM MATCH CHECKLIST**
- **All figures:**
  - [ ] usetex with CM-sans plain text and CM-serif math/ticks;
  - [ ] ticks out, left and bottom;
  - [ ] named k/r/b/g colours with fixed semantics (black net, red Lindblad,
    blue corotation);
  - [ ] tab:gray dashed reference lines.
- **Fig 1:**
  - [ ] top horizontal viridis bar ±2 with its label above;
  - [ ] (a)(b)(c) above the maps;
  - [ ] left/bottom-spine maps to ±2;
  - [ ] planet dot at (1, 0);
  - [ ] insets only in (a), (b).
- **Fig 2:** [ ] log–log, m from 2 to 150; legend `Numerical`/`Analytical`
  in sans.
- **Fig 3:** [ ] column titles with p/q/h values; ko- / rx / b-- encoding;
  half-clipped end markers.
- **Fig 4:** [ ] tab10 colours; `0.0T` … `0.75T` labels; coloured r_p dots;
  y label with a wide gap before `[F_{J,0}/a]`.
- **Fig 5:**
  - [ ] panel titles instead of y labels;
  - [ ] jet colourbar titled `e`;
  - [ ] legend text `T_LB` (sic);
  - [ ] italic panel letters.
- **Fig 6:** [ ] two insets with offset tick text `+1.546`; three light-red
  dashed lines.
- **Figs 7/8/C1–C4:**
  - [ ] figure-top 3-column legend;
  - [ ] p titles on top, q labels rotated at the right;
  - [ ] letters (a)–(p): lower-left for torques, upper-left for timescales;
  - [ ] transonic and zero dashed lines.
- **Fig 9:** [ ] hspace 0; top t/yr axis; legend in (b).
- **Fig 10:** [ ] inferno scatter tracks starting black; top colourbar 0–350;
  inset upper left.
- **Fig 11:** [ ] spelled-out legend names; (a)/(b) at upper right; CN06
  Q_e colours red/blue/orange.
- **Fig A1:** [ ] colourbar label below the bar; dashed red/blue slice
  markers on the map edges; twin top axis on the slice panel.
- **Fig B1:** [ ] no legend; OLR dashed line; hspace 0.

---

## 10. Defects and caveats found while inspecting

1. **Fig 5(c) legend reads `$T_\mathrm{LB}$` [SOLID]**, while the caption and
   the text call it T_L. Reproduce the glyph as published, or record the
   deviation.
2. **Fig 4 uses `T` for the orbital period** in the legend (`0.25T`), while
   the y label uses T for torque (notation clash in the published figure).
3. **Fig 10 colour scale:** the colourbar stops at 350, but the inset gives
   t(a = 0.995) ≈ 365 T_{a0} for e_0 = 0.02. The last part of that track is
   colour-saturated [PRELIMINARY: inset read by eye].
4. **Fig 8 / C3 / C4 blue series:** the near-identity across the 16 panels
   is REAL. It was checked against `anc/eccentric/tau_e_inv.out` (index
   order [q][p][h][e], 4×4×3×13 = 624 rows). It must not be "fixed".
5. **Fig B1 has no legend**, and its caption calls q = 1.005 black and
   q = 0.995 green "almost overlapping", which matches the render.
6. **Fig 2 analytic curve:** it needs the PL00 shifted Lindblad radii, which
   FR24 does not give (brief).
