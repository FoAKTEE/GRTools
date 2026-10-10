# FIGURE STYLE MANUAL — arXiv:2506.19917 (Fairbairn & Dittmann 2025, "Pushing the limits of eccentricity in planet-disc interactions", MNRAS) — FD25

VLM-inspected from the source figures in
`/data/haiyangw/claude/emri-disk/ref-paper/arXiv-2506.19917/src/Figures/`:
- 14 matplotlib PDFs, rendered at 200 dpi, with crops at 400–500 dpi;
- 3 PNG rasters, viewed downsampled and at full resolution.

They were cross-checked against the compiled PDF `2506.19917.pdf` (pages
6–22, 150 dpi).

| Fig | file |
|---|---|
| 1 | `sigma_map_lr.png` |
| 2 | `sigma_grid_exps_lr.png` |
| 3 | `sigma_diff.png` |
| 4 | `trq_grid_p1.5q0.0iso.pdf` |
| 5 | `trq_grid_p0.5q0.0iso.pdf` |
| 6 | `wiggle_modes.pdf` |
| 7 | `torque.pdf` |
| 8 | `taua.pdf` |
| 9 | `taue.pdf` |
| 10 | `amf_deposit_18k.pdf` |
| 11 | `theory_comparison.pdf` |
| 12 | `corotation_e0.01.pdf` |
| 13 | `corotation_e0.3.pdf` |
| A1 | `dTdr_convergence.pdf` |
| A2 | `net_T_convergence.pdf` |
| A3 | `dTdr_convergence_sim.pdf` |
| A4 | `dTdr_err.pdf` |

Captions are from `src/Sections/4_results.tex`, `5_discussion.tex` and
`app_modal_convergence.tex`. Mission brief:
`emri-disk/results/emri-eccentric/paper_EMRI-ECC/evidence/acquire-2506.19917/BRIEF.md`.
Second author A. J. Dittmann is also first author of arXiv:2608.17003.

Tags: [SOLID] = metadata / pdffonts / measured pixels. [PRELIMINARY] = by eye.

READ THE RELEVANT SECTION BEFORE GENERATING ANY REPRODUCTION FIGURE; then
VLM-compare per the `vlm-figure-repro` skill.

---

## 0. Shared style conventions (whole paper)

**Toolchain [SOLID]:**
- All figures are **matplotlib 3.8.0**, with **`text.usetex=True`,
  `font.family='serif'`**. Text and math are Computer Modern serif (CMR10,
  CMMI, CMSY).
- A few strings fall back to **CM Sans** (CMSS8):
  - the `Disco−Linear` box of Fig 3;
  - the legend entries `PL00` and `MTI11` of Fig 11;
  - the `max` subscript in the A1/A2 colourbar label `$|m-l|_\mathsf{max}$`.
- The column headers of Fig 1 set `Disco` in typewriter (`\texttt`).
- The PNG maps are matplotlib rasters:
  - Fig 1: 3401 × 3617 px @ 500 dpi (6.8 × 7.2 in);
  - Fig 2: 2735 × 1326 px @ 400 dpi;
  - Fig 3: 3376 × 2523 px @ 1000 dpi.

**Two tick styles [SOLID for direction and sides; lengths PRELIMINARY]:**
- **Maps (Figs 1–3):**
  - major ticks only, on ALL four sides;
  - drawn **crossing the spine** (`direction='inout'`, length ≈ 8 pt), with
    no minor ticks;
  - labels on the left and bottom only;
  - panels inside a block share their edges (`wspace = hspace = 0`).
- **Line plots (all others):**
  - long **outward** major ticks (≈ 5–6 pt) on the left and bottom only;
  - minor ticks only on log axes;
  - full black box (lw ≈ 1.2), no grid.

**Palette [SOLID, measured hex]:** the SEMANTIC pair used throughout the
results is **simulation (Disco) = navy #202d85** and **linear theory = amber
#ffab00**:
- linear theory: dotted in the profile plots, `'s-'` small squares in the
  e-scans;
- simulation: solid in the profile plots, large `'*'` stars in the e-scans.

Other colours:
- reference prescriptions: **PL00 = green #008000** (dashed); **MTI11 =
  rust #a93400** (dash-dot);
- radial markers in Fig 10: rust #a93400 = planet peri/apocentre, green
  #008000 = a_p ± l_sh;
- time slices in Figs 12–13 (lines at alpha 0.8): t = 2 T_p rust #a93400,
  25 T_p green #008000, 100 T_p navy #202d85, 190 T_p violet #cb6ded;
- appendix: A3 red #ff0000 dashed; A4 linear = matplotlib `'orange'`
  #ffa500 dotted (NOT the #ffab00 amber of the main text); A1/A2 discrete
  viridis;
- Fig 6: black, red #ff0000, `'orange'` #ffa500, gray; nodal guides blue at
  alpha 0.4 (#9999ff) dotted.

**Reference lines [SOLID]:**
- zero lines are dashed grey: #808080, lw ≈ 2.5, in Figs 7–9 and A2–A4;
- #a5a5a5 (gray at alpha 0.7) in Figs 4, 5, 10;
- planet extent lines (a(1 ± e)) are black dashed in Figs 4–5.

**Label boxes [SOLID]:**
- **Line plots:** rounded boxes, face `gainsboro` #dcdcdc, no edge. They
  carry the disc stamp `$p = 1.5\ q = 1.0$` (no comma), `$e = 0.01$`, and the
  panel letters `$(a)$`, `$(b)$` (italic).
- **Maps:** rounded boxes, white face, black edge (lw ≈ 1.5).

Legends are framed (rounded, light-grey edge).

**Colormap for δΣ maps [SOLID hex; identification PRELIMINARY]:**
- The ramp is a diverging **dark maroon → off-white → dark navy**.
- Measured stops (t = 0 … 1): #3c0912, #7e0e29 (0.1), #b33826 (0.2),
  #dcb2a4 (0.4), **#f1eceb (0.5)**, #aac2cb (0.6), #1670bc (0.8), #293f9a
  (0.9), #181c43 (1.0).
- These end and centre colours are those of **cmocean `balance`, reversed**
  (`cmocean.cm.balance_r`): negative δΣ is red, positive is blue.
- cmocean is not installed on this host. Verify the identification there
  before relying on it, or embed the measured stops as a
  `LinearSegmentedColormap`.

**Colour limits (measured from tick pixel positions) [SOLID]:**
- Fig 1: **±2.80**, ticks −2 … 2;
- Figs 2 and 3: **±5.45**, ticks −4 … 4 step 2.

These are NOT round numbers; reproduce them literally.

**Figure sizes [SOLID, PDF page boxes; inches]:**

| figure(s) | size (in) |
|---|---|
| 4, 5 | 7.3 × 8.6 (two-column) |
| 6 | 3.48 × 3.30 |
| 7, 8, 9 | 3.5–3.56 × 4.92 |
| 10 | 6.95 × 5.30 (two-column) |
| 11 | 3.63 × 3.42 |
| 12, 13 | 3.9–4.0 × 3.96–4.0 |
| A1 | 4.97 × 8.27 |
| A2 | 4.07 × 5.55 |
| A3 | 3.68 × 4.08 |
| A4 | 3.49 × 4.08 |

Single-column figures are included at `\linewidth`, so fonts are ≈ 1:1 on the
page (ticks ≈ 8 pt, labels ≈ 9–10 pt [PRELIMINARY]).

Reference rcParams block:
```python
import matplotlib as mpl, numpy as np
mpl.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.size': 9,
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.major.size': 5.5, 'ytick.major.size': 5.5,
    'axes.linewidth': 1.2, 'axes.grid': False, 'legend.frameon': True,
})
SIM, LIN = '#202d85', '#ffab00'          # simulation navy, linear amber
PL00, MTI11 = '#008000', '#a93400'
BOX = dict(boxstyle='round', fc='gainsboro', ec='none')
BALANCE_R_STOPS = ['#3c0912', '#7e0e29', '#b33826', '#ca765b', '#dcb2a4', '#f1eceb',
                   '#aac2cb', '#5c9dbb', '#1670bc', '#293f9a', '#181c43']   # t = 0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1
CMAP = mpl.colors.LinearSegmentedColormap.from_list('balance_r_meas',
         list(zip([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0], BALANCE_R_STOPS)))
# for maps: ax.tick_params(direction='inout', top=True, right=True, length=8)
```

---

## 1. Plot type A — δΣ maps: sim vs linear, peri vs apo (Figs 1, 2, 3)

Quantity: `$(M_\mathrm{p}/M_\mathrm{th})^{-1}\delta\Sigma/\Sigma$`.
- **Frame:** star-centred, planet on the +x axis.
- **Field of view:** a square [−2.5, 2.5]² (Fig 1) or [−2, 2]² (Figs 2–3).
  The disc fills the whole square, with no white corners. There is a small
  white central hole (r < 0.05).
- **No planet marker** is drawn.

**Fig 1 layout [SOLID]:**
- 4 rows (e = 0.01, 0.06, 0.3, 0.5) × 4 columns in TWO blocks of two:
  - [Disco | Linear] at **pericentre**;
  - a small gap;
  - [Disco | Linear] at **apocentre**.
- Within a block: wspace = hspace = 0.
- **Headers:** boxed `Disco` (typewriter) / `Linear` above each column, and
  boxed `Pericentre` / `Apocentre` centred under each block.
- **Row stamps:** boxed `$e = 0.01$` etc. at the lower left of the first
  column.
- **Colourbar:** ONE horizontal bar across the top, ticks on its top side,
  label above `$(M_\mathrm{p}/M_\mathrm{th})^{-1}\delta\Sigma/\Sigma$`,
  ±2.80.
- **Axes:** `$x$` / `$y$` labels; ticks −2 … 2 every 1, labelled only on the
  outer left column / bottom row.

**Fig 2:**
- 1 × 3 square panels, separated by small gaps.
- Boxed stamps top-left: `$p = 1.5\ q = 0.0$`, `$p = 1.5\ q = 1.0$`,
  `$p = 0.5\ q = 0.0$`.
- Top horizontal colourbar ±5.45 with the same label.
- y label on the first panel only.

**Fig 3:**
- A single square map with a vertical colourbar at the right, ±5.45, label
  rotated, reading bottom → top.
- Boxed sans text top-left: `Disco−Linear`.

**Dominant visual STRUCTURE:**
- **Background:** near-white. Arms are thin blue (+) lines with a faint red
  (−) companion.
- **e = 0.01 (Fig 1 row 1):** ONE trailing spiral from a kinked wake at
  (1, 0). The outer arm wraps clockwise to r ≈ 1.8, then a second loop
  reaches r ≈ 2.4. The inner arm winds into a ring at r ≈ 0.45, with tight
  noisy rings inside. Disco and Linear look near-identical, and peri vs apo
  look the same.
- **e = 0.06:** same topology. The wake near (0.95–1.05, 0) is doubled (a
  second short crest just inside). The outer arms have stronger red/blue
  contrast in Linear.
- **e = 0.3:** a crescent-shaped double-crest wake hanging below the planet
  position: from (0.75, 0) to (0.4, −1.3) at pericentre, and near (1.3, 0)
  pointing up-left at apocentre. Many concentric arcs fill the disc. Linear
  shows more fine ring structure (dense red/navy speckle) than Disco.
- **e = 0.5:** the wake is a faint crescent near (0.5, −1.5) at pericentre
  and a thin arc ending at (1.5, 0) at apocentre. The disc is filled with
  many low-contrast concentric arcs.
- **Fig 2:** three nearly identical single-arm patterns at e = 0.06
  pericentre. The (1.5, 1.0) disc has visibly more outer-arm contrast near
  the frame edge.
- **Fig 3:** almost everything cancels to off-white. Residuals are thin
  red/blue lines along the outer arm at r ≈ 1.6–2.0, a thin ring at
  r ≈ 0.5, and the inner noisy rings (r < 0.45). Only a faint residual is
  left at the planet's wake.

**Recipe (Fig 1):**
```python
fig = plt.figure(figsize=(6.8, 7.2))
outer = fig.add_gridspec(1, 2, wspace=0.05)
for b, phase in enumerate(['peri', 'apo']):
    g = outer[b].subgridspec(4, 2, wspace=0, hspace=0)
    for i, e in enumerate([0.01, 0.06, 0.3, 0.5]):
        for j, src in enumerate(['Disco', 'Linear']):
            ax = fig.add_subplot(g[i, j]); ax.set_aspect('equal')
            im = ax.pcolormesh(X, Y, dS[(phase, e, src)], cmap=CMAP, vmin=-2.8, vmax=2.8,
                               shading='auto', rasterized=True)
            ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5)
            ax.tick_params(direction='inout', top=True, right=True, length=8)
# headers/stamps: ax.text(..., bbox=dict(boxstyle='round', fc='white', ec='k', lw=1.5))
# fig.colorbar(im, ax=..., orientation='horizontal', location='top', ticks=[-2,-1,0,1,2])
```
Save as PNG (500 dpi; the original is a raster). Linear columns come from
GRTools; Disco columns are SIM, reference only (§9).

---

## 2. Plot type B — dT/dr profiles, linear vs sim, with wiggle zooms (Figs 4, 5)

**Layout [SOLID]:**
- 5 rows (e = 0.01, 0.06, 0.12, 0.30, 0.50).
- Each row: a WIDE left panel (≈ 60% of the width) plus two small stacked
  zoom panels at the right. Total size 7.3 × 8.6 in.
- **Left panel:**
  - `$dT/dr\ [q_\mathrm{r}^2 h_\mathrm{p}^{-3} a_\mathrm{p}^{-1}]$` vs
    `$r/a_\mathrm{p}$`, 0.1–4.0 (ticks 0.5 … 4.0 step 0.5);
  - x label in every row (it sits close under the ticks);
  - black dashed vertical lines at r = a(1 − e) and a(1 + e).
- **Zoom panels:**
  - inner disc: r/a_p 0.1–0.65, y ≈ −0.15 … 0.25 (ticks 0.00, 0.25);
  - outer disc: r/a_p 1.35–4.0, y ≈ −0.13 … 0.11 (ticks −0.1, 0.0, 0.1);
  - each has a grey #a5a5a5 thick dashed zero line.
- **Series:** simulation = navy #202d85 solid, lw ≈ 2; linear = amber
  #ffab00 DOTTED, lw ≈ 2, drawn on top.
- **Legend:** framed, in the row-1 left panel, lower right:
  `Linear theory`, `Simulation`.
- **Row stamp:** gainsboro box at the lower left: `$e = 0.01$` … `$e = 0.50$`
  (two decimals).

**Dominant visual STRUCTURE (Fig 4, p = 1.5, q = 0):**
- **e = 0.01:** the classic one-sided pair, −6 at 0.95 and +7.8 at 1.04.
- **e = 0.06:** +5.2 at 0.97, a sharp dip, then +2.4 at 1.08, with −2.2 at
  0.91.
- **e = 0.12:** REVERSED, with a +15.5 peak at 0.89 and −20 at 1.10.
- **e = 0.30:** +3.3 spike at 0.70, small ripples, then −7.2 at 1.30.
- **e = 0.50:** many spikes (±1) between 0.5 and 1.3, then −2.9 at 1.50.
- **Extrema follow the orbit:** they move with the dashed peri/apo lines.
- **Linear vs sim:** they coincide everywhere except in the e = 0.50 zooms
  (inner r ≈ 0.5–0.65 and outer r ≈ 1.35–2.0), where the dotted amber
  departs from navy.
- **Zooms:** a characteristic wiggle train, ±0.05 outer bumps at r ≈ 1.85,
  2.4, 2.9, 3.35, 3.8, and an inner ladder growing to +0.24 at r ≈ 0.51.

**Fig 5 (p = 0.5, q = 0):** the same layout, with these differences:
- e = 0.01: −6.5 / +7.4;
- e = 0.06: two narrow spikes, +3.3 at 0.97 and +2.4 at 1.08, around a deep
  −2.4 trough at 1.03;
- e = 0.12: +14 / −22;
- e = 0.30: +2.3 / −9;
- e = 0.50: −4.3 at 1.5.

---

## 3. Plot type C — mode spectrum amplitude and phase (Fig 6)

- **Layout:** two stacked panels, hspace = 0, shared x `$m$` 0–80 (ticks 0,
  20, …, 80).
- **(a)** `$|\delta\Sigma_{mm}|/\delta\Sigma_\mathrm{lin}$`, 0–0.12 (ticks
  0.00, 0.05, 0.10).
- **(b)** `$\arg[\delta\Sigma_{mm}\exp(im\phi_\mathrm{lin})]/\pi$`, −1 … 1.
- **Series** (framed legend in (a), upper right; label format
  `$e = 0.0\ \ r/r_\mathrm{p} = 3.0$`):

  | (e, r/r_p) | line |
  |---|---|
  | (0.0, 3.0) | black solid |
  | (0.0, 4.0) | red #ff0000 dashed |
  | (0.06, 3.0) | orange #ffa500 dashed |
  | (0.06, 4.0) | gray solid |

- **Guides:** blue alpha-0.4 dotted vertical lines at m ≈ 21.5, 50.5, 80
  (the nodes).
- **Letters:** gainsboro-boxed italic `$(a)$`, `$(b)$` at the lower right.
- **Structure:**
  - **(a):** the e = 0 curves (black, red) coincide. They rise from 0.035 at
    m = 1 to a peak of 0.105 at m ≈ 8, then decay to ≈ 0.01 at 80.
  - The e = 0.06 curves (orange, gray) coincide: a peak of 0.095 at m ≈ 6,
    then **zeros at m ≈ 21.5 and 50.5** with lobes between (≈ 0.03 at
    m ≈ 33; 0.012 at m ≈ 63).
  - **(b):** the e = 0 phases rise smoothly from −1 to ≈ 0, then drift; red
    is above black after m ≈ 20.
  - The e = 0.06 phases jump by ≈ π at each node (m ≈ 21.5 and 50.5).
    - Orange drops to ≈ −0.95 at 21.5, wraps from −1 to +1 at m ≈ 33, then
      returns onto the black curve at 50.5.
    - Gray drops to ≈ −0.8 at 21.5, then jumps onto the red curve at 50.5
      and drops again at 80.

---

## 4. Plot type D — integrated quantities vs e, three discs (Figs 7, 8, 9)

**Layout [SOLID]:**
- Three stacked panels, hspace = 0, shared x `$e$` 0–0.7 (ticks 0.0 … 0.7).
- Rows: `$p = 1.5\ q = 1.0$`, `$p = 1.5\ q = 0.0$`, `$p = 0.5\ q = 0.0$`
  (gainsboro box: lower right in Fig 7, upper right in Figs 8 and 9).
- **y labels:**
  - Fig 7: `$T\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}]$`;
  - Fig 8: `$\tau_a^{-1}/\tau_0^{-1}$`;
  - Fig 9: `$\tau_e^{-1}/\tau_0^{-1}$`.
- **Series:**
  - linear = amber `'s-'` (ms ≈ 3.5, lw ≈ 1.5). Dense in e: 0, 0.01, …,
    0.06, then ≈ 0.005 steps to 0.1, then 0.11, 0.12, 0.15, 0.18, 0.21,
    0.24, 0.27, 0.30, 0.4, 0.5, 0.6, 0.7;
  - simulation = navy `'*'` stars, ms ≈ 13, no line. About 11 points:
    0.01, 0.06, ≈ 0.07–0.08, 0.12, 0.18, 0.24, 0.30, 0.4, 0.5, 0.6, 0.7
    [PRELIMINARY].
- Grey #808080 thick dashed zero line in each panel.
- Framed legend in the top panel (`Linear theory`, `Simulations`): upper
  right in Fig 7, right-middle in Fig 8, upper middle in Fig 9.
- **Fig 9 insets:** in each panel, an inset (upper right, about 45% × 45%
  of the panel) on **log–log** axes, x 0.1–1.0 and y ≈ 1e-1 … 1e2. It shows
  the same two series plus a thick grey solid guide line labelled
  `$\propto e^{-3}$`.

**Dominant visual STRUCTURE:**
- **Fig 7:**
  - T is positive at low e (top 0.25 → peak 0.45 at e ≈ 0.06; middle
    0.18 → 0.35; bottom flat 0.12).
  - It plunges through 0 at e ≈ 0.08 to a minimum at e ≈ 0.12 (−0.6 / −0.55
    / −0.78), then recovers slowly toward −0.1 / −0.08 / −0.2 at e = 0.7.
  - The stars sit BELOW the amber line for e ≥ 0.4.
- **Fig 8:**
  - a sharp peak at e ≈ 0.07–0.08 (≈ 1.75 / 1.65 / 0.98), then a fall to a
    low plateau;
  - for **p = 0.5 the rate goes negative** from e ≈ 0.12 (minimum ≈ −0.15 at
    e ≈ 0.15–0.18) and returns to ≈ 0 at 0.7;
  - a small local **bump at e = 0.24** (amber 0.24 and 0.27 nearly equal;
    star at 0.24 above its neighbours).
- **Fig 9:**
  - a peak of ≈ 100 at e ≈ 0.08 in all three discs, then a decay to ≈ 0 by
    e ≈ 0.3;
  - insets: amber and stars fall roughly parallel to the e⁻³ guide, and the
    stars are slightly above amber at e ≥ 0.5.

---

## 5. Plot type E — AMF and deposition profiles (Fig 10)

- **Layout:** 4 rows (e = 0.01, 0.06, 0.12, 0.30) × 2 columns.
- **Left column:** `$F_J\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}]$` vs
  `$r/a_\mathrm{p}$` 0.4–2.0. Linear = amber dotted; simulation = navy solid.
  Framed legend in row 1, lower right.
- **Right column:** `$dF_\mathrm{dep}/dr\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}a_\mathrm{p}^{-1}]$`
  vs `$r/a_\mathrm{p}$` 0.4–3.0. SIMULATION ONLY (navy).
- **Markers in all panels:**
  - rust #a93400 dashed verticals at a(1 ± e);
  - green #008000 dashed verticals at a_p ± l_sh (≈ 0.555 and 1.445);
  - #a5a5a5 dashed zero line.
- **Stamps:** gainsboro `$e = 0.01$` … in the left panels, lower left.
- **Structure:**
  - **Left:** linear F_J is flat outside the orbit (≈ 0.58 / 0.77 for
    e = 0.01; 1.22 / 1.56; 2.55 / 2.0; 1.03 / 0.75). The simulation follows
    linear near the orbit but **decays away from the orbit on both sides**,
    falling below linear. At e = 0.01 the split is at ≈ a_p ± l_sh (the green
    lines). At e ≥ 0.06 it starts closer to the orbit (r ≈ 0.8 and 1.2).
  - e = 0.12 shows a central hump (3.4 at r ≈ 1); e = 0.30 shows a flat
    top 1.12 between the peri/apo lines.
  - **Right:** deposition is negative inside and positive outside, with
    extremes ≈ −0.45 / +0.39 (0.01), −1.8 / +1.75 (0.06), −3.8 / +2.0
    (0.12), −0.9 / +0.33 (0.30). It is ≈ 0 between the peri/apo lines.

---

## 6. Plot type F — log–log comparison with prescriptions (Fig 11)

- **Layout:** two stacked panels, hspace = 0, log–log, shared x `$e$` ≈
  0.05–1.0 (labelled 10⁻¹, 10⁰).
- **Top:** `$\tau_e^{-1}/\tau_0^{-1}$` 1e-2 … 1e3. **Bottom:**
  `$\tau_a^{-1}/\tau_0^{-1}$` ≈ 2e-2 … 1e1.
- **Series:**
  - linear = amber `'s-'`; simulation = navy stars;
  - `PL00` = green dashed; `MTI11` = rust dash-dot;
  - a grey solid `$\propto e^{-3}$` guide (top panel only).
- **Legend and stamp:** framed legend in the top panel, lower left (`PL00`
  and `MTI11` in CM Sans). Gainsboro stamp `$p = 1.5\ q = 1.0$` at the top
  right.
- **Structure:**
  - **Top:** all curves fall at high e. MTI11 starts highest (≈ 500 at
    e = 0.065) and **plunges near e → 1**. PL00 lies above linear by ≈ ×2
    for e ≳ 0.2. Amber and the stars hug each other.
  - **Bottom:**
    - amber and stars peak at ≈ 1.7 at e ≈ 0.07–0.08, fall to ≈ 0.17 by
      e = 0.15, then show a plateau and bump to e ≈ 0.25, then a decline
      to 2.4e-2 (amber) at 0.7;
    - the stars stay higher (6e-2 at 0.7);
    - PL00 dips to 0.25 at e ≈ 0.11, plateaus at ≈ 0.5, and **diverges near
      e → 1**;
    - MTI11 declines from 2.7 to 0.33, with an upturn at e → 1.

---

## 7. Plot type G — coorbital AMF time evolution (Figs 12, 13)

- **Layout:** two panels (NOT shared, with a gap).
- **(a)** `$F_J\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}]$` vs `$r/a_\mathrm{p}$`
  0.97–1.03:
  - four simulation time slices, alpha 0.8, lw ≈ 2: `$t = 2T_\mathrm{p}$`
    rust, `25T_p` green, `100T_p` navy, `190T_p` violet;
  - plus `Linear theory`, amber dotted (Fig 12 only);
  - framed legend: lower left (Fig 12), upper left (Fig 13).
- **(b)** `$T_{r\sim a_\mathrm{p}}\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}]$` vs
  `$t/T_\mathrm{p}$` 0–200, black line, with dashed verticals at
  t = 2, 25, 100, 190 in the slice colours.
- **Letters:** gainsboro `$(a)$`, `$(b)$` at the lower right.
- **Structure, Fig 12 (e = 0.01):**
  - amber has a **vertical step at r = 1.000** (0.035 → 0.095);
  - the sim slices are smooth S-shaped curves through ≈ 0.08 at r = 1;
  - 2 T_p (rust) is a shallow U (minimum 0.07);
  - 25 and 100 T_p dip deepest (minimum ≈ 0.035 at r ≈ 0.996);
  - 190 T_p is intermediate;
  - (b): T drops from 0.09 to 0.02 within t ≈ 3, declines to ≈ 0.003 by
    t ≈ 100, then rises to 0.047 at 200.
- **Structure, Fig 13 (e = 0.30):**
  - F_J 1.00–1.10, smooth rising lines with only a weak S near r = 1;
  - order green (25) > navy (100) > violet (190) ≫ rust (2, ≈ 1.00–1.035);
  - no linear curve (off-scale, as the caption says);
  - (b): flat at ≈ −0.008, range ±0.05.

---

## 8. Plot type H — convergence diagnostics (Figs A1–A4)

**A1 and A2: modal-truncation families.**
- **Colourbar:** discrete **viridis with 9 bins**, `viridis(k/40)` for
  k = 0, 5, …, 40, horizontal on top. Ticks at the bin centres 0 … 40, label
  above `$|m-l|_\mathsf{max}$`.
- **Curves:** only k = 5 … 40 are drawn; the k = 0 bin colour #440154 is
  absent from the panels [SOLID, pixel count].
- **A1:**
  - 4 stacked panels at r/a_p 0.4–1.6. The first three show dT/dr for the
    three discs (gainsboro stamps lower left, y ticks −2, 0).
  - The 4th is a **symlog** y axis `$\Delta dT/dr$` (±10⁰, linthresh ≈ 1e-3)
    with red #ff0000 = `$|m-l|_\mathsf{max} = 80$` and gray =
    `$m_\mathsf{max} = 220$`, legend lower right.
  - Structure: the curves converge toward the yellow 40 curve. The big
    trough at r ≈ 1.49 deepens from −0.3 (k = 5) to −3 (top) / −2.7 / −4
    (k = 40).
  - Red residual: ±0.1, oscillating densely between 0.8 and 1.3, ≈ +0.1 at
    r = 1.3–1.5. Gray residual: ≲ 1e-3 everywhere.
- **A2:**
  - 3 stacked panels `$T[q_\mathrm{r}^2\,h_\mathrm{p}^{-3}]$` vs `$e$`
    0–0.7, `'s-'` per k, grey dashed zero, gainsboro stamps lower right.
  - The bottom panel adds red × (`$|m-l|_\mathsf{max} = 80$`) and gray +
    (`$m_\mathsf{max} = 220$`) at e = 0.5.
  - Structure: the k = 5 curve is far shallower (minimum −0.25 / −0.25 /
    −0.33); the curves converge for k ≥ 25. The red ×/gray + at e = 0.5
    land on the k = 40 curve (≈ −0.21).

**A3 and A4: simulation resolution and linear-vs-sim residuals.**
- **Layout:** 3 stacked panels, small gaps.
  - (a) `$dT/dr$` with two curves;
  - (b) `$\Delta dT/dr$`;
  - (c) cumulative `$\Delta T\,[q_\mathrm{r}^2 h_\mathrm{p}^{-3}]$`.
- **Axes:** A3 r/a_p 0.1–2.0; A4 0.1–4.0.
- **Series:**
  - A3: black solid `Fiducial resolution` + red dashed `Double resolution`;
  - A4: black `Fiducial resolution` + `'orange'` dotted `Linear theory`;
  - framed legends in (a).
- Grey dashed zero lines; gainsboro letters at the lower right.
- **Structure:**
  - **A3:** curves are indistinguishable (trough −4.3 at 1.5); (b) spikes
    ±0.13 at r ≈ 0.6–1.2; (c) |ΔT| ≲ 7.5e-4, ending ≈ −7.5e-4 at 2.
  - **A4:** (b) spikes ±0.25 at r ≈ 0.5–1.5; (c) ΔT ≈ 0 to r ≈ 0.9, then a
    monotone fall to ≈ −0.06 by r ≈ 3.5 (the nonlinear deposition gap).

---

## 9. Per-figure table

Provenance:
- **LIN** = linear theory, the target of reproduction with GRTools.
- **SIM** = Disco hydro. Reference only: do not re-run by default. Plot it
  only if digitised or supplied, and label it as such.
- **FIT** = published prescription evaluated from formulas (PL00 = FR24 eqs
  52–53; MTI11 needs arXiv:1106.0417).
- **GUIDE** = analytic guide line.
- **DATA** = linear curves published in `ref-code/callumf8-linear_ecc_planet`
  (e ∈ {.01, .06, .12, .3}, no e = .5). These are golden comparisons only.

| Fig | panels | type | what is plotted | provenance per series | structure to check |
|---|---|---|---|---|---|
| 1 | 4×4 (2 blocks) + top cbar | A | δΣ/Σ, (p,q) = (1.5,0), e = .01/.06/.3/.5, peri & apo | cols 1, 3 SIM; cols 2, 4 LIN | ±2.80 balance_r; one arm at low e → crescent wakes at e ≥ .3; boxed headers |
| 2 | 1×3 + top cbar | A | Disco δΣ/Σ, e = .06, pericentre, 3 discs | SIM | ±5.45; three near-identical single-arm maps |
| 3 | 1 + right cbar | A | Disco − Linear δΣ, (1.5,0), e = .06 | SIM − LIN | near-total cancellation; outer-arm and r ≈ 0.5 residual rings |
| 4 | 5 rows × (1 wide + 2 zooms) | B | dT/dr, (1.5,0), e = .01, .06, .12, .3, .5 | navy SIM; amber LIN (DATA for e ≤ .3) | sign reversal by e = .12; extrema at a(1 ± e); e = .5 zoom mismatch |
| 5 | as 4 | B | as 4 for (0.5,0) | as 4 | as 4, with the values of §2 |
| 6 | 2 stacked | C | \|δΣ_mm\|/δΣ_lin and phase vs m, e = 0, .06 at r/r_p = 3, 4 | LIN (needs C+24 normalisation) | nodes at m ≈ 21.5, 50.5 with π phase jumps for e = .06 |
| 7 | 3 stacked | D | net T(e), 3 discs | amber LIN; navy stars SIM | positive → minimum at e ≈ .12 → slow recovery; stars below amber at e ≥ .4 |
| 8 | 3 stacked | D | τ_a⁻¹(e) | LIN; SIM | peak at .07–.08; p = .5 negative for e > .12; bump at .24 |
| 9 | 3 stacked + 3 log insets | D | τ_e⁻¹(e) | LIN; SIM; GUIDE e⁻³ | peak ≈ 100 at e ≈ .08; inset slopes parallel to e⁻³ |
| 10 | 4×2 | E | F_J (lin vs sim) and dF_dep/dr (sim), e = .01, .06, .12, .3 | amber LIN (DATA); navy SIM | sim F_J decays outside ≈ a ± l_sh; deposition sign −/+ |
| 11 | 2 stacked, log–log | F | τ_e⁻¹, τ_a⁻¹ vs e for (1.5,1.0) | LIN; SIM; FIT PL00 (green), FIT MTI11 (rust); GUIDE | MTI11/PL00 behaviour at e → 1; lin/sim τ_a peak at .07–.08 |
| 12 | 2 | G | coorbital F_J at t = 2, 25, 100, 190 T_p; T_{r∼a}(t), e = .01 | coloured SIM; amber LIN | linear step at r = 1; T(t) drop–dip–rise |
| 13 | 2 | G | as 12, e = .30 | SIM only (linear off-scale) | slice order 25 > 100 > 190 ≫ 2; T ≈ −0.008 |
| A1 | 4 stacked + top cbar | H | dT/dr for \|m−l\|_max = 5 … 40, e = .5; symlog residuals | LIN | 9-bin viridis bar, 8 curves; trough at 1.49 deepening with k |
| A2 | 3 stacked + top cbar | H | T(e) for \|m−l\|_max = 5 … 40; extended points at e = .5 | LIN | convergence for k ≥ 25; ×/+ on the k = 40 curve |
| A3 | 3 stacked | H | Disco fiducial vs double-resolution dT/dr, residual, cumulative | SIM | \|ΔT\| ≲ 7.5e-4 |
| A4 | 3 stacked | H | Disco fiducial vs linear dT/dr, residual, cumulative | SIM; LIN (orange #ffa500) | ΔT → −0.06 beyond r ≈ 1.5 |

**Per-figure VLM MATCH CHECKLIST**
- **All figures:**
  - [ ] usetex CM serif;
  - [ ] navy #202d85 = sim and amber #ffab00 = linear, never swapped;
  - [ ] gainsboro rounded stamps (line plots) / white black-edged boxes
    (maps);
  - [ ] outward long ticks (line plots) / inout all-side ticks (maps).
- **Figs 1–3:**
  - [ ] measured balance_r colours with limits ±2.80 (Fig 1) and ±5.45
    (Figs 2–3);
  - [ ] no planet marker;
  - [ ] Fig 1 block gap and the `Disco`/`Linear`/`Pericentre`/`Apocentre`
    boxes.
- **Figs 4–5:**
  - [ ] wide + 2-zoom row geometry;
  - [ ] dotted amber over solid navy;
  - [ ] black dashed a(1 ± e) lines;
  - [ ] the gainsboro `e = 0.01` stamp format.
- **Fig 6:** [ ] hspace 0; four series with the stated styles; blue alpha
  dotted node lines; π-jumps.
- **Figs 7–9:**
  - [ ] 3 stacked panels with hspace 0;
  - [ ] large navy stars;
  - [ ] amber 's-' with dense low-e sampling;
  - [ ] Fig 9 log insets with the `∝ e^{-3}` guide.
- **Fig 10:** [ ] rust/green dashed verticals; the right column is
  simulation only.
- **Fig 11:** [ ] log–log; PL00 green dashed, MTI11 rust dash-dot; sans
  legend names.
- **Figs 12–13:** [ ] alpha-0.8 time-slice colours rust/green/navy/violet,
  matched by the (b) dashed verticals.
- **A1/A2:** [ ] 9-bin viridis bar with 8 curves; symlog residual panel
  (A1).
- **A3/A4:** [ ] (a)/(b)/(c) stacks; A4 uses `'orange'` #ffa500, not amber.

---

## 10. Defects and caveats found while inspecting

1. **Figs 4–5 line style [SOLID]:** the caption says linear theory is
   "*yellow dashed*", but the figure (and its legend handle) draws an amber
   **dotted** line. Follow the figure.
2. **Fig 3 sign [SOLID]:** the in-figure box reads `Disco−Linear` (simulation
   minus linear). The caption only says "comparing", and the mission BRIEF
   row 3 says "linear minus simulation". The brief has the order reversed.
3. **Colour limits differ between maps of the same quantity [SOLID]:** Fig 1
   is ±2.80 and Figs 2–3 are ±5.45. Do not share one norm.
4. **A1 colourbar has an unused bin [SOLID]:** there are 9 bins (0 … 40),
   but only the 5 … 40 curves are drawn. This is consistent with the
   caption; reproduce the 9-bin bar.
5. **τ_a bump location [PRELIMINARY]:** the figure shows the local bump at
   e = 0.24 (amber 0.24 ≈ 0.27 within ≈ 0.005). The text has ".27" (4_results
   l.139) and ".24" (5_discussion l.26), so the figure supports .24.
6. **Mixed amber shades:** #ffab00 in the main text vs `'orange'` #ffa500 in
   Figs 6 and A4. They are visually near-identical, but keep the distinction
   in a faithful reproduction.
7. **Colormap identification:** the measured stops match cmocean `balance`
   reversed at its ends and centre. Confirm with cmocean (not installed here)
   before replacing the measured stop list.
