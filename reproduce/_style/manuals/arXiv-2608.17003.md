# FIGURE STYLE MANUAL — arXiv:2608.17003 (Dittmann, Hegade K. R., Fairbairn, "Extreme-Mass-Ratio Inspirals in Gaseous Disks")

VLM-inspected from the reference figures at
`/data/haiyangw/claude/emri-disk/ref-paper/arXiv-2608.17003/src/figures/`
(12 matplotlib PDFs + 2 PNGs), rendered at 160 dpi
(`scratchpad/wave0/fig2608/`) and 300 dpi (`scratchpad/wave1/fig2608hi/`).
Captions quoted from `main.tex` (figure environments, lines 853–1012) and
`appendix.tex` (lines 68–160). Parameter tables from the wave-0 brief §4.
Colour hex codes below were MEASURED from the 300-dpi renders (dominant
saturated pixels), colormap identifications were verified by sampling the
matplotlib colormaps (plasma/viridis/Set2) — tagged [SOLID]; line widths,
font sizes and marker sizes are pixel estimates — tagged [PRELIMINARY].

READ THE RELEVANT SECTION OF THIS MANUAL **BEFORE** GENERATING ANY
REPRODUCTION FIGURE. After generating, VLM-compare against the reference
per the `vlm-figure-repro` skill (3-step recurring procedure).

Macros used in captions: `\aEMRI` = $a_\bullet$, `\pdisk` = $\mathfrak{p}$,
`\qdisk` = $\mathfrak{q}$, `\bsoft` = $b$, `\mninefive` = $m_{95}$.

---

## 0. Shared style conventions (whole paper)

All 12 PDF figures are **matplotlib 3.10.x** (pdfinfo Creator: v3.10.1 for
chi_rp, convergence, dT_dr_versus_ap, pole_errs, rp_h0, soft_differences;
v3.10.8 for the rest) [SOLID]. The two spiral figures are PNG rasters
(spiral_grid.png 3276×1949 px @500 dpi; spirals_two_rows.png 3226×2001 px)
[SOLID].

**Typography [SOLID]:** real LaTeX (`text.usetex=True`) — embedded fonts are
Computer Modern Type 1 (CMR10/7/8/6, CMMI10/7/8/6, CMSY10/7, plus **EUFM**
Euler-Fraktur for $\mathfrak{p},\mathfrak{q}$). All labels, tick labels,
legend entries and in-axes annotations are CM serif; math in CM italic.
Fraktur $\mathfrak{p},\mathfrak{q}$ appear in axis labels, legends and
panel annotations — requires `\usepackage{amssymb}`/`amsfonts` in
`text.latex.preamble` (`\mathfrak`). Roman words inside math use
`\mathrm{}` (e.g. $T_\mathrm{GR}/T_\mathrm{Newt}$, $T_\mathrm{L}$,
$T_{m,\mathrm{ref}}$, $r_\mathrm{ISCO}$, $m_\mathrm{max}$).
Font sizes [PRELIMINARY]: tick labels ≈ 9–10 pt, axis labels ≈ 11–12 pt,
legend ≈ 9–10 pt, panel annotations ≈ 10 pt (figures are rendered at
~3.1–3.3 in width and included at `0.99\linewidth` ≈ 3.4 in, so sizes are
≈ 1:1 on the page).

**Figure size [SOLID, from PDF page boxes]:** single-column PDFs are
≈ 218–235 pt × 163–186 pt (≈ 3.1–3.25 in × 2.25–2.6 in, i.e. `figsize≈(3.3,
2.4)` + `bbox_inches='tight'`); the two-panel stacked figures
(chi_rp, rp_h0) are ≈ 225–228 × 186 pt; `convergence.pdf` is portrait
230 × 259 pt (≈ 3.2 × 3.6 in). Spiral PNGs are two-column (`figure*`,
`0.99\linewidth`), aspect ≈ 1.65:1.

**Axes/ticks [SOLID]:** ticks point **inward** on **all four sides**
(top and right spines carry ticks, labels only bottom/left); **minor ticks
on** (both axes, also on log axes); black spines ≈ 0.8–1 pt; **no grid**
except Fig 11 (light-gray major grid, both axes) and Fig 8 (light-gray
vertical gridlines only). This is the SciencePlots-`science`-like rc set:
`xtick.direction=ytick.direction='in'`, `xtick.top=ytick.right=True`,
`xtick.minor.visible=ytick.minor.visible=True`. Minus signs are proper
LaTeX minus (−) in the PDFs.

**Palette [SOLID, measured]:** the paper uses THREE recurring series
palettes — reproduce them literally:

1. *"Pantone-like" categorical list* (Figs 3, 5, 9, 10, 11, 12, 14), in
   cycle order:
   `#0f4c81` (classic blue) → `#e2583e` (orange-red) → `#009473` (emerald
   green) → `#bf1932` (true red) → `#5f4b8b` (ultra-violet) → `#a47864`
   (mocha brown). Set as `axes.prop_cycle` for these figures. Blue is always
   series 1, orange-red series 2.
2. *plasma sample* for the $h_0$ family (Figs 7, 13):
   `plt.cm.plasma(np.linspace(0.1, 0.9, 5))` =
   `['#41049d','#8f0da4','#cc4778','#f2844b','#fcce25']` for
   $h_0 = 0.05, 0.075, 0.1, 0.125, 0.15$ (dark purple → yellow), with the
   marker sequence square `s`, star `*`, circle `o`, triangle-right `>`,
   hexagon `h` (hollow in Fig 7 and 13).
3. *viridis sample* for the $a_\bullet$ family (Fig 4):
   `plt.cm.viridis([0.8, 0.6, 0.4, 0.2])` =
   `['#7ad151','#22a884','#2a788e','#414487']` for
   $a_\bullet = 10, 30, 100, 300\,M$ with linestyles `-`, `--`, `-.`, `:`;
   the $a_\bullet = 1000M$ (Newtonian-limit) curve is thin solid **black**.
4. Fig 6 ($\chi$ family) uses its own triple: `#e70234` (red, $\chi=0.9$,
   square), `#580772` (purple, $\chi=0$, star), `#2584da` (blue,
   $\chi=-0.9$, circle) — hollow markers.
5. Fig 8 violins: matplotlib **Set2** first three colours `#66c2a5`
   (teal, $a_\bullet=15M$), `#fc8d62` (orange, 30M), `#8da0cb` (blue-grey,
   120M).
6. Spiral maps: diverging **`RdBu`** (NOT reversed: negative $\delta U$ red,
   positive blue), symmetric limits, `extend='both'`.

**Legends [SOLID]:** two styles coexist — (a) matplotlib default framed
legend (`fancybox=True`, rounded corners, light-grey edge `0.8`, white face,
`framealpha≈0.8`): Figs 4, 5, 6, 7, 11, 12, 13; (b) `frameon=False`: Figs 3,
8, 9, 10. Legend handles for marker-only series are
single markers (`numpoints=1`/scatter handles); for Fig 8 the handles are
rounded filled colour patches with no edge. Multi-column legends (`ncol=2`
or 3) are used when 5–6 entries must fit: Figs 4, 10, 12 (ncol=3), 13.

**Panel annotations [SOLID]:** parameter stamps are plain text placed
inside the axes at the top-left (`ax.text(0.03, 0.95, ..., transform=
ax.transAxes, va='top')`), CM roman with inline math, e.g.
`$\mathfrak{p}=1.5,\ \mathfrak{q}=1.0$` (Figs 6, 7),
`$(\mathfrak{p},\mathfrak{q},h_0,\chi)=(0,0,0.15,-0.9)$` (Figs 1, 2),
`distribution over $h_0,\chi$` (Fig 8). No bounding box, except in the
spiral PNGs where a white box hides the field under the text.

**Lines/markers [PRELIMINARY]:** line plots use lw ≈ 1.2 pt (Fig 4) to
≈ 2 pt (Fig 3, the blue solid); scatter/marker plots are **marker-only (no
connecting lines)**; filled markers ≈ 7–9 pt diameter in the sparse surveys
(Figs 5, 12, 14), ≈ 4–5 pt in Fig 9/10, ≈ 2.5–3 pt in the dense Fig 11;
hollow markers (`mfc='none'`, mew ≈ 1.2) ≈ 7 pt in Figs 6, 7, 13.

**Stacked panels [SOLID]:** Figs 6, 7, 11 use two vertically stacked axes
`sharex=True`, equal heights, small `hspace` (≈ 0.08), x-label only on the
bottom panel, separate y-labels; legend in the top panel only.

**Colorbars [SOLID]:** only in Figs 1–2: one tall vertical colorbar shared
by all six panels at the right, triangular extensions at both ends,
9 labelled ticks at multiples of 0.154 (−0.616 … 0.616), label rotated
vertical (`δU M/m•`).

Reference rcParams block (use for every reproduction of this paper):
```python
import matplotlib as mpl, numpy as np
mpl.rcParams.update({
    'text.usetex': True, 'font.family': 'serif',
    'text.latex.preamble': r'\usepackage{amsmath,amssymb}',
    'font.size': 10, 'axes.labelsize': 12, 'legend.fontsize': 9.5,
    'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'axes.grid': False, 'legend.fancybox': True, 'legend.framealpha': 0.8,
    'axes.prop_cycle': mpl.cycler(color=['#0f4c81','#e2583e','#009473',
                                         '#bf1932','#5f4b8b','#a47864']),
    'figure.figsize': (3.3, 2.4), 'savefig.bbox': 'tight',
})
PLASMA5 = mpl.colormaps['plasma'](np.linspace(0.1, 0.9, 5))
VIRIDIS4 = mpl.colormaps['viridis']([0.8, 0.6, 0.4, 0.2])
SET2_3 = mpl.colormaps['Set2'].colors[:3]
```

---

## 1. Fig 1 — `spiral_grid.png` (`fig:spirals`): linear disk response, 2×3 survey

**Caption (main.tex:853–858):** The linear response of various accretion
disks to an EMRI at $a_\bullet = 15M$, visualized via the resulting enthalpy
perturbation. The spiral arms are truncated at the ISCO (black dash-dot
line), and wind more tightly in thinner disks. In disks with $\partial_r
c_s<0$, the inner spiral arm winds more loosely, while the outer spiral arm
winds more tightly.

**Data / parameters (brief §4 row 1):** $\delta U\,M/m_\bullet$ (enthalpy
perturbation, synthesized from the m-modes, $m_{max}\gtrsim 100$) on the
equatorial plane $(r\cos\phi, r\sin\phi) \in [-40,40]M$; $a_\bullet=15M$;
columns $(\mathfrak{p},\mathfrak{q},h_0) = (0,0,0.15), (0,0,0.1),
(2,1,0.1)$; rows $\chi=-0.9$ (top), $\chi=+0.9$ (bottom). Domain
$r_{in}=r_{ISCO}$ (8.72M for χ=−0.9; 2.32M for χ=+0.9) to $r_{out}=5a_\bullet$
(clipped by the ±40 frame, the visible disk edge is the circle $r\simeq 40$).

**Recorded attributes [SOLID unless noted]:**
- Layout: `subplots(2, 3, sharex=True, sharey=True)`, square panels
  (`aspect='equal'`), very tight spacing (wspace/hspace ≈ 0.04), thin black
  frames, ticks in on all sides with minor ticks. x ticks −20, 0, 20 on the
  bottom row only; y ticks −20, 0, 20 on the left column only. Axis labels:
  x `$r\cos(\phi)$` (every bottom panel), y `$r\sin(\phi)$` (every left
  panel) — parentheses present in this figure.
- Field: `pcolormesh`/`imshow` of $\delta U M/m_\bullet$ with `cmap='RdBu'`,
  `vmin=-0.77, vmax=+0.77` [PRELIMINARY on 0.77; ticks are
  ±0.154·k, k=0..4 hence vmax = 5×0.154 = 0.77], colorbar `extend='both'`,
  ticks labelled to 3 decimals (0.616, 0.462, 0.308, 0.154, 0.000, −0.154 …
  −0.616). Colorbar: single tall bar at right spanning both rows, label
  `$\delta U\,\frac{M}{m_\bullet}$` rotated 90° (vertical, reads bottom→top),
  placed right of the bar.
- Region outside the plotted disk ($r>40$ corners) is pure white; the disk
  interior background is near-white with a faint pale-pink/blue mottling
  (near-zero field). Blue = positive (compression) spiral arm, a thin red
  trough on its inner side.
- Centre: small **black filled disk** of radius $r_+ = M(1+\sqrt{1-\chi^2})
  \approx 1.44M$ (horizon) — appears as a dot ≈ 1.5M in radius [PRELIMINARY
  identification as $r_+$]. **Black dash-dot circle** at $r_{ISCO}$ (clearly
  visible at ≈ 8.7M in the top row; tiny at 2.3M in the bottom row, hugging
  the black dot). Legend entry `$r_\mathrm{ISCO}$` with a dash-dot black
  handle, lower-left of panel (0,0) only, unframed.
- Per-panel text stamp top-left: `$(\mathfrak{p},\mathfrak{q},h_0,\chi)=
  (0,\,0,\,0.15,\,-0.9)$` etc., ≈ 10 pt, inside axes, roman with math.
- Visible physics structure: ONE dominant trailing spiral arm launched at
  the EMRI (on the +x axis at $r=15$: a dark-blue wake blob at (15, 0)),
  outer arm winds outward; inner arm winds inward and is truncated at the
  ISCO circle (top row) or continues as many tight windings down to the
  horizon (bottom row, χ=0.9 where ISCO is at 2.3M). Winding is tighter in
  thinner disks (middle/right columns have 2–3 turns inside r=40 vs ~1 turn
  for h0=0.15). Right column (q=1): inner arm looser, outer arm tighter.

**Reproduction recipe:**
```python
fig, axs = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True, sharey=True,
                        gridspec_kw=dict(wspace=0.04, hspace=0.04))
for ax, (p,q,h0,chi) in zip(axs.flat, cases):
    im = ax.pcolormesh(X, Y, dU*M/m, cmap='RdBu', vmin=-0.77, vmax=0.77,
                       rasterized=True, shading='auto')
    ax.add_patch(Circle((0,0), r_plus, color='k'))            # horizon
    ax.add_patch(Circle((0,0), r_isco, fill=False, ls='-.', color='k', lw=1.2))
    ax.set_aspect('equal'); ax.set_xlim(-40,40); ax.set_ylim(-40,40)
    ax.text(0.03, 0.96, rf'$(\mathfrak{{p}},\mathfrak{{q}},h_0,\chi)=({p},\,{q},\,{h0},\,{chi})$',
            transform=ax.transAxes, va='top')
axs[0,0].plot([], [], 'k-.', label=r'$r_\mathrm{ISCO}$'); axs[0,0].legend(loc='lower left', frameon=False)
cb = fig.colorbar(im, ax=axs, extend='both', pad=0.02, fraction=0.03,
                  ticks=np.linspace(-0.616, 0.616, 9))
cb.set_label(r'$\delta U\,\frac{M}{m_\bullet}$')
```
Save as PNG (~500 dpi; the original is a raster).

**VLM MATCH CHECKLIST — Fig 1:**
- [ ] 2×3 square panels, shared axes, one shared RdBu colorbar with ±arrow ends and 0.154-spaced ticks
- [ ] ISCO dash-dot circle + central black horizon dot; r_ISCO legend in panel (0,0) only
- [ ] one trailing spiral launched at (+15, 0); blue crest with red inner trough
- [ ] top row (χ=−0.9): inner arm terminates at ISCO r≈8.7; bottom row (χ=+0.9): many tight inner windings down to ≈2M
- [ ] tighter winding in h0=0.1 vs 0.15; q=1 column: looser inner / tighter outer arm
- [ ] parameter stamp top-left of every panel; `r cos(φ)` / `r sin(φ)` labels; ticks −20,0,20

---

## 2. Fig 2 — `spirals_two_rows.png` (`fig:spirals_conv`): convergence of the spiral with $m_{max}$

**Caption (main.tex:861–866):** The convergence of the spiral patterns in a
$\mathfrak{p}=2,\mathfrak{q}=1,h_0=0.1$ disk with respect to the maximum
azimuthal mode number $m$ for an EMRI located at $a_\bullet = 15M$. The top
(bottom) panel shows the spiral waves for $\chi=-0.9$ ($\chi=0.9$). Fourier
ringing and blurring of the co-orbital region are obvious for in the first
two columns. In the $\chi=0.9$ case … the single spiral arm decoheres into
secondary and tertiary spirals … bottom-right panel …

**Data / parameters:** same field as Fig 1 for $(\mathfrak{p},\mathfrak{q},
h_0)=(2,1,0.1)$, $a_\bullet=15M$; columns $m_{max}=10, 30, 120$; top row
$\chi=-0.9$ with frame $|x|,|y|\le 40$; bottom row $\chi=+0.9$ **zoomed** to
$|x|,|y|\le 10$ (NOT shared y between rows).

**Recorded attributes [SOLID]:**
- Layout 2×3 as Fig 1 but rows have different limits: top ±40 (ticks −40,
  −20, 0, 20, 40 on y; −20, 0, 20 on x), bottom ±10 (ticks −10,−5,0,5,10 on
  y; −5, 0, 5 on x). Column titles above the top row: `$m_\mathrm{max}=10$`,
  `$m_\mathrm{max}=30$`, `$m_\mathrm{max}=120$`. Axis labels `$r\cos\phi$`,
  `$r\sin\phi$` (NO parentheses in this figure).
- Same RdBu colorbar (±0.77, 0.154 ticks, extend both), label written
  `$\delta U\,M/m_\bullet$` (slash form, not \frac).
- Parameter stamp only in the LEFT column: top-left
  `$(\mathfrak{p},\mathfrak{q},h_0,\chi)=(2,\,1,\,0.1,\,-0.9)$` (top) and
  `…(2,1,0.1,0.9)` (bottom) — with a white background box (visible in the
  bottom-left panel where it masks the field). `$r_\mathrm{ISCO}$` dash-dot
  legend lower-left of panel (0,0) only.
- Bottom row (zoom): large black filled central disk (horizon $r_+\approx
  1.44M$ fills ≈ 15% of the half-width), ISCO dash-dot circle at 2.32M just
  outside it, concentric tight windings (blue/red alternating rings) between
  r≈2.3 and 6, the main arm leaving the frame at top-right; in the $m_{max}
  =10$ panel strong Fourier ringing (concentric ripples, blurred
  co-orbital arc); at 120 the arm is crisp and splits into secondary/tertiary
  crests (multiple peaks/troughs) in the inner region.
- Top row: m_max=10 shows faint radial ripples and a blurred co-orbital
  region near (15,0); 30 and 120 nearly identical.

**Reproduction recipe:** as Fig 1 with per-row `set_xlim/ylim` (±40 / ±10),
`axs[0,j].set_title(rf'$m_\mathrm{{max}}={mmax}$')`, stamp only on
`axs[:,0]` with `bbox=dict(fc='white', ec='none', pad=1)`.

**VLM MATCH CHECKLIST — Fig 2:**
- [ ] column titles m_max = 10/30/120; top row ±40 frame, bottom row ±10 zoom
- [ ] ringing/blur visible at m_max=10, gone at 120; bottom-right shows split multi-crest inner spiral
- [ ] black horizon disk + ISCO dash-dot circle (large in bottom row)
- [ ] shared RdBu colorbar with arrows, label `δU M/m•`
- [ ] parameter stamps (white-boxed) in left column only

---

## 3. Fig 3 — `dT_dr_vs_dFdr.pdf` (`fig:angular_momentum_balance`): torque density vs AMF density, with inset

**Caption (main.tex:903–908):** Angular momentum balance in disk-EMRI
interactions. The red-dashed curves are used to denote the AMF density and
the blue solid curves are the torque density. Observe that they are perfectly
on top of each other except near corotation $r\sim a_\bullet$ where angular
momentum is deposited. The inset shows the spike in the AMF density near
corotation.

**Data / parameters (brief row 3):** $d\mathfrak{T}/dr$ (eq. 72b integrand
summed over m) and $dF_J/dr$ (eq. 72a) in units $F_{J,1}=F_{J,0}/a_\bullet$
vs $r/a_\bullet$; $(\mathfrak{p},\mathfrak{q})=(0,0)$, $\chi=0$, $h_0=0.05$,
$a_\bullet=15M$ (ISCO at $6M/15M = 0.4$).

**Recorded attributes [SOLID]:**
- Single axes, x `$r/a_\bullet$` linear 0.4–3.0 (ticks 0.5 … 3.0 step 0.5,
  minor ticks), y label `Torque Balance` (roman text), linear −7 … 8
  (ticks −5.0, −2.5, 0.0, 2.5, 5.0, 7.5).
- Series: **blue solid** `#0f4c81` lw≈2 `$d\mathfrak{T}/dr\ (F_{J,1})$`;
  **orange-red dashed** `#e2583e` lw≈1.8, long dashes,
  `$dF_J/dr\ (F_{J,1})$` drawn on top of the blue so both are visible as an
  alternating blue/red line; **black dash-dot vertical line** at $r/a_\bullet
  = 0.4$ (`ISCO`, `axvline(ls='-.', color='k')`).
- Shape: flat ≈0 with small wiggles for r/a• 0.4–0.85 and 1.2–3 (tiny
  bumps near 0.5, 0.6, 1.7, 2.1, 2.5); large antisymmetric double peak around
  corotation: trough ≈ −6.3 at r/a•≈0.96, crest ≈ +6.6 at ≈1.04; the red
  dashed AMF density has an additional **narrow vertical spike at r/a•=1.00**
  shooting off the top of the axes (clipped at y≈8).
- Legend: upper-right, **unframed**, three entries (blue solid, red dashed,
  black dash-dot).
- **Inset** (`inset_axes`/`fig.add_axes([0.55, 0.17, 0.35, 0.33])` in
  axes-fraction coordinates, lower-right quadrant of the main axes, thin
  black frame, ticks in): x 0.99–1.01 (ticks 0.99, 1.00, 1.01), y ≈ −3 … 48
  (ticks 0, 20, 40); shows the blue/red overlapped line as a gently rising
  straight line through 0 and the **red dashed spike** reaching ≈ 45 at
  exactly 1.000 (width ≈ 2×10⁻⁴ a•, the Landau-regularised corotation
  delta). The inset has no title, its fonts are slightly smaller (≈ 8 pt).

**Reproduction recipe:**
```python
fig, ax = plt.subplots(figsize=(3.3, 2.4))
ax.plot(r, dTdr, color='#0f4c81', lw=2.0, label=r'$d\mathfrak{T}/dr\ (F_{J,1})$')
ax.plot(r, dFdr, color='#e2583e', lw=1.8, ls=(0,(5,3)), label=r'$dF_J/dr\ (F_{J,1})$')
ax.axvline(r_isco/a, color='k', ls='-.', lw=1.2, label='ISCO')
ax.set_xlim(0.3, 3.0); ax.set_ylim(-7, 8); ax.set_xlabel(r'$r/a_\bullet$'); ax.set_ylabel('Torque Balance')
ax.legend(loc='upper right', frameon=False)
axin = ax.inset_axes([0.55, 0.15, 0.36, 0.34]); axin.plot(...same two series...)
axin.set_xlim(0.99, 1.01); axin.set_ylim(-3, 48); axin.set_xticks([0.99,1.0,1.01]); axin.set_yticks([0,20,40])
```

**VLM MATCH CHECKLIST — Fig 3:**
- [ ] blue solid and red dashed coincide everywhere except a red-only spike at r/a•=1
- [ ] trough ≈ −6 at 0.96, crest ≈ +6.5 at 1.04 (antisymmetric S-shape), flat wiggly tails
- [ ] black dash-dot ISCO line at 0.4; unframed 3-entry legend upper-right
- [ ] inset lower-right, x 0.99–1.01, y to ≈45, spike at 1.000

---

## 4. Fig 4 — `dT_dr_versus_ap.pdf` (`fig:dtdr_ap`): torque-density profiles vs $a_\bullet$

**Caption (main.tex:915–922):** Torque density profiles near corotation for
various $a_\bullet/M$ for $\chi=\mathfrak{p}=\mathfrak{q}=0$. In the nearly
Newtonian limit ($a_\bullet=1000M$) the peaks of the torque density are
offset by about $\pm h_0 a_\bullet$ from corotation; in the strongly
relativistic case ($a_\bullet=10M$), the peaks of the torque density shift
inward by $\sim a_\bullet h_0/5$ and increase in magnitude by a factor of
$\sim 1.5$.

**Data / parameters (brief row 4):** $d\mathfrak{T}/dr$ in $F_{J,1}$ vs
$(r/a_\bullet-1)/h_0 \in [-2,2]$; $a_\bullet \in \{10,30,100,300,1000\}M$;
$\chi=\mathfrak{p}=\mathfrak{q}=0$, $h_0=0.05$.

**Recorded attributes [SOLID]:**
- Single axes; x `$(r/a_\bullet-1)/h_0$` linear −2…2 (ticks −2,−1,0,1,2);
  y `$d\mathfrak{T}/dr\,(F_{J,1})$` linear −7.5…7.5 (ticks −7.5…7.5 step
  2.5).
- 5 lines, no markers: $a_\bullet=10M$ light green `#7ad151` solid;
  30M teal `#22a884` dashed; 100M blue-teal `#2a788e` dash-dot; 300M indigo
  `#414487` dotted; 1000M **black thin solid** (lw≈0.8). lw ≈ 1.2 for the
  coloured ones. (= `viridis([0.8,0.6,0.4,0.2])` + black.)
- Shape: antisymmetric S-curve through the origin; 1000/300/100 almost
  coincide with extrema ≈ ±5.0 at x≈±1; 30M slightly larger (≈ ±5.7 at
  ≈±0.95); 10M much larger (≈ ±7.5) with peaks shifted inward to x≈±0.8 and
  the curve crossing the others at |x|≈1.4.
- Legend: upper-left, framed (fancybox, rounded, light-grey edge), **2
  columns**, entries `$a_\bullet$=10$M$` … `$a_\bullet$=1000$M$` (typeset
  as `$a_\bullet=10M$`), order 10, 30, 100 | 300, 1000.

**Reproduction recipe:** `for a, c, ls in zip([10,30,100,300], VIRIDIS4,
['-','--','-.',':']): ax.plot(x, dT, color=c, ls=ls, lw=1.2, label=...)`;
`ax.plot(x, dT1000, 'k-', lw=0.8, label=r'$a_\bullet=1000M$')`;
`ax.legend(loc='upper left', ncol=2)`; `ax.set_xlim(-2,2); ax.set_ylim(-7.5,7.5)`.

**VLM MATCH CHECKLIST — Fig 4:**
- [ ] 5 S-curves through origin; four near-coincident (±5), green 10M curve ±7.5 with inward-shifted peaks
- [ ] viridis-green→indigo colour order with solid/dashed/dash-dot/dotted; black thin 1000M
- [ ] framed 2-column legend upper-left; axes ±2 × ±7.5

---

## 5. Fig 5 — `T_GR_versus_T_Newt.pdf` (`fig:TGR_vs_TNewt`): relativistic/Newtonian torque ratio vs spin

**Caption (main.tex:927–932):** Ratio of relativistic torque to the Newtonian
torque for disk profiles for a scale height of $h_0=0.05$ and $a_\bullet=12M$
as a function of the spin of the SMBH. Observe that disks with negative
(positive) surface density slope $\mathfrak{p}>0$ ($\mathfrak{p}<0$) can
decrease (increase) the relativistic torque compared to the Newtonian
calculation.

**Data / parameters (brief row 5):** $T_{GR}/T_{Newt}$ vs $\chi \in
\{-0.9,-0.7,\dots,0.9\}$ (10 points) for $(\mathfrak{p},\mathfrak{q}) \in
\{(-1.5,0),(-1.5,1),(0,0),(1.5,0),(1.5,1)\}$; $h_0=0.05$, $a_\bullet=12M$;
Newtonian denominator from the Newtonian-limit solver (eq. 50).

**Recorded attributes [SOLID]:**
- Single axes; x `$\chi$` −1…1 (ticks −0.5, 0.0, 0.5), y
  `$T_\mathrm{GR}/T_\mathrm{Newt}$` 0.15–1.6 (ticks 0.25 … 1.50).
- 5 marker-only series, large filled markers (≈ 8 pt):
  (−1.5, 0.0) blue `#0f4c81` circle `o`; (−1.5, 1.0) orange-red `#e2583e`
  star `*`; (0.0, 0.0) green `#009473` diamond `D`; (1.5, 0.0) red `#bf1932`
  triangle-down `v`; (1.5, 1.0) purple `#5f4b8b` square `s`.
- Legend lower-right, framed fancybox, entries written
  `$\mathfrak{p}=-1.5,\mathfrak{q}=0.0$` etc. (one decimal).
- Trends: p=−1.5 series ≈1.4–1.5 decreasing mildly with χ (stars 1.5→1.32,
  circles 1.37→1.15, circles peak at χ≈−0.7); p≥0 series below 1, rising
  with χ (diamonds 0.58→0.93; squares 0.55→0.89; triangles 0.33→0.85).

**Reproduction recipe:** `ax.plot(chi, ratio, marker, ms=8, ls='none',
color=c, label=...)` for the 5 series; `ax.legend(loc='lower right')`.

**VLM MATCH CHECKLIST — Fig 5:**
- [ ] two clusters: p=−1.5 (circle, star) at 1.15–1.5 falling with χ; p≥0 (diamond, triangle-down, square) at 0.3–0.95 rising with χ
- [ ] large filled markers, exact colour/marker pairing as above; framed legend lower-right

---

## 6. Fig 6 — `chi_rp.pdf` (`fig:torque_chi_ap`): total torque vs $a_\bullet$ for three spins, 2 panels

**Caption (main.tex:934–939):** Variation in the total torque as a function
of $a_\bullet$ for different values of $\chi$ for a disk of scale height
$h_0=0.05$. Relativistic effects are important as we approach $a_\bullet\sim
\mathcal{O}(50M)$; the retrograde orbits ($\chi<0$) face stronger relativistic
effects since the ISCO is closer to the EMRI orbit compared to the prograde
($\chi>0$) and the non-spinning ($\chi=0$) cases.

**Data / parameters (brief row 6):** $T$ in $F_{J,0}h_0$ vs $a_\bullet/M \in
[10,1000]$ (≈ 40 log-spaced points by visual count [PRELIMINARY]), $\chi \in \{0.9, 0, -0.9\}$; top panel
$(\mathfrak{p},\mathfrak{q})=(1.5,1)$, bottom $(0,0)$; $h_0=0.05$.

**Recorded attributes [SOLID]:**
- Two stacked panels `sharex=True`, equal heights; x **log** `$a_\bullet/M$`
  $10^1$–$10^3$ (ticks 10¹, 10², 10³ + log minor ticks); y linear
  `$T\,(F_{J_0}h_0)$` (note subscript written `F_{J_0}` here): top 0.8–2.1
  (ticks 1.0, 1.5, 2.0), bottom 0.4–1.2 (ticks 0.50, 0.75, 1.00).
- Hollow markers only (`mfc='none'`, mew≈1.3, ms≈7): $\chi=0.9$ red
  `#e70234` square `s`; $\chi=0.0$ purple `#580772` star `*` (the 5-point
  star drawn hollow); $\chi=-0.9$ blue `#2584da` circle `o`.
- Panel stamps top-left: `$\mathfrak{p}=1.5,\ \mathfrak{q}=1.0$` (top),
  `$\mathfrak{p}=0.0,\ \mathfrak{q}=0.0$` (bottom). Legend framed fancybox
  in the top panel, right-centre/lower-right (`loc='center right'`), entries
  `$\chi=0.9$`, `$\chi=0.0$`, `$\chi=-0.9$`.
- Trends: all three converge to a plateau (≈1.97 top, ≈1.13 bottom) for
  $a_\bullet\gtrsim 100M$; at small $a_\bullet$ the blue (χ=−0.9) drops
  steeply (to 0.8 top / 0.4 bottom at 10–11M), purple intermediate, red
  highest (1.6–1.8 top / 1.0 bottom at 10M) with slight scatter at 10–30M.

**Reproduction recipe:**
```python
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(3.3, 2.6), gridspec_kw=dict(hspace=0.08))
for chi, c, mk in [(0.9,'#e70234','s'), (0.0,'#580772','*'), (-0.9,'#2584da','o')]:
    ax1.plot(a, T_top[chi], mk, mfc='none', mec=c, mew=1.3, ms=7, label=rf'$\chi={chi}$')
ax1.set_xscale('log'); ax1.legend(loc='center right'); ax1.text(0.03,0.93, r'$\mathfrak{p}=1.5,\ \mathfrak{q}=1.0$', transform=ax1.transAxes, va='top')
```

**VLM MATCH CHECKLIST — Fig 6:**
- [ ] 2 stacked panels, log-x 10–1000, hollow square/star/circle in red/purple/blue
- [ ] plateau at large a•, blue circles drop lowest at a•→10M, red highest
- [ ] stamps p,q top-left of each panel; framed legend in the top panel only

---

## 7. Fig 7 — `rp_h0.pdf` (`fig:torque_ap_h0`): total torque vs $a_\bullet$ for five $h_0$, 2 panels

**Caption (main.tex:944–949):** Variation in the torque with $h_0$ and
$a_\bullet$ for disks around a SMBH of spin $\chi=0$. Thicker disks are
impacted strongly by the inner regions of the disk leading to some
oscillatory behavior. Observe that the torque monotonically increases until
$h_0=0.125$ and then decreases for very thick disks $h_0=0.15$.

**Data / parameters (brief row 7):** as Fig 6 but $\chi=0$ and $h_0 \in
\{0.05, 0.075, 0.1, 0.125, 0.15\}$; top $(1.5,1)$, bottom $(0,0)$.

**Recorded attributes [SOLID]:**
- Same 2-panel layout/axes as Fig 6; y top 1.25–2.15 (ticks 1.25, 1.50,
  1.75, 2.00), bottom 0.8–1.3 (ticks 1.0, 1.2).
- Hollow markers, plasma colours: $h_0=0.050$ dark purple `#41049d` square
  `s`; 0.075 purple `#8f0da4` star `*`; 0.100 magenta-pink `#cc4778` circle
  `o`; 0.125 orange `#f2844b` triangle-right `>`; 0.150 yellow `#fcce25`
  hexagon `h`. Legend in the top panel, right-centre, framed, labels
  `$h_0=0.050$` … `$h_0=0.150$` (three decimals), order 0.050 → 0.150.
- Trends: top panel: all rise from ≈1.4–1.5 at 10M to ≈1.95–2.05; yellow
  hexagons highest and oscillating at 20–40M; bottom panel: plateaus
  ordered 0.125 (≈1.21) > 0.1 (≈1.18) > 0.15 ≈ 0.075 ≈ 0.05 (≈1.14–1.15);
  yellow/orange oscillatory at 15–40M.

**VLM MATCH CHECKLIST — Fig 7:**
- [ ] plasma purple→yellow hollow markers (square, star, circle, >, hexagon)
- [ ] bottom panel plateau ordering: orange (0.125) on top, pink next, yellow/purples bunched below
- [ ] oscillations of thick-disk series at a• ≈ 15–40M; stamps p,q; framed legend top panel

---

## 8. Fig 8 — `fig_pq_m95h0_violin.pdf` (`fig:m_95_plot`): violin plot of $m_{95}h_0$

**Caption (main.tex:958–963):** Mode number at which 95% of the total
cumulative torque is attained for different values of
$(\mathfrak{p},\mathfrak{q})$. The distribution over $h_0$ and $\chi$ is
obtained by scanning $h_0\in[0.05,0.15]$ and $\chi\in[-0.9,0.0,0.9]$.
Observe that $m_{95}h_0$ is approximately constant and the largest variation
is observed when the EMRI is closest to the SMBH.

**Data / parameters (brief row 8, main.tex:992):** for each
$(\mathfrak{p},\mathfrak{q}) \in \{(0,0),(1.5,0),(2,0),(1.5,1),(2,1)\}$ and
each $a_\bullet\in\{15,30,120\}M$, the sample of $m_{95}h_0$ over the scan
$h_0\in[0.05,0.15]$ (several values) × $\chi\in\{-0.9,0,0.9\}$;
$m_{95}$ = smallest m at which the cumulative torque reaches 95% of the
total; $b=0.6$.

**Recorded attributes — VIOLIN SPECIFICS [SOLID]:**
- Single axes; x categorical, 5 group centres at 0…4 labelled
  `$(0,0)$, $(1.5,0)$, $(2,0)$, $(1.5,1)$, $(2,1)$`, x-label
  `$(\mathfrak{p},\mathfrak{q})$`; y `$m_{95}h_0$` linear **2.00–3.50**
  (ticks 2.00, 2.25, … 3.50; two decimals). Ticks in, minor ticks on y.
- **Light-grey vertical gridlines** at the 5 group centres
  (`ax.xaxis.grid(True, color='0.85', lw=0.8)` / `ax.grid(axis='x')`),
  drawn behind the violins; no horizontal grid.
- Per group, **three violins** offset at x = centre + {−0.3, 0, +0.3}, each
  `widths≈0.25`: $a_\bullet=15M$ Set2 teal `#66c2a5` (left), 30M orange
  `#fc8d62` (middle), 120M blue-grey `#8da0cb` (right). Bodies filled,
  `alpha≈0.8–0.9` (they appear slightly lighter than pure Set2, e.g. measured
  `#7ccbb2`, `#fc9d79`, `#9daed2`), with a **thin black edge**
  (`pc.set_edgecolor('k'); pc.set_linewidth(0.6)`).
- Inner statistics: `violinplot(showmeans=True, showextrema=True,
  showmedians=False)` styled black: a thin black vertical line spanning
  min→max (`cbars`), short black horizontal caps at min and max
  (`cmins`/`cmaxes`, width ≈ violin half-width) and a short black horizontal
  tick at the mean (`cmeans`) which with the vertical bar reads as a small
  "+" inside each body. KDE bodies (default Gaussian KDE, `points=100`).
- Group shapes: (0,0) small and low (2.15–2.5, means ≈2.3/2.25/2.24 for
  15/30/120M); (1.5,0) taller (teal 2.3–3.0 mean 2.55; orange 2.3–2.6;
  blue 2.25–2.45); (2,0) tallest teal (2.25–3.15, mean 2.65), orange
  2.25–2.9, blue 2.15–2.6; (1.5,1) high & narrow (teal 2.75–3.15 mean 2.94;
  orange 2.72–2.93; blue 2.7–2.83); (2,1) teal 2.75–3.2 mean 2.96, orange
  2.72–3.05, blue 2.6–2.9. Teal (15M) always widest spread; blue (120M)
  narrowest — the "largest variation closest to the SMBH".
- Text stamp top-left: `distribution over $h_0,\chi$` (roman + math).
- Legend upper-right, **unframed**, handles = rounded filled patches
  (`Patch`/`FancyBboxPatch` with `boxstyle='round'`, no edge) in the three
  Set2 colours, labels `$a_\bullet=15M$`, `$a_\bullet=30M$`,
  `$a_\bullet=120M$`.

**Reproduction recipe:**
```python
fig, ax = plt.subplots(figsize=(3.3, 2.5))
cols = mpl.colormaps['Set2'].colors[:3]; off = [-0.3, 0, 0.3]
for j, (ap, c) in enumerate(zip([15,30,120], cols)):
    parts = ax.violinplot([samples[pq][ap] for pq in PQ], positions=np.arange(5)+off[j],
                          widths=0.25, showmeans=True, showextrema=True)
    for b in parts['bodies']: b.set_facecolor(c); b.set_edgecolor('k'); b.set_linewidth(0.6); b.set_alpha(0.85)
    for k in ('cbars','cmins','cmaxes','cmeans'): parts[k].set_color('k'); parts[k].set_linewidth(0.8)
ax.set_xticks(range(5)); ax.set_xticklabels([r'$(0,0)$', r'$(1.5,0)$', r'$(2,0)$', r'$(1.5,1)$', r'$(2,1)$'])
ax.set_xlabel(r'$(\mathfrak{p},\mathfrak{q})$'); ax.set_ylabel(r'$m_{95}h_0$'); ax.set_ylim(2.0, 3.5)
ax.xaxis.grid(True, color='0.85', lw=0.8); ax.set_axisbelow(True)
ax.text(0.03, 0.95, r'distribution over $h_0,\chi$', transform=ax.transAxes, va='top')
from matplotlib.patches import FancyBboxPatch
handles=[mpl.patches.Patch(facecolor=c, edgecolor='none', label=rf'$a_\bullet={ap}M$') for ap,c in zip([15,30,120],cols)]
ax.legend(handles=handles, loc='upper right', frameon=False, handlelength=1.6)
```

**VLM MATCH CHECKLIST — Fig 8:**
- [ ] 5 groups × 3 violins (teal/orange/blue-grey Set2, left→right = 15/30/120M), black edges, black min–max bar with caps and a mean tick ("+")
- [ ] y range 2.0–3.5; (0,0) lowest ≈2.25, (1.5,1)/(2,1) highest ≈2.75–3.2; teal always widest
- [ ] vertical grey gridlines at group centres only; stamp "distribution over h0, χ" top-left; unframed patch legend upper-right

---

## 9. Fig 9 — `num_vs_analytical.pdf` (`fig:num_vs_analytical`): analytical/numerical Lindblad torque ratio

**Caption (main.tex:965–970):** Comparison between the Lindblad torque
obtained from linear theory and analytical results from
[HegadeKR:2025dur, HegadeKR:2025rpr] in Schwarzschild spacetime for
$h_0=0.05$. The analytical torque cutoff parameter has been set to
$m_\mathrm{max}=m_{45}$. The difference between the analytical results and the
numerical result increases as we approach the ISCO.

**Data / parameters (brief row 9):** ratio $T_L^{analytical}/T_L^{num}$ vs
$a_\bullet \in [12,100]M$ (≈19 points by visual count [PRELIMINARY], step ≈ 4.9M: 12.0, 17.3, …, 100),
$(\mathfrak{p},\mathfrak{q})=(0,0)$ and $(0,1)$; $h_0=0.05$, $\chi=0$;
analytical torque = relativistic pressureless Lindblad torque of
arXiv:2510.03564 eqs (18b),(20) with $m_{max}=m_{45}\sim 0.35/h_0$
(eq. 80); numerical $T_L$ from eq. (75).

**Recorded attributes [SOLID]:**
- Single axes; x `$a_\bullet\ (M)$` linear 0–103 (ticks 20, 40, 60, 80,
  100), y `$T_\mathrm{L}$ (analytical/num)` linear 0.86–1.53 (ticks 0.9,
  1.0, …, 1.5).
- Two marker-only series, filled, ms≈5: (0,0) blue `#0f4c81` diamond `D`;
  (0,1) orange-red `#e2583e` filled plus `P`. Legend upper-right, unframed,
  labels `$(\mathfrak{p},\mathfrak{q})=(0,0)$`, `$(\mathfrak{p},\mathfrak{q})=(0,1)$`.
- Trends: blue 1.52 at 12M, 1.1 at 17M, monotonically → 0.92 at 100M
  (crosses 1.0 at ≈ 32M); orange 0.93 at 12M, small bump 0.96 at 17M, then
  slowly → 0.875 at 100M.

**VLM MATCH CHECKLIST — Fig 9:**
- [ ] blue diamonds falling from 1.5 to 0.92; orange pluses 0.88–0.96 nearly flat, below blue
- [ ] unframed legend upper-right; axes 0–100 × 0.86–1.53

---

## 10. Fig 10 — `torque_reversal_in_extreme_disks.pdf` (`fig:LR_torque_and_reversal`): per-mode Lindblad torque

**Caption (main.tex:1007–1012):** Fraction of negative per-mode Lindblad
torque for different disk profiles in Schwarzschild spacetime with
$a_\bullet=15M$. The torques for $\mathfrak{p}=0$ is always positive while
the lower m-mode yield a negative torque for $\mathfrak{p}=5$. Observe that
the number of modes with negative torque density increase as we lower $h_0$.

**Data / parameters (brief row 10):** $T_{L,m}$ (per-mode Lindblad torque,
eq. 75 per m) in $F_{J,0}h_0$ vs $m \in [4, 40]$ (integer m, first plotted
point m=4), for $(\mathfrak{p},h_0) \in \{(5,0.03),(5,0.05),(0,0.03),
(0,0.05)\}$; $\chi=0$, $\mathfrak{q}=0$, $a_\bullet=15M$.

**Recorded attributes [SOLID]:**
- Single axes; x `$m$` 0–40 (ticks 0,10,20,30,40), y
  `$T_{\mathrm{L},m}\,(F_{J,0}h_0)$` −0.125…0.145 (ticks −0.10 … 0.10 step
  0.05). **Black solid horizontal line at y=0** (`axhline(0, color='k',
  lw=1.5)`).
- Four marker-only series, small filled markers (ms≈4): (p=5,h0=0.03) blue
  `#0f4c81` triangle-down `v`; (p=5,h0=0.05) orange-red `#e2583e`
  triangle-up `^`; (p=0,h0=0.03) purple `#5f4b8b` square `s`; (p=0,h0=0.05)
  brown `#a47864` diamond `D`.
- Legend: top, **2 columns**, unframed, labels `$\mathfrak{p}=5,h_0=0.03$`,
  `$\mathfrak{p}=5,h_0=0.05$` | `$\mathfrak{p}=0,h_0=0.03$`,
  `$\mathfrak{p}=0,h_0=0.05$` (column-major order as listed).
- Trends: brown diamonds start 0.14 at m=4 decaying to 0.005 at m=40;
  purple squares ≈0.085 at m=4–6, decaying to 0.01; orange triangles-up
  from −0.12 at m=4 cross zero at m≈10 and rise to ≈0.025 then decay to
  0.018; blue triangles-down from −0.085 at m=4 cross zero at m≈15, reach
  ≈0.02 by m≈25.

**VLM MATCH CHECKLIST — Fig 10:**
- [ ] black zero line; two positive decaying series (brown diamond, purple square) and two negative→positive series (orange ^ crossing ~m=10, blue v crossing ~m=15)
- [ ] small filled markers, no lines; unframed 2-column legend at top; x 0–40

---

## 11. Fig 11 — `convergence.pdf` (`fig:convergence`): per-mode torque convergence with resolution

**Caption (appendix.tex:68–73):** Convergence of the per-mode torque as
functions of resolution. The top (bottom) panel shows the absolute
(relative) errors from the finest resolution value. The parameters for the
convergence study are: $a_\bullet=25M, \chi=0, \mathfrak{p}=0,
\mathfrak{q}=0$ and $h_0=0.05$.

**Data / parameters (brief row 11):** $|T_m-T_{m,\mathrm{ref}}|$ (top) and
$|T_m-T_{m,\mathrm{ref}}|/|T_{m,\mathrm{ref}}|$ (bottom) vs $m \in [0,200]$
for $(n_L=n_R, n_C) = (30,3000), (60,6000), (120,12000)$, reference
$(240,240,24000)$.

**Recorded attributes [SOLID]:**
- Portrait figure (≈3.2×3.6 in), two stacked panels `sharex=True`, equal
  heights, hspace ≈ 0.06; x `$m$` 0–200 (ticks 0, 50, …, 200) labelled on
  the bottom only; **both y log**: top `$|T_m - T_{m,\mathrm{ref}}|$`
  10⁻¹⁵…10⁻⁶ (decade ticks 10⁻¹⁴, 10⁻¹², 10⁻¹⁰, 10⁻⁸, 10⁻⁶ labelled every
  2 decades), bottom `$|T_m-T_{m,\mathrm{ref}}|/|T_{m,\mathrm{ref}}|$`
  2×10⁻⁸…10⁻³ (labelled 10⁻⁷ … 10⁻⁴).
- **Light-grey major grid** in both panels (both axes), behind data.
- Three marker-only series, small filled markers (ms≈3): (30, 3000) blue
  `#0f4c81` star `*`; (60, 6000) orange-red `#e2583e` square `s`; (120,
  12000) green `#009473` circle `o`. Legend in the top panel upper-right,
  framed fancybox, labels `$n_L=n_R=30,\ n_C=3000$` etc.
- Trends: top — all decrease with m; blue forms a smooth band 10⁻⁷→10⁻¹²
  (flat shoulder 10⁻⁹ at m≈50–80, then steady decline), orange scattered
  10⁻⁸→10⁻¹³ with a kink up at m≈130, green 10⁻⁸→10⁻¹⁴ noisy (hitting the
  10⁻¹⁴–10⁻¹⁵ floor at m>150). Bottom — blue relative error rises from
  ≈10⁻⁵ (m<50) to 5×10⁻⁴ at 200 (smooth rising band), orange ≈10⁻⁶–10⁻⁵
  scattered then rising after m≈130 to 10⁻⁴, green noisy scatter
  3×10⁻⁸–10⁻⁵ with no trend.

**VLM MATCH CHECKLIST — Fig 11:**
- [ ] portrait 2-panel, both log-y with grey grid; x 0–200
- [ ] blue star smooth band (highest), orange square, green circle noisiest/lowest; top decreasing, bottom blue rising to ≈10⁻³.3
- [ ] framed legend top panel with n_L=n_R, n_C labels; small markers

---

## 12. Fig 12 — `pole_errs.pdf` (`fig:pole_disp_errors`): pole-displacement convergence

**Caption (appendix.tex:131–136):** Convergence of the total torque as a
function of the pole displacement parameter for different values of disk
thickness. Our reference solution is computed at $\epsilon_p=10^{-6}$. The
disk and EMRI parameters for this study are: $\mathfrak{p}=\mathfrak{q}=
\chi=0$ and $a_\bullet=100M$. Observe that the relative errors scale as
$\mathcal{O}(\epsilon_p)$ for all values of $h_0$ studied here.

**Data / parameters (brief row 12):** $|T-T_\mathrm{ref}|/T_\mathrm{ref}$
vs $h_0 \in \{0.15, 0.13, 0.11, 0.09, 0.07, 0.05\}$ for $\epsilon_p \in
\{10^{-3},10^{-4},10^{-5}\}$; ref $\epsilon_p=10^{-6}$; grid $n_L=n_R=240$,
$n_C=2400$ (appendix A:162).

**Recorded attributes [SOLID]:**
- Single axes; **x axis REVERSED**: `$h_0$` runs 0.15 (left) → 0.05 (right)
  (ticks 0.15, 0.13, 0.11, 0.09, 0.07, 0.05; `ax.invert_xaxis()`); y **log**
  `$|T-T_\mathrm{ref}|/T_\mathrm{ref}$` 6×10⁻⁶…3×10⁻³ (labelled 10⁻⁵, 10⁻⁴,
  10⁻³).
- Three marker-only series, **large filled markers** (ms≈10): $\epsilon_p=
  10^{-3}$ blue `#0f4c81` star `*`; $10^{-4}$ orange-red `#e2583e` square
  `s`; $10^{-5}$ green `#009473` circle `o`. Legend at the top spanning the
  width, **3 columns**, framed fancybox, labels typeset literally as
  `$\epsilon_p$=1e−03`, `…1e−04`, `…1e−05` (Python float formatting `1e-03`
  inside math, i.e. `rf'$\epsilon_p={eps:.0e}$'`).
- Trends: three nearly horizontal rows at ≈1.1×10⁻³, 1.2×10⁻⁴, 1.2×10⁻⁵,
  each drifting slightly down toward h0=0.05 (≈0.85×10⁻³, 0.9×10⁻⁴,
  0.85×10⁻⁵).

**VLM MATCH CHECKLIST — Fig 12:**
- [ ] reversed x (0.15→0.05), log y, three horizontal rows one decade apart
- [ ] big star/square/circle in blue/orange/green; framed 3-column legend on top with `1e-03` style labels

---

## 13. Fig 13 — `soft_differences.pdf` (`fig:variation_with_respect_to_b`): torque vs softening $b$

**Caption (appendix.tex:143–149):** The dependence of the torque on the
point-mass softening prescription for $\mathfrak{p}=\mathfrak{q}=\chi=0$,
and $a_\bullet=100M$. The total torque decreases as we increase the
softening parameter $b$.

**Data / parameters (brief row 13):** $T$ in $F_{J,0}h_0$ vs $b \in
\{0.3,0.4,\dots,0.8\}$ for $h_0 \in \{0.05,0.075,0.1,0.125,0.15\}$;
$\mathfrak{p}=\mathfrak{q}=\chi=0$, $a_\bullet=100M$.

**Recorded attributes [SOLID]:**
- Single axes; x `$b$` 0.3–0.8 (ticks 0.3…0.8), y `$T\,(F_{J,0}h_0)$`
  0.96–1.39 (ticks 1.0, 1.1, 1.2, 1.3).
- Hollow markers, plasma encoding identical to Fig 7 (0.050 purple square,
  0.075 purple star, 0.100 pink circle, 0.125 orange `>`, 0.150 yellow
  hexagon), ms≈7, mew≈1.3. Legend upper-right, framed fancybox, **2
  columns**, order 0.150, 0.125, 0.100 | 0.075, 0.050 (descending), labels
  `$h_0=0.150$` (three decimals).
- Trends: all decrease roughly linearly with b; orange (0.125) highest
  (1.39→1.045), pink (0.1) next (1.35→1.02), yellow (0.15) third
  (1.33→0.985), purple square and star (0.05, 0.075) coincide lowest
  (1.28→0.97).

**VLM MATCH CHECKLIST — Fig 13:**
- [ ] 5 decreasing hollow-marker series; orange > above pink o above yellow hexagon; purple square+star overlapping lowest
- [ ] framed 2-column legend upper-right in descending h0 order; axes 0.3–0.8 × 0.96–1.39

---

## 14. Fig 14 — `m95_vs_b.pdf` (`fig:m95_vs_b`): $m_{95}h_0$ vs softening $b$

**Caption (appendix.tex:153–160):** The dependence of $m_{95}$ on the
point-mass softening parameter $b$ for $\mathfrak{p}=\mathfrak{q}=\chi=0,
a_\bullet=30M$ for $h_0=0.05$. Decreasing $b$ increases the strength of the
point-mass potentially, leading an increase in $m_{95}$.

**Data / parameters (brief row 14):** $m_{95}h_0$ vs $b \in \{0.3,\dots,
0.8\}$ (6 points), single series; $\mathfrak{p}=\mathfrak{q}=\chi=0$,
$a_\bullet=30M$, $h_0=0.05$.

**Recorded attributes [SOLID]:**
- Single axes; x `$b$` 0.3–0.8, y `$m_{95}\,h_0$` 1.95–2.9 (ticks 2.0, 2.2,
  …, 2.8). One series: blue `#0f4c81` large filled circles (ms≈8), no line,
  no legend. Values ≈ 2.9, 2.7, 2.45, 2.25, 2.1, 1.95 (monotone decreasing,
  quantised because $m_{95}$ is an integer × 0.05).

**VLM MATCH CHECKLIST — Fig 14:**
- [ ] six large blue filled circles, monotone decreasing 2.9→1.95, no legend, no lines
- [ ] axes 0.3–0.8 × ≈1.95–2.9

---

## 15. Cross-figure consistency rules (apply to every reproduction)

1. Fonts: usetex CM serif + `\mathfrak` Fraktur for p, q. Ticks in, all
   four sides, minor ticks on. No grid unless Figs 8 (x-only) / 11 (both).
2. Colour encodings are SEMANTIC and fixed across figures: blue `#0f4c81`
   = first/base series, `#e2583e` = second; plasma(0.1–0.9) = $h_0$ ladder;
   viridis = $a_\bullet$ ladder; Set2 = $a_\bullet\in\{15,30,120\}$ in the
   violin; RdBu (unreversed, symmetric ±0.77) = $\delta U$ maps.
3. Survey figures are **marker-only** (no connecting lines); hollow markers
   when many points overlap (Figs 6, 7, 13), filled otherwise.
4. Torque unit labels are verbatim: `$T\,(F_{J,0}h_0)$` (Fig 10, 13; written
   `F_{J_0}` in Figs 6–7), `$d\mathfrak{T}/dr\,(F_{J,1})$` (Figs 3, 4).
5. Panel parameter stamps top-left inside the axes; legends framed-fancybox
   by default, unframed where noted.
6. Never fabricate the field or data to hit the reference: a reproduction
   figure is produced from GRTools pipeline output and then VLM-compared
   using the per-figure checklist above.
