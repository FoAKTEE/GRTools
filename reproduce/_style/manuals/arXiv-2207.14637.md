# FIGURE STYLE MANUAL — arXiv:2207.14637 (Fairbairn & Rafikov 2022, "Density waves in protoplanetary discs excited by eccentric planets: linear theory", MNRAS 517, 2121) — FR22

VLM-inspected from the arXiv PDF mirror
`/data/haiyangw/claude/emri-disk/ref-paper/arXiv-2207.14637/2207.14637.pdf`
(pages 2, 7, 8, 9, 10 rendered at 150 dpi, figure crops at 300 dpi) and the
EPS headers in `src/Figures/` (`0.1_sigma_maps.eps`, `0.25_sigma_maps.eps`,
`0.5_sigma_maps.eps`, `cartoon.eps`, `convergence_m.eps`, `convergence_e.eps`).
Captions are from `src/Sections/2_framework.tex` and `4_results.tex`. Mission
brief: `emri-disk/results/emri-eccentric/paper_EMRI-ECC/evidence/acquire-2207.14637/BRIEF.md`.

Tags: [SOLID] = read from file metadata or measured from pixels (colourbar
sampled against the matplotlib colormaps, hex codes of saturated pixels);
[PRELIMINARY] = estimated by eye from the render (sizes, line widths, curve
values).

READ THE RELEVANT SECTION OF THIS MANUAL **BEFORE** GENERATING ANY
REPRODUCTION FIGURE. After generating, VLM-compare against the reference per
the `vlm-figure-repro` skill (3-step recurring procedure).

---

## 0. Shared style conventions (whole paper)

**Toolchain [SOLID]:** Figs 2–6 are **matplotlib 3.3.4** EPS files
(`%%Creator` header). The text font is **DejaVu Sans**, embedded as a Type-3
font. There is NO LaTeX: math uses the matplotlib mathtext default (DejaVu
italic). For example, `$T = 0.00T_p$`, `$(\phi-\phi_p)/\pi$`, `$log_{10}r/a$`
(rendered all-italic, with no roman "log"), `$m_{\rm max}$` and
`$\log_{10}[\Delta(m_{\rm max}, m_{\rm max}+10)]$` (the y labels use a roman
"log"). Fig 1 is different: a **cairo 1.17.6** vector drawing (Inkscape-style)
with CM-like italic labels.
This paper is the ONLY one of the three eccentric papers drawn in DejaVu Sans.
Do not apply usetex here.

**Axes [SOLID]:** these are the plain matplotlib defaults. Ticks point
**outward**, on the left and bottom only. There are no minor ticks and no grid,
and the frame is a full black box. Exception: the polar-map row of Figs 2–4
has no box. Its first panel shows only a detached left and bottom spine, and
the other three panels have no axes at all (see §2).

**Colours [SOLID, measured]:** the paper uses pure named matplotlib colours:
`'k'` #000000, `'r'` #ff0000, `'b'` #0000ff (Figs 5–6). The wavelet overlays
are dashed `'cyan'` #00ffff and `'red'` #ff0000. The planet marker is a black
filled circle.

**Colormap [SOLID]:** **viridis**, symmetric limits **±2.0**. A colourbar
sample matches viridis with a mean ΔRGB of 0.0097 (the next best, cividis,
is 0.29). There are no extend arrows.

**Figure sizes [SOLID, EPS BoundingBox]:**
- map figures: 754 × 571 pt (10.5 × 7.9 in), included at `\textwidth`;
- convergence figures: 268–278 × 241–242 pt (≈ 3.8 × 3.35 in), included at `\columnwidth`;
- cartoon: 368 × 365 pt (square), included at `0.8\columnwidth`.

Reference rcParams block (matplotlib defaults, stated explicitly):
```python
import matplotlib as mpl
mpl.rcParams.update({
    'text.usetex': False, 'font.family': 'sans-serif',      # DejaVu Sans
    'mathtext.fontset': 'dejavusans',
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.top': False, 'ytick.right': False, 'axes.grid': False,
    'legend.frameon': True, 'legend.fancybox': True, 'legend.framealpha': 0.8,
})
CMAP_MAP = 'viridis'; VLIM = 2.0
```

---

## 1. Plot type A — geometry schematic (Fig 1)

- A square vector cartoon drawn with no axes frame.
- **Disc:** a pink/magenta radial-gradient disc, darkest at the centre and
  fading to white at the edge.
- **Outer boundary:** a grey dashed circle.
- **Axes arrows:** grey x and y arrows through the centre.
- **Orbit and bodies:**
  - black solid ellipse = the planet orbit, with the star at the origin;
  - black filled dot = star `$M_*$`, at the origin;
  - smaller black dot = planet `$M_p$`, on the ellipse at upper right.
- **Angles and radii:**
  - black dotted radial lines run to an arbitrary point (`$r$`, arrowhead at
    upper left) and to the planet (`$r_p$`);
  - **red** arcs mark `$\phi$` (large arc) and `$\phi_p$` (small arc);
  - a grey shaded wedge marks `$\phi_p$`.
- **Labels:** CM-style italics.
- **Recipe:** draw it as patches. It carries no data, so reproduction is
  optional.

**Visual structure to match:** the star is drawn at the origin, slightly
right of the ellipse centre (the cartoon is not to scale: the offset is
smaller than a focal offset would be). There are two red angle arcs, a dashed
outer circle on a pink gradient, and the planet on the upper-right part of
the ellipse.

---

## 2. Plot type B — 3 × 4 orbital-phase sequence of δΣ maps (Figs 2, 3, 4)

**Quantity:** `$(M_p/M_{th})^{-1}\,\delta\Sigma/\Sigma$` from linear theory
(sum over (m, l), m ≤ 150, |l−m| ≤ 8). The colourbar carries **no label**: the
normalisation is stated only in the caption [SOLID].

**Layout [SOLID]:**
- 3 rows × 4 columns, with columns = orbital phase `T = 0.00, 0.25, 0.50,
  0.75 T_p`;
- column titles above the top row only, mathtext `$T = 0.00T_p$` (two
  decimals, no space before `T_p`);
- ONE tall vertical colourbar at the far right, spanning all three rows, with
  ticks −2.0 … 2.0 every 0.5 (9 ticks, one decimal);
- small gaps between panels (wspace ≈ 0.08, hspace ≈ 0.15 [PRELIMINARY]).

**Row 1 — polar-geometry map:**
- **Frame:** a frame corotating with the planet's guiding centre. The planet's
  azimuth is therefore NOT fixed: it is on the +x axis at T = 0 and 0.5, above
  it at T = 0.25 and below it at T = 0.75.
- **Projection:** Cartesian x, y in units of a, square aspect.
- **Disc annulus:** drawn from r = 0.1 (a small white hole at the centre) to
  **r = 2** (outer edge). The region outside the disc is white.
- **Axes:**
  - only the first panel shows axes: a detached ("despined", offset) left and
    bottom spine with ticks −2 … 2 every 1, and no axis labels;
  - the other three panels are axis-off.
- **Planet marker:** black filled circle, ms ≈ 8 [PRELIMINARY], at
  (r_p cos φ_p, r_p sin φ_p).
  - At T = 0 it is at pericentre, x = a(1 − e) = 0.9 / 0.75 / 0.5 for
    e = 0.1 / 0.25 / 0.5.
  - At T = 0.5 it is at apocentre, x = 1.1 / 1.25 / 1.5.

**Row 2 — unwrapped map:**
- `pcolormesh` over x = `$log_{10}r/a$` ∈ [−0.5, 0.5] (ticks −0.4 … 0.4,
  step 0.2) and y = `$(\phi-\phi_p)/\pi$` ∈ [0, 2] (ticks 0.0 … 2.0, step 0.5).
- The planet sits at y = 0, so its marker is clipped to a half-disc on the
  bottom spine, at x = log10(r_p/a).
- y label on the first column only; x tick labels on every panel; NO x label
  on this row.
- Full black box frame.

**Row 3 — row 2 plus the ZZ22 wavelet overlay:**
- Identical map, with dashed lines on top (lw ≈ 2, dash ≈ 6/3
  [PRELIMINARY]):
  - **cyan** #00ffff = inward-propagating wavelet fronts;
  - **red** #ff0000 = outward-propagating fronts.
- x label `$log_{10}r/a$` under every panel of this row.

**Recipe:**
```python
fig = plt.figure(figsize=(10.5, 7.9))
gs = fig.add_gridspec(3, 5, width_ratios=[1,1,1,1,0.06])
for j, T in enumerate([0.0, 0.25, 0.5, 0.75]):
    ax = fig.add_subplot(gs[0, j]); ax.set_aspect('equal')
    ax.pcolormesh(X, Y, dS[j], cmap='viridis', vmin=-2, vmax=2, shading='auto', rasterized=True)
    ax.plot(xp[j], yp[j], 'ko', ms=8); ax.set_title(rf'$T = {T:.2f}T_p$')
    if j == 0: sns_like_despine(ax, offset=5)   # left/bottom spines only, ticks -2..2
    else: ax.axis('off')
    for i in (1, 2):
        ax = fig.add_subplot(gs[i, j])
        im = ax.pcolormesh(logr, phi/np.pi, dSu[j], cmap='viridis', vmin=-2, vmax=2, shading='auto')
        ax.plot(np.log10(rp[j]), 0, 'ko', ms=8, clip_on=True)
        if i == 2:
            for c, tr in [('cyan', inward[j]), ('red', outward[j])]:
                for seg in tr: ax.plot(*seg, '--', color=c, lw=2)
            ax.set_xlabel(r'$log_{10}r/a$')
        if j == 0: ax.set_ylabel(r'$(\phi - \phi_p)/\pi$')
fig.colorbar(im, cax=fig.add_subplot(gs[:, 4]), ticks=np.arange(-2, 2.01, 0.5))
```

**Dominant visual STRUCTURE (what a reproduction MUST show):**
- Background: the disc is uniformly teal (≈ 0, viridis mid #21918c). The arms
  are thin bright-yellow (+) crests with a dark-purple (−) trough on one side.
- **Inner disc, r ≲ 0.45a:** a dense set of tight, noisy concentric rings
  (unresolved high-m inner waves), plus a pronounced ring at r ≈ 0.45–0.5.
- **e = 0.1 (Fig 2):**
  - Row 1 shows ONE dominant trailing spiral launched at the planet. The
    outer arm runs clockwise and outward, reaching r ≈ 1.6–1.9 on the far
    side. The inner arm winds into the ring at r ≈ 0.45.
  - Rows 2 and 3 show the inner wake as a bright line rising from the planet
    up-left to (−0.45, ≈ 1), with 1–2 parallel inner lines near x ≈ −0.45 …
    −0.5. The outer wake enters from the top near x ≈ 0–0.15 and falls
    steeply to x ≈ 0.35 at y = 0.
  - At T = 0 a closed **loop / bifurcation** ("bubble") sits near
    (0.05–0.15, 1.7–2.0).
  - The red overlay is kinked: a shallow segment from (0, 2) to
    (0.17, 1.65), then steep down to (0.35, 0). A second red line follows
    x ≈ 0.48. The cyan overlay traces the main inner wake.
- **e = 0.25 (Fig 3):**
  - Row 1: at T = 0 there is a doubled (bifurcated) bright arm just outside
    the planet, two parallel arcs going down. Two or more outer arcs are
    visible at large r.
  - Rows 2 and 3: 3–4 parallel inner lines on the left. The outer wakes
    cross each other, and the red overlay has multiple kinks and jagged
    near-vertical segments at x ≈ 0.35–0.45.
- **e = 0.5 (Fig 4):**
  - Row 1: the wake amplitude is visibly weaker and fragmented.
    - At T = 0 the planet (0.5, 0) sits on the bright inner ring, and the
      strongest feature is a bright patch near (0.5, −1.1), DETACHED from the
      planet.
    - At T = 0.25 / 0.5 / 0.75 the planet carries only a faint thin arm or a
      diffuse blob.
  - Rows 2 and 3: about 6–8 parallel inner lines. The red overlay becomes a
    dense fan of many dashed lines between x ≈ 0.1 and 0.5, and the cyan
    overlay has several lines.
  - Planet x positions in rows 2 and 3: −0.30 at T = 0 and +0.18 at T = 0.5.
- **In every map:** spiral sense = trailing for counter-clockwise disc
  rotation (outer arm lags clockwise).

---

## 3. Plot type C — modal-convergence curves (Figs 5, 6)

**Layout [SOLID]:**
- A single axes, ≈ 3.8 × 3.35 in, full box, ticks out (left and bottom).
- Default framed legend (rounded corners, white face, grey edge) with labels
  `$e = 0.1$`, `$e = 0.25$`, `$e = 0.5$`.

**Series [SOLID]:** three solid lines with **'x' markers** (`'-x'`), lw ≈ 1.5,
ms ≈ 7 [PRELIMINARY]: black = e = 0.1, red = 0.25, blue = 0.5.

**Fig 5:**
- **x axis:** `$m_{\rm max}$`, linear, 0–140, ticks 0, 25, …, 125.
- **Data points:** 15 points at m_max = 0, 10, …, 140. The FIRST marker is
  drawn at x = 0 [SOLID, measured from the crop].
- **y axis:** `$\log_{10}[\Delta(m_{\rm max}, m_{\rm max}+10)]$`. This is a
  LINEAR axis holding log values: −6.7 … −2.7, ticks −6.5 … −3.0, step 0.5.
- **Legend:** upper right.

**Fig 6:**
- **x axis:** `$|l-m|_{\rm max}$` 0–7, 8 points at 0 … 7.
- **y axis:** `$\log_{10}[\Delta(|l-m|_{\rm max}, |l-m|_{\rm max}+1)]$`, a
  linear axis from −3.64 to −2.78 (ticks −3.6 … −2.8, step 0.1).
- **Legend:** lower left.

**Dominant visual STRUCTURE:**
- **Fig 5:**
  - All three curves start at ≈ −2.75 … −2.8 and fall by about 3–4 decades.
  - Black (e = 0.1) has a two-step shape: a steep drop between m_max = 50 and
    60 and a kink up at 130, ending lowest at ≈ −6.65 at 140.
  - Blue (0.5) drops fastest early: −4.2 at 40, then −5.0 at 80. It has a
    local bump up at m_max = 130 (−5.4).
  - Red (0.25) is the smoothest, ending ≈ −5.7.
  - Red lies above blue for m_max ≳ 30.
- **Fig 6:**
  - Black peaks at |l−m|_max = 1 (−2.82), then falls almost linearly to
    −3.6 at 7.
  - Red also peaks at 1 (−2.83), then decays slowly to −3.07.
  - Blue starts lowest (−3.21), peaks at 2 (−3.05), then decays slowly to
    −3.2.
  - The black curve crosses the blue between 3 and 4.

---

## 4. Per-figure table

Provenance classes: **LIN** = linear theory (reproducible from the paper's
equations with GRTools). **HEUR** = the ZZ22 Huygens-wavelet heuristic
(semi-analytic, not hydro). **SCHEM** = schematic. There is no hydro data and
no digitised overlay anywhere in FR22.

| Fig | file | panels | plot type | what is plotted | provenance per series | structure to check |
|---|---|---|---|---|---|---|
| 1 | cartoon.eps | 1 | A schematic | orbit geometry: star, ellipse, r, φ, r_p, φ_p | SCHEM | pink gradient disc, dashed outer circle, two red angle arcs, planet upper-right |
| 2 | 0.1_sigma_maps.eps | 12 (3×4) + cbar | B phase sequence | (M_p/M_th)⁻¹δΣ/Σ, e = 0.1, T = 0, .25, .5, .75 T_p | rows 1–2 LIN; row 3 LIN map + HEUR cyan/red dashed | one trailing arm; T = 0 loop near (0.1, 1.8) in rows 2–3; kinked red overlay; planet at log r = −0.046 (T = 0) / +0.041 (T = 0.5) |
| 3 | 0.25_sigma_maps.eps | 12 + cbar | B | same, e = 0.25 | as Fig 2 | bifurcated arm at T = 0; 3–4 inner lines; jagged red segments 0.35–0.45 |
| 4 | 0.5_sigma_maps.eps | 12 + cbar | B | same, e = 0.5 | as Fig 2 | weak, fragmented wake; bright patch detached from planet at T = 0; dense red fan; planet log r = −0.30 / +0.18 |
| 5 | convergence_m.eps | 1 | C convergence | log10 Δ(m_max, m_max+10) vs m_max, \|l−m\|_max = 8 | LIN (3 series k/r/b) | 15 points from 0 to 140; black ends lowest ≈ −6.65; blue bump at 130 |
| 6 | convergence_e.eps | 1 | C convergence | log10 Δ(\|l−m\|_max, +1) vs \|l−m\|_max, m_max = 150 | LIN (3 series) | black & red peak at 1, blue at 2; black steep linear fall to −3.6 |

**Per-figure VLM MATCH CHECKLIST**
- Fig 1:
  - [ ] square cartoon, no frame; ellipse + star at origin + planet dot;
    dotted radii; red φ, φ_p arcs; dashed boundary; pink radial gradient.
- Figs 2–4:
  - [ ] DejaVu Sans fonts (no CM).
  - [ ] Layout: 3 × 4 panels; phase titles `T = 0.00T_p` … `0.75T_p`; one
    tall right viridis colourbar ±2.0 (0.5 ticks, no label).
  - [ ] Axes:
    - polar row: only panel (0, 0) has detached left and bottom axes
      (−2 … 2), the rest are axis-off;
    - unwrapped rows: x ∈ [−0.5, 0.5], y ∈ [0, 2]; x labels only on row 3;
      y labels only on column 1.
  - [ ] Planet markers: black dot on the orbit in row 1; a half-dot on the
    bottom spine at log10(r_p/a) in rows 2–3.
  - [ ] Overlays and structure: cyan (inward) and red (outward) dashed
    overlays only in row 3; the per-e structure of §2.
- Figs 5–6:
  - [ ] Series and legend: black/red/blue `-x` lines; framed legend (upper
    right in Fig 5, lower left in Fig 6).
  - [ ] Axes: log10 values on a linear y axis; the exact abscissae (0 … 140
    step 10; 0 … 7).
  - [ ] Structure: crossing pattern and end values as in §3.

---

## 5. Defects and caveats found while inspecting (checked on the rendered pages)

1. **Fig 2 caption/panel-count mismatch [SOLID].** The caption (PDF p. 7;
   `4_results.tex` l.13) lists the phases "(0.0, 0.25, 0.50) T_p", but the
   figure has FOUR columns titled T = 0.00, 0.25, 0.50, 0.75 T_p. Figs 3–4
   inherit the caption. A reproduction must draw four columns.
2. **Fig 5 first abscissa [PRELIMINARY on its meaning].** The first plotted
   point sits at m_max = 0, so it reads as Δ(0, 10). The text only discusses
   m_max = 10 … 80. Reproduce the 15 abscissae 0 … 140 as drawn, and record
   how the m_max = 0 point is defined before gating against it.
3. **No colourbar label [SOLID]:** the normalisation (M_p/M_th)⁻¹ appears only
   in the caption. Do not add a label in a faithful reproduction.
4. **Wavelet overlay:** row 3 is NOT data from this paper's equations. It is
   the ZZ22 heuristic (wavelets launched every T/800 over 5 orbits; brief §
   Figures). Any reproduction of row 3 needs that construction as a separate
   node.
