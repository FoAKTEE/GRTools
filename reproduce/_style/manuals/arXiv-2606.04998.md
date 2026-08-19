# FIGURE STYLE MANUAL — arXiv:2606.04998 (Vu et al., "Self-force calculations with numerical relativity methods")

VLM-inspected from the reference PDFs at
`/data/haiyangw/claude/emri-disk/ref-paper/arxiv-2606.04998/src/`
(ssf_domain.pdf, matrix.pdf, precond.pdf, ssf_convergence.pdf, ssf_modes.pdf,
linsolv.pdf, ssf_error_budget.pdf). Captions from `main.tex` (lines 528–709).

READ THE RELEVANT SECTION OF THIS MANUAL **BEFORE** GENERATING ANY
REPRODUCTION FIGURE. After generating, VLM-compare against the reference PDF
per the `vlm-figure-repro` skill.

---

## 0. Shared style conventions (whole paper)

The figures split cleanly into two families:

**A. Quantitative data plots (Figs 4, 6, 7)** — clearly matplotlib:
- Typography: LaTeX Computer Modern serif everywhere (axis labels, ticks,
  legend), i.e. `text.usetex=True` or `mathtext.fontset='cm'` +
  `font.family='serif'`.
- Data colors: **grayscale only** (black → light gray). No color at all in
  the line plots.
- Markers: filled black circles `'o'` when iterations/discrete spins are
  plotted (Figs 6, 7); no markers on convergence-family curves (Fig 4).
- Reference lines: horizontal **black dashed** line for "Solver tolerance",
  labeled with plain text placed just below-left of the line *inside* the
  axes (no arrow, no box) — used identically in Figs 4 and 6.
- Grid: light-gray thin gridlines on (both axes), standard mpl z-order
  (below data).
- Legends: framed (rectangular box, white background, thin dark edge),
  placed upper-right inside the axes.
- Layout: vertical 2-panel stacks **sharing the x-axis** with the x-label
  only on the bottom panel; top (main) panel taller than the bottom
  (diagnostic) panel, roughly 3:1 (Fig 4) or 3:2 (Fig 6).
- Log y-axes almost everywhere; custom human tick labels with units where
  appropriate ("1 s", "10 s"; "Walltime: 0:00:06").
- Figure widths: single-column (`0.5\textwidth`) for Figs 4, 6; full
  two-column width (`\textwidth`) for Figs 1, 3, 5, 7.

**B. Schematics / 3D renders (Figs 1, 2, 3, 5)**:
- Fig 1 is a rendered 3D surface schematic, Fig 3 a pastel flowchart, Fig 2
  a B/W matrix-spy + domain sketch, Fig 5 triple 3D stacked wireframes.
- **Red is the reserved accent color** for the worldtube / regularized
  region across the paper: red rectangle + red text in Fig 1, red boxes and
  red dotted interface lines in Fig 2. Reuse red ONLY for that meaning.
- Serif/math (CM-style italic) labels even inside schematics (r, z², r₊,
  Ψ_m, x_p …); bold upright labels for the null coordinates **t, u, v**.
- White backgrounds throughout; no figure titles except the per-panel math
  titles of Fig 5.

---

## 1. Fig 1 — `ssf_domain.pdf` (`fig:domain`): computational domain

**Caption (main.tex:529):** Computational domain for a=0.5M, r0=10M, m=1
composed of discontinuous Galerkin elements. Domain from Kerr horizon r₊ to
outer boundary R=10⁴M, and z²=0 (equator) to z²=1 (poles). Plotted is |Ψ_m|;
in the regularized region (red rectangle) |Ψ_m^R| = |Ψ_m − Ψ_m^P| instead.
Null slicing in the v and u regions (boost H=±1). Inset magnifies the region
around the puncture x_p.

**What it is:** a landscape 3D-rendered schematic (not an axes plot — no
ticks, no frame) of the field magnitude as a **height-mapped, color-mapped
surface over the curvilinear (r, z²) grid**, with white lines tracing DG
element boundaries on the surface.

**Recorded attributes:**
- Plot type: 3D surface (`plot_surface`-like) with superimposed white
  element-edge lines; plus a detached small flat plane (the excised
  regularized/worldtube patch) drawn floating below/beside the main surface,
  outlined by a **red rectangle**, with a black dot labeled *x_p*.
- Colormap: terrain-like blue→teal→sand→brown ramp — deep blue (far field,
  flat) → teal/cyan → tan/orange → dark brown/maroon at the peak near the
  particle. Best matplotlib identifications: `gist_earth` / `terrain` /
  cmocean `topo`-family. NOT viridis (there is a distinct blue-to-brown
  hue reversal).
- Annotations: bold **t, u, v** with curved black arrows marking the three
  slicing regions along the top-left edge; rotated label "Kerr horizon r₊"
  along the left surface edge; double-headed arrow "z² ≡ cos(θ)²" and arrow
  "r" under the front-left corner; red text
  "Regularized field Ψ_m^R = Ψ_m − Ψ_m^P" pointing at the red-framed patch.
- Inset: bottom-right, thick black rectangular frame, showing a flat dark
  navy plane with **white grid lines at several nested refinement levels**
  (2:1 nested squares shrinking toward the puncture) — communicates hp-AMR.
- No axes, no colorbar, white background.

**Reproduction recipe (matplotlib):**
1. `fig = plt.figure(figsize=(~12, 6))`; main `Axes3D` with
   `ax.set_axis_off()`; camera azimuth roughly from lower-left, low
   elevation (~20–30°).
2. Surface: z = |Ψ_m| (strongly peaked at the particle), over an (r, z²)
   logical grid; `plot_surface(..., cmap='gist_earth', linewidth=0)`, then
   overplot white element-boundary curves (`ax.plot` white, lw≈1.5) at DG
   element edges; far region flat and blue.
3. Cut a rectangular hole at the worldtube; draw the regularized patch as a
   separate small, nearly-flat surface offset below, framed by red line
   segments; black dot + `$\boldsymbol{x}_p$`.
4. Add text/arrow annotations (`annotate` with arc arrows) for t/u/v,
   horizon, z², r; red bold text for the regularized-field formula.
5. Inset: second axes (`fig.add_axes`) with dark navy background and white
   nested grid lines, thick black spines.
- Match check: warped colored surface with white DG lines + flat blue far
  zone + detached red-framed patch + AMR inset. If any of these four
  structures is missing it is NOT a match.

---

## 2. Fig 2 — `matrix.pdf` (`fig:matrix`): DG operator sparsity + element numbering

**Caption (main.tex:625):** *Top:* matrix representation of the DG operator
at the coarsest grid level (equatorial symmetry). *Bottom:* numbering of the
DG elements corresponding to rows/columns. Each element has 4×4 grid points.
Jump conditions at the worldtube boundary (red rectangle) and between v, t,
u domains (red dotted lines).

**Recorded attributes:**
- Two vertically stacked sub-figures, single-column width.
- **Top: spy / block-sparsity plot**, 20×20 element blocks (each block
  16×16 entries), square aspect. Axis tick labels 1…20 on top edge (with
  header "Ω_k" at the corner) and left edge; small serif numerals.
  Diagonal blocks: dense diagonal-line texture with solid **black** squares
  (element self-coupling); off-diagonal neighbor couplings appear as light
  **gray hatched/checkered** squares scattered symmetric-ish about the
  diagonal. A short line segment with annotation "16 entries" below the
  matrix as a scale bar. **Red rectangles** frame the row band and column
  band of elements 3–6 (worldtube-adjacent) along the top and left tick
  rows.
- **Bottom: domain-decomposition sketch** (pure 2D schematic, thin black
  rectangles): horizontal r-axis marked r₊, r_v, r_p, r_u, …, R=10⁴; vertical
  z² axis with ticks 0, 0.25, 1. Twenty numbered boxes in two rows
  (bottom row z²∈[0,0.25]: elements 1,2,3,5,7,9,10 …; top row z²∈[0.25,1]),
  radially widening toward R=10⁴ (elements 9,10,19,20 much wider). Region
  labels v, t, u with brace/arrow spans on top. **Red solid rectangle**
  around worldtube elements (3,4,5,6); **red dotted vertical lines** at the
  v–t and t–u interfaces.
- Colors: black/gray/white + red accents only. Serif labels.

**Reproduction recipe:**
1. `fig, (ax_top, ax_bot) = plt.subplots(2, 1)` with the top axes square
   (`ax_top.set_aspect('equal')`), height ratio ≈ 2:1.
2. Top: build the 320×320 boolean sparsity of the coarsest DG matrix and
   `spy` it (or `pcolormesh` with a 2-level gray map); overlay black filled
   rectangles for dense sub-blocks; draw red `Rectangle` patches (fill=False)
   around rows/cols of the worldtube elements; ticks at block centers
   labeled 1..20, tick labels on top+left; add "16 entries" scale annotation.
3. Bottom: draw numbered rectangles with `Rectangle(fill=False)` + centered
   text numerals; annotate r-axis anchors (r₊, r_v, r_p, r_u, R=10⁴) and z²
   ticks (0, 0.25, 1); red rectangle + red dotted `axvline`-style segments.
- Match check: block-structured sparsity with dense diagonal + sparse
  neighbor blocks, 1..20 labels, red highlight bands, and a bottom sketch
  whose numbering is consistent with the matrix ordering.

---

## 3. Fig 3 — `precond.pdf` (`fig:precond`): preconditioner flowchart

**Caption (main.tex:632):** Preconditioning of the iterative linear solver:
each GMRES iteration preconditioned with one multigrid V-cycle built by AMR;
coarsest level solved directly (bottom solver); finer levels run 2 additive
Schwarz smoothing iterations solving element-centered subdomains in
parallel; subdomains solved directly (scalar) or by a nested GMRES
preconditioned by Laplacian approximations (gravity).

**Recorded attributes:**
- Pure flowchart (diagram-tool aesthetic, e.g. draw.io; sans-ish text in the
  boxes), full text width, white background.
- Box color code (pastel fills, thin dark outlines):
  - **Orange/peach** = iterative Krylov ("GMRES iteration", "Serial GMRES").
  - **Light blue** = smoothing/Schwarz ("Pre smoother: 2x Schwarz", "Post
    smoother: 2x Schwarz", "Additive Schwarz smoother", "Element-centered
    subdomain solver").
  - **Green** = direct solves ("Bottom solver: Direct LU (cached)",
    "Direct LU (cached)").
  - **Pink/red with vertical hatch segmentation** = "10 Laplacians" (drawn
    as a box subdivided into ~10 vertical slices).
- Left block: the V-cycle — GMRES box on top, arrows down through
  pre-smoothers to the green bottom solver, then up through post-smoothers;
  a **dotted rectangle** groups the finest-level pre/post smoothers with the
  side label "Finer levels added by AMR"; caption text "1x multigrid
  V-cycle" beneath. Repetition of GMRES iterations conveyed by
  **gradient-faded orange boxes** ("…") on either side of the GMRES box.
- Right block: zoom-out detail (connected by thin diverging lines from the
  smoother box): "Additive Schwarz smoother" over several parallel
  "Element-centered subdomain solver" boxes (one faded to imply more), each
  branching to alternatives labeled **Scalar** → green Direct LU, "OR",
  **Gravity** → orange Serial GMRES → pink "10 Laplacians" → green Direct
  LU (cached).
- Plain black straight arrows; up-and-down double arrows between smoother
  and subdomain boxes.

**Reproduction recipe:** simplest as matplotlib patches
(`FancyBboxPatch`/`Rectangle` + `annotate` arrows) or graphviz with fixed
positions; keep the 4-color semantic code (orange=Krylov, blue=Schwarz,
green=direct LU, pink=Laplacian approx), the dotted AMR grouping box, the
faded-edge boxes for "repetition", and the left(V-cycle)/right(Schwarz
detail) two-cluster layout joined by thin zoom lines.
- Match check: V-cycle down/up shape on the left, Schwarz detail with the
  Scalar-vs-Gravity OR-branch on the right, correct color semantics.

---

## 4. Fig 4 — `ssf_convergence.pdf` (`fig:convergence`): DG convergence

**Caption (main.tex:687):** Convergence of the DG method for modes up to
m=20 (a=0.5M, r0=10M). Exponential convergence despite the non-smooth
solution, under the hp-AMR scheme. Each mode solves in a few seconds on an
Apple M2 Pro.

**Recorded attributes:**
- Single-column, portrait; **two vertically stacked panels sharing x**,
  height ratio ≈ 3:1; x only labeled at the bottom:
  "Grid points per dimension (N_points)^{1/2}", linear, range ~20–80,
  ticks 20, 40, 60, 80.
- Top panel: y = "Absolute error |F_r^m − F_{r,ref}^m|", **log scale**
  10⁻¹² … 10⁻⁵ (every decade labeled). **21 solid lines, no markers**, in a
  grayscale ramp: black for m=0 → very light gray for m=20 (e.g.
  `plt.cm.gray(np.linspace(0, 0.8, 21))` or grays via
  `cmap='Greys_r'`-style mapping). Legend upper-right, framed, with only
  the two extreme entries: "m=0" (black), "m=20" (light gray).
  **Black dashed horizontal line at 1e-11** with inside-axes text label
  "Solver tolerance" just below it at the left.
  Visual structure: steep, roughly straight (≈exponential) descent on the
  semilog-y axes, lines crossing/jittering, flooring at the tolerance
  around N^{1/2} ≳ 60–70.
- Bottom panel: y = "Cumulative walltime", log scale, custom tick labels
  "1 s" and "10 s". Same 21-line grayscale family, tightly bundled, gently
  rising from ~0.2 s to ~5–6 s across the x-range (one lighter outlier
  starting higher at small N).
- Light gray grid both panels; CM serif fonts.

**Reproduction recipe:**
```python
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(5, 6.5),
                               gridspec_kw={'height_ratios': [3, 1]})
colors = plt.cm.gray(np.linspace(0.0, 0.75, 21))   # m=0 black ... m=20 light
# ax1: semilogy error vs sqrt(Npoints); ax1.axhline(1e-11, ls='--', c='k')
# ax1.text(..., 'Solver tolerance'); legend with only m=0 and m=20 proxies
# ax2: semilogy cumulative walltime; custom yticks [1,10] -> ['1 s','10 s']
```
- Match check: 21 grayscale lines, exponential (straight-on-semilog) drop
  over ~6–7 decades, dashed tolerance floor at 1e-11, walltime sub-panel of
  seconds-scale bundled curves.

---

## 5. Fig 5 — `ssf_modes.pdf` (`fig:modes`): stacked m-mode planes, 3 spins

**Caption (main.tex:694):** m-mode solutions for a scalar charge at the
ISCO for three spins. Each plane is an m-mode |Ψ_m| up to m=20, first few
modes highlighted in color; regularized field inside the worldtube; domain
z²∈[0,1], to R=10⁴M, clipped at r = 2 r0.

**Recorded attributes:**
- **Three 3D panels in a row** + a color-swatch legend column at far right.
  Panel titles (math serif, above each box):
  `a = −0.998, r0 = 9` | `a = 0, r0 = 6` | `a = 0.998, r0 = 2`.
- Each panel: a tall 3D box (thin black frame) containing **21
  semi-transparent wireframe surfaces stacked by amplitude**: each surface
  is log₁₀|Ψ_m| over (r, z²). Vertical axis ticks 0 … −5 drawn on BOTH the
  left and right vertical box edges; shared axis label "log₁₀|Ψ|" once at
  the far left of the whole figure. x-axis "r" below each box with
  per-panel ranges (≈[3,18] ticks 5,10,15; ≈[2,12] ticks 2..12; ≈[1.2,4]
  ticks 1.5..4).
- Color code (legend right, filled rounded rectangles):
  m=0 brick red, m=1 golden yellow, m=2 teal/green, m=3 medium blue,
  m=4 salmon/light red, m≥5 black/dark gray (fading lighter for higher m).
  (Reads like a custom categorical palette, roughly
  ['#b0413e', '#eec33f', '#39b48a', '#3d7dbd', '#e2695c', '#333333'].)
- Visual structure: m=0 plane on top nearly flat ≈10⁰; successive modes
  lower and increasingly peaked/oscillatory; colored ribbons dive sharply
  (downward spikes) inside the worldtube region near r0; the m≥5 gray
  stack forms a dense oscillating "curtain" toward small amplitudes;
  right-hand side of each surface clipped at r = 2 r0.
- White background; no grid; serif math typography.

**Reproduction recipe:**
1. `fig, axes = plt.subplots(1, 3, subplot_kw={'projection': '3d'})`,
   wide figure (~14×6); elevation ~15–20°, azimuth ~ −60°.
2. For each m: `plot_wireframe(r, z², log10|Ψ_m|, color=c_m, alpha≈0.5–0.7,
   rstride/cstride small)`; alpha decreasing with m for m≥5 grays.
3. Set zlim [−5.5, 0.5], z-ticks 0..−5; suppress panes
   (`ax.xaxis.pane.set_visible(False)` etc.), keep box edges; clip r at
   2 r0.
4. Legend: manual `Patch` handles for the 6 color classes, placed right of
   the last panel.
- Match check: three boxes, 21 stacked sheets per box with the 5 colored +
  gray-curtain ordering top-to-bottom, worldtube dips near r0, panel titles
  with the exact (a, r0) triples.

---

## 6. Fig 6 — `linsolv.pdf` (`fig:linsolv`): linear-solver convergence

**Caption (main.tex:701):** Convergence of the preconditioned linear
solver (m=2, a=0.5M, r0=10M). Top: GMRES iterations for consecutive AMR
levels, each converging in <10 preconditioned iterations; first AMR level
converges immediately (direct bottom solve). Bottom: walltime per
iteration; first iteration of each new AMR level is slower (builds/caches
subdomain matrices).

**Recorded attributes:**
- Single-column, **two stacked panels sharing x**, ratio ≈ 3:2.
  x = "Cumulative linear solver iteration", 0–32, ticks every 4.
- Top: y = "Absolute residual", log, ticks 10⁻¹, 10⁻⁴, 10⁻⁷, 10⁻¹⁰, 10⁻¹³
  (every **3** decades). Data: single black solid line style with filled
  circle markers, drawn as **~8 disconnected descending segments** — each
  AMR level restarts near 10⁻⁴…10⁻² and falls below tolerance in ≤7
  iterations; the very first segment drops from ~10⁰ to ~10⁻¹⁵ in one step.
  Black dashed horizontal "Solver tolerance" line near 10⁻¹²; text label
  inside axes below-left of the line. Framed legend upper right:
  "GMRES iterations" (line + circle marker sample).
- Bottom: y = "Walltime per iteration [s]", log, labeled tick 10⁻¹.
  Same black line+circles, sawtooth rising pattern (0.02–0.7 s), spikes on
  the first iteration after each AMR level. Right-aligned annotation above
  the panel: "Walltime:  0:00:06".
- Light gray grid; CM serif.

**Reproduction recipe:**
```python
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(5, 5.5),
                               gridspec_kw={'height_ratios': [3, 2]})
# per AMR level l: ax1.semilogy(iters_l, resid_l, 'ko-')  (separate calls
#   -> disconnected segments); ax1.axhline(tol, ls='--', c='k')
# ax1.legend(['GMRES iterations']); yticks every 3 decades
# ax2.semilogy(all_iters, walltime, 'ko-'); ax2.set_xticks(range(0,33,4))
# ax2 title-side text: fig walltime as H:MM:SS, right-aligned above axes
```
- Match check: disconnected restart segments (NOT one continuous curve),
  each dipping through the dashed tolerance line; sawtooth walltime panel.

---

## 7. Fig 7 — `ssf_error_budget.pdf` (`fig:error_budget`): error budget

**Caption (main.tex:708):** Error budget in the extracted self-force.
Upper panels: value of the self-force; lower panels: relative error. The
numerical-resolution floor lies below the error to the reference
[Warburton & Barack 2010], so accuracy is truncation-dominated (sum over
m-modes), not resolution-limited.

**Recorded attributes:**
- Full-width figure, **2×2 grid**: columns = r0 = ISCO (left), r0 = 10
  (right); rows = small top strip + large bottom panel (height ratio
  ≈ 1:4). Columns share y-ranges; x shared within each column.
- x = "Spin a" with ticks −0.998, −0.5, 0, 0.5, 0.998 (7 data points per
  series at a ≈ −0.998, −0.7, −0.5, 0, 0.5, 0.7, 0.998).
- Top strips: y = |F_r|, log (10⁻⁴…10⁻² labeled at left); one black solid
  line with circle markers per column. Left column rises steeply toward
  a→1 (ISCO shrinks); right column shallow/dipping near a≈0.5.
- Bottom panels: y = "Relative error |F_r − F_r^ref|/|F_r|", log
  10⁻⁷…10⁻³. Three black series with circle markers:
  - solid = "Error to reference"
  - dashed = "Truncation at m ≤ 20"
  - dotted = "Numerical resolution"
  with **layered gray shading** (fill_between from each lower curve down),
  light gray under the dashed curve, darker gray under the dotted curve —
  the stacked-fill "budget" look. Ordering everywhere:
  solid > dashed > dotted.
- Panel tags "r0 = ISCO" / "r0 = 10" as plain text top-left inside the
  bottom panels. Framed legend in the RIGHT bottom panel, upper right.
- Peak of all error curves near a = 0.5; grayscale only; light grid; CM
  serif.

**Reproduction recipe:**
```python
fig, axs = plt.subplots(2, 2, sharex='col', figsize=(11, 4.8),
                        gridspec_kw={'height_ratios': [1, 4]})
# axs[0,i]: semilogy |Fr| vs a, 'ko-'
# axs[1,i]: three series 'ko-', 'ko--', 'ko:' + fill_between(a, curve, ymin,
#   color='0.9' / '0.75', zorder below lines)
# x ticks: [-0.998, -0.5, 0, 0.5, 0.998] with those literal labels
# text tags r0=..., legend only in axs[1,1]
```
- Match check: 1:4 row-height split, three-tier solid/dashed/dotted
  hierarchy with gray budget fills, truncation curve hugging the reference
  error and resolution floor well below both.

---

## 8. Checklist before claiming any reproduction "matches"

1. Same panel count, arrangement, and height ratios.
2. Same axis quantities, scales (log/linear), ranges, and tick labeling
   style (incl. unit ticks like "1 s").
3. Same series encoding: grayscale ramps / line-style hierarchy / marker
   usage / categorical palette of Fig 5 / red-accent semantics.
4. Same qualitative data structure (exponential slopes, restart segments,
   budget ordering, stacked-plane ordering, worldtube dips).
5. Legend and annotation placement equivalent (framed, upper-right;
   inside-axes tolerance labels).
Only then, and only with a recorded VLM comparison verdict, may a figure be
declared reproduced (see `.claude/skills/vlm-figure-repro/SKILL.md`).
