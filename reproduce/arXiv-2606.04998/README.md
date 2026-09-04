# arXiv:2606.04998 — Vu, Nishimura, Osburn, Thompson, Kidder, Upton & Wardell, "Self-force calculations with numerical relativity methods"

Mirrors (read-only, emri-disk repo): tex `ref-paper/arxiv-2606.04998/src/main.tex`, the
production input file `ref-paper/arxiv-2606.04998/src/anc/ScalarSelfForce.yaml`, the
reference figure PDFs (`ssf_domain.pdf matrix.pdf precond.pdf ssf_convergence.pdf
ssf_modes.pdf linsolv.pdf ssf_error_budget.pdf`) in the same directory, and the spectre
`SolveScalarSelfForce` executable under `ref-code/spectre` @ d3cf5af0.

Physics lives in `core/EllipticCore` (the complex-parametric 2D hp-DG elliptic solver stack:
mesh, mortars, matrix-free operator, boundary-data hooks, Schwarz, multigrid, GMRES) and
`capabilities/ScalarSelfForce` (the m-mode Kerr scalar self-force: orbits, the flux-form
operator, the vtu null slicing, the Detweiler–Whiting puncture and effective source, the
block lattice, the worldtube formulation, the jump/boundary tables and the operator
assembly). **Drivers here are THIN** — zero physics, per `reproduce/README.md`.

Style manual for every figure of this paper: `reproduce/_style/manuals/arXiv-2606.04998.md`
(a copy; the canonical, evolving version is emri-disk
`progress/elliptic-ssf/FIGURE_STYLE_MANUAL.md`). The recurring VLM procedure is mandatory —
read the manual (`reproduce/_style/manuals/arXiv-2606.04998.md`), render, VLM-check the render
alone, VLM-compare against the reference PDF, decide script-vs-structure, iterate.

> **Parameters.** The paper's figure captions name no parameters. Every driver here uses the
> production example of the paper's own input file, `anc/ScalarSelfForce.yaml:8-12`:
> `spin = 0`, `orbital_radius = 10 M`, `m_mode = 2`, `inner_radius = 2 M` (= r₊ at a = 0),
> `outer_radius = 10⁴ M`, `worldtube_radius = 4.8 M`, `worldtube_dz2 = 0.25`, `P = 4`,
> `L = 0`, with `ImposeEquatorialSymmetry: True`. That is an ASSUMPTION where the figure does
> not pin the value, and each driver states which of its registered numbers depend on it. Where a
> caption DOES name a parameter, the caption wins: Fig. 6 says "the m = 2 mode of the a = 0.5 M,
> r0 = 10 M configuration" (main.tex:701), so `fig_linsolv.jl` runs at `spin = 0.5` and takes
> `inner_radius = r₊(a = 0.5)` accordingly.

> **Why `InitialCoarseLevels: L` means "no coarser levels" at `L = 0`.** Primary source:
> spectre `src/ParallelAlgorithms/LinearSolver/Multigrid/Tags.hpp:30-37` — *"Maximum number of
> coarser levels in the initial multigrid hierarchy. Set to `'0'` to create no coarser levels
> initially. Set to `'Auto'` to coarsen all the way up to single-element blocks."*

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
| 1 | `fig:domain` — the DG computational domain: `\|Ψ_m\|` height- and colour-mapped over the (r, z²) element grid, `\|Ψ^ℛ_m\|` in the worldtube | `fig_domain.jl` | `expected/fig_domain.json` | heavy | **solid** | **match** | A single 3D panel plus an inset: |Ψ_m| height- **and** colour-mapped over the (r, z²) DG element grid with white element edges and the Legendre-Gauss collocation grid, |Ψ^ℛ_m| inside the red-framed worldtube, the v/t/u slicing arc with H = ∓1, and the six-step hp-AMR cascade at x_p magnified bottom-right. Seven AMR levels, 20 → 60 elements and 320 → 3762 dofs, GMRES 1, 7, 6, 5, 5, 4, 4 (the counts and the level map moved with the all-dimension 2:1 closure and the coarse-mesh chain; see the note below). Structure reproduced: the two |Ψ_m| lobes of the published render — the tan ridge in the polar row over the worldtube's radial span and the sharper wedge just outside the tube in r at the equator — agree in count and position, as do the flat navy near-horizon band, the low navy tube patch and the teal wavezone plateau. `--check` 52/52 at rtol 1e-10 with max|Δ| = 0, `-t1` ≡ `-t4`. ~500 s at -t4, of which the sweep is 223.5 s — **heavy**. |
| 2 | `fig:matrix` — DG operator sparsity at the coarsest grid level + element numbering | `fig_matrix.jl` | `expected/fig_matrix.json` | heavy | **solid** | **match** | Two panels: the B/W spy of the 320×320 coarsest-level matrix over a light-gray element-block checkerboard, and the domain sketch numbering the 20 elements that index its rows and columns. Structure reproduced exactly: 20 diagonal blocks + 33 symmetric off-diagonal block pairs, of which the 8 pairs across the four 2:1 interfaces (2↔3, 2↔4, 7↔9, 8↔9, 12↔13, 12↔14, 17↔19, 18↔19 in the paper's numbering) are fully dense 16×16 and the other 25 carry 64 entries. Adjacency graph-isomorphic to the domain topology in BOTH directions (0 matrix-only, 0 topology-only pairs). Reference block-occupancy extracted from `matrix.pdf` at 600 dpi and compared block-by-block: identical. 398.8 s at -t4 — **heavy**, although the numerics are ~10 s: the cost is the first-call specialization of the ComplexF64 `MeshOperator` (104 s) and of the first Plots/GR render. |
| 4 | `fig:convergence` — the DG convergence of every `F_r^m` up to m = 20 over the hp-AMR sequence, + the cumulative walltime | `fig_convergence.jl` | `expected/fig_convergence.json` | heavy | **preliminary** | **match** | Two panels for all 21 modes at a = 0.5 M, r₀ = 10 M: `\|F_r^m − F^m_{r,ref}\|` against `sqrt(N_points)` with each mode's own highest-resolution level as its reference, and the cumulative walltime below. Presentation reproduces item by item (2 panels 3:1 sharing x, 21 grayscale lines at the published ramp `m/23` with no markers, log y over ~8 decades, the framed upper-right `m=0`/`m=20` legend, the dashed "Solver tolerance" line with its inside label, "1 s"/"10 s" ticks, light-gray grid, CM serif) and every mode converges exponentially. The **abscissa reproduces exactly**: the driver runs two lattices from one pool and plots the published one — full-z, `ImposeEquatorialSymmetry: False`, as spectre's own `ScalarCircularOrbit.yaml` regression sets it and NOT as the shipped `anc/ScalarSelfForce.yaml:36` does — whose grid-point sequence `480, 975, 1604, 2427, 3468, 4751, 6300` equals the published vertices level for level (GATE 5, all seven). The half-domain lattice (`true`, 320 → 3762) is kept as a registered and gated secondary series, never drawn. Against the published envelope: rates [−0.1395, −0.1020] (published [−0.1421, −0.1043]), R² ∈ [0.9143, 0.9932] median 0.9781 with 10/21 ≥ 0.98 (published 0.7996 … 0.9929, median 0.9833, 9/21), 6.02 … 8.38 decades dropped (published 6.06 … 8.12), last plotted error 5.6e-14 … 7.6e-13 (published 8.9e-14 … 7.9e-13). One STRUCTURAL deviation remains, escalated: the published grids **depend on m** (975 for m ≤ 4 against 970 for m ≥ 5) while ours are gated m-independent, which the shipped criteria (`RefineAtPuncture` ∘ `IncreaseResolution`) cannot produce — now a 0.5 % effect rather than the 50–78 % lattice mismatch it was hiding behind. See the Fig. 4 section for the evidence and the `TruncationError` hypothesis. `--check` 86/86; 760 s at -t4 on a 5-worker pool — **heavy**. |
| 6 | `fig:linsolv` — GMRES residual staircase over AMR grids + walltime per iteration | `fig_linsolv.jl` | `expected/fig_linsolv.json` | heavy | **solid** | **partial** | Two panels for the m = 2 mode at a = 0.5 M, r0 = 10 M on the **full-z** lattice (the published one; the yaml says otherwise — see below). Eight AMR grids, 30 → 119 elements, 480 → 8139 dofs. The grid sequence reproduces the paper's own convergence-figure abscissa — 480, 975, 1604, 2427, 3468, 4751, 6300 — **exactly, all seven, zero deviation**. Iterations `1, 7, 6, 5, 4, 3, 3, 3` (32 cumulative) against the published `1, 7, 6, 5, 4, 3, 3, 2` (31), **seven of eight grids exact**; ‖b‖₀ 1.39146 against 1.39401 (−0.18 %); per-iteration reduction 28.3, 44.5, 51.3, 82.6, 260.4, 192.2, 234.3 against 26.0, 42.7, 48.0, 82.3, 244.6, 191.9, 567.0. All eight gates pass. Remaining holes: the last grid takes 3 iterations against 2, and the walltime panel's amplitude (eager vs lazy build). |
| 5 | `fig:modes` — 21 stacked `log₁₀\|Ψ_m\|` planes over the (r, z²) domain, at three spins | `fig_modes.jl` | `expected/fig_modes.json` | heavy | **preliminary** | **match** (3 of 3 panels) | Three 3D panels — `a = −0.998` and `a = 0` at their ISCOs (8.9944, 6 M) and `a = +0.998` at the **published** `r₀ = 2 M` — each carrying 21 semi-transparent wireframe planes of `log₁₀\|Ψ_m\|` over (r, z²), `Ψ^ℛ_m` inside the worldtube and `Ψ_m` outside, clipped at `r = 2 r₀`, m = 0…4 in the paper's categorical colours and m ≥ 5 in the published `m/23` grey ramp. NOT a mode-sum curve: the caption plots the FIELD. All three match item by item, including the tick sets (5/10/15, 2…12, 1.5…4.0), the two m = 0 dives per panel and the grey curtain's cliff at `r₀ + w_r`. The horizon order hold is chosen by a MEASURED rule — sweep with the package default, read `horizon_axis_profile(...).layer_ratio` for every mode, re-sweep held only if the maximum over m exceeds 1.05 — which leaves panels 1-2 on the default (0.1105, 0.3108) and holds panel 3 (3.889 → **0.5776**); GATE 7 gates the plotted value and the panel states both numbers and the setting. **preliminary** for one measured reason, flagged on the figure: at `a = 0.998, r₀ = 2` GMRES stagnates at m = 19 and 20 (final residual 9.3e-6 against the 1e-12 tolerance), so those two planes are drawn dashed and named, `pK_converged` carries the verdict per mode, and panel 3's decay claim breaks only at 19 → 20. 66 elements and 4814 dofs (4604 held) at the 8th grid from 320 at level 0. Four sweeps (three panels + panel 3's held re-sweep) through ONE `mode_pool` of 5 single-threaded workers: 1286 s — **heavy**. |
| 7 | `fig:error_budget` — the F_r error budget: value, and resolution / truncation / reference error against Warburton–Barack, over spin at the ISCO and at r₀ = 10 M | `fig_error_budget.jl` | `expected/fig_error_budget.json` | heavy | **preliminary** | **partial** | 2×2: a `\|F_r\|` strip over a budget panel, per column (r₀ = ISCO, r₀ = 10 M), nine spins −0.998 … +0.998 — the published sample set, read off the PDF's own vector paths, not the five tick labels. Three grayscale curves with the layered gray budget fills: solid "Error to reference" (only where WB print a cell — one on the ISCO column, seven on r₀ = 10; the published solid curve stops at \|a\| = 0.9 for the same reason), dashed "Truncation at m ≤ 20", dotted "Numerical resolution". What is plotted is the paper's own quantity, the BARE partial sum with no fitted tail (three independent confirmations in the section below); the Heffernan-tail F_r is registered beside it. MEASURED: the published dashed curve is reproduced to **0.3–3.6 %** at 17 of 18 configurations, the published solid curve to **within 23 %** at all 8 printed WB points (ratios 0.775–0.988), the main.tex:726 offset lands at 2.28–3.75 against the published 2.81–4.63, and the unregularized F_t matches the printed Table II to ≤ 1e−6 at all 8. **main.tex:723's ordering is a LATTICE property**: resolution below truncation at 15/18 on the paper's full-z lattice against 0/9 on the shipped half domain at the same 8 AMR grids — the controlled variation that closes the row-46 finding. **preliminary** for two measured reasons, both ring-marked or named: at a = +0.998 at its ISCO the mode sum has not converged at m = 20 (last term 2.4 % of the total), and at \|a\| = 0.998, r₀ = 10 the resolution error is 74–100× the published one with the layer_ratio instrument reading only 0.098. 27 sweeps (18 + 9 secondary, plus one horizon-hold re-sweep) through ONE `mode_pool` of 6 single-threaded workers: 4342 s — **heavy**. |

### What "the coarsest grid level" is

`LinearSolver.Multigrid.InitialCoarseLevels: L` (yaml:188) counts the COARSER levels only, so
the production hierarchy has `L + 1` levels; at the paper's `L = 0` the coarsest level IS the
initial mesh — **20 elements × 4×4 grid points = 320 dofs**, which is exactly what the
reference figure's tick labels 1…20 show. The 5×2 block lattice carries the yaml's
`RefinedLevels` table (yaml:78-110): the source column pre-refined to `(L+1, L+1)`, the
`r₀+w`…`r_u` column split angularly, the `r_u`…`R` column split radially. Coarsening all the
way to the irreducible block structure instead would give 10 elements / 160 dofs; that number
is registered as `nelements_block_structure` / `ndofs_block_structure` so the two readings
stay distinguishable.

### Element numbering — a convention difference, gated

The paper numbers elements by block (the yaml flat block order) and then, within a block, by
the RADIAL segment index major and the ANGULAR index minor. GRTools' `Mesh2D` enumerates the
transpose, so its elements 4, 5 (and 14, 15) are the paper's 5, 4 (15, 14). Both panels are
drawn in the PAPER's order — a symmetric relabeling `P·A·Pᵀ` of the same operator, computed
from the mesh's own `ElementId`s by the sort key `(block, radial segment, angular segment)`.
The permutation is registered as `element_mesh_index`; every gate is permutation-invariant
and is asserted on the mesh-order data.

### The sparsity threshold

Exact-nonzero counting gives 10112 entries, 576 of which are roundoff: on the coarse side of
a 2:1 interface the mortar projections satisfy `Σ R·P = I` analytically, so the 144
off-pattern entries of each of those four diagonal blocks cancel. Measured there, the largest
off-pattern entry is 1.3e-17 … 9.6e-17 of its own block maximum against a smallest in-pattern
entry of 4e-3 … 7e-3 — a gap of 10¹³ … 10¹⁴.

The structural count is **9536** for every relative threshold from **1e-7 to 1e-16 — ten
decades measured** (9528 at 3e-7, 9480 at 1e-6, 9838 at 1e-18 where roundoff enters); the
driver asserts the inner window 1e-8 … 1e-16 at five thresholds as `nnz_threshold_plateau`
and uses `1e-12·max|A|` for the figure and the block table. The plateau stops around 1e-7 on
the high side because `max|A|` is set by the horizon column while `A^rr` varies over decades
across the compactified domain, so the outer wavezone blocks' genuine structural entries are
small in *global* relative terms.

For the comparison with the published figure this matters in one way only: the reference does
**not** draw those four coarse-side diagonal blocks dense — it draws them with the same
112-entry tensor pattern as every other diagonal block — which is consistent with the 9536
structural entries of our own count, and is what `compare_blocks.py` measures block by block.

## Fig. 1 — six h-steps, the drawn radial coordinate, and what Plots cannot do

### Why `iterations = 6`

main.tex:558, about this figure specifically: *"The result of this h refinement strategy after
six refinement steps is shown in \cref{fig:domain}, with the region around the puncture magnified
in the inset."* A *refinement step* is one `amr_step!` — spectre's `Amr.Iterations` counts
**solves**, and its own trigger (`EveryNIterations {N: 1, Offset: 1}`, yaml:157-160) skips
refinement on the initial guess, so `n` solves are separated by `n − 1` refinements — hence
`solve_mode_amr(...; iterations = 6)`: seven levels, 0…6, not the yaml's production
`Amr.Iterations: 8`. The published inset confirms it twice over: its cascade shows six
2:1 steps into the puncture corner, and its coarse elements carry a 9-cell collocation grid, i.e.
`N = P + 6 = 10`.

### The abscissa is the paper's own visualisation coordinate

A domain running to `R = 10⁴ M` cannot be drawn against linear `r`. The paper's run observes
`RadiallyCompressedCoordinates: {InnerRadius: *t_to_u_transition, OuterRadius: *outer_radius,
Compression: *outer_shell_distribution}` (yaml:222-225), and spectre implements it
(`src/Domain/RadiallyCompressedCoordinates.cpp:27-52` with `src/Domain/CoordinateMaps/Interval.cpp:120-129`)
as the inverse of an `Interval` map from `[r_in, r_in·log₁₀R]` onto `[r_in, R]` — at the
production point

    r_c(r) = 97.77904 − 1790.4487/r  on  [24.4, 10⁴] → [24.4, 97.6],    r_c(r) = r  below r_in.

That map and its return trip are **`radially_compressed` / `radially_compressed_inv`** in
`capabilities/ScalarSelfForce/src/Domain/config.jl`, taken from `cfg.observers.rcc_*` and
gated in `test_runtime_config.jl` against an independent transcription of the two C++
expressions (`Interval.cpp:123-131` for the map, `:81-89` for the inverse, plus the
`compressed_outer_radius = r_in·log₁₀R` of `RadiallyCompressedCoordinates.cpp:34-39`). The
driver only *calls* them: the 66 golden samples that are laid out uniformly in the drawn
coordinate are pinned against the package, never against arithmetic the driver owns. **Framing:** the published render itself stops at the u-block element boundary
`r = 48.68 M`; measured on `ssf_domain.pdf` the published worldtube patch sits at 4.5–21.5 % of
the drawn ribbon, which the compressed coordinate predicts as 5.6–21.9 % if the ribbon is cut
there (against 3.5–13.5 % for the full `r_c ∈ [1.87, 97.6]`). The driver frames the same way —
`R_FRAME` is *derived* as 1.068 × the u block's own first element boundary
(`map_forward(blockmap(u_block, 1), 0)` = 48.681 M), so the outermost wavezone element runs
off the right edge exactly as it does in the paper. Nothing is removed from the solve or from the golden.

### What is plotted

`sol.u` everywhere, which by the row-42 `ModeSolution` contract is `Ψ^ℛ_m` inside the worldtube
and `Ψ_m` outside it — verbatim what the caption asks for. `reconstruct_full_field` (`Ψ_m`
everywhere, log-divergent at `x_p`) is computed only for the gates and is never drawn.

### Rendering — an honest list of Plots.jl-vs-ParaView differences

The reference is a ParaView render annotated in Keynote. Plots.jl/GR has no 3D-surface primitive
that carries a per-quad colour, a hard field discontinuity at the worldtube boundary and white
element edges at once, so the driver does its own axonometric perspective projection and paints
depth-sorted quadrilaterals into one `:shape` series (`fill_z` carrying a (value, shade)
categorical palette). Consequences, none structural:

* the golden PNG is held to the shared **200 KB** reproduce budget by one knob, `PX_PER_UNIT`
  (the raster resolution of the frame). Point sizes are absolute while the annotation offsets
  are in data coordinates, so `SC = PX_PER_UNIT/429` scales type and stroke with it and the
  layout is invariant under the knob. The canvas is derived from the frame, so
  `aspect_ratio = :equal` letterboxes nothing. Measured at 272 px/unit: **191 839 B**. Palette
  reduction is *not* the lever here — cutting the (value, shade) palette from 768 entries to 192
  moved the file by 1 % (326 666 → 323 626 B), because GR antialiases every one of the 4457 quad
  seams and every grid line, so the byte count follows the resolution, not the colour count;
* `--check` renders to a **scratch path** and leaves the tracked PNG alone, so a check can never
  rewrite a golden artifact; the size verdict is folded into the gate list *before* `gates_pass`
  and `conclude!`, so an oversized figure is a FAIL row with a full comparison table rather than
  a stacktrace. The vector PDF is written only on `--regen`/`--plot-only` and is kept out of the
  index by `expected/.gitignore` beside the goldens — the root `.gitignore` is untouched;
* `--plot-only` redraws from a **surface cache** (`expected/fig_domain_surface.jls`, untracked,
  53 888 B) that `--regen` writes: 57 element sample lattices, ~8k values, which do not belong in
  a 200 KB golden. It turns a cosmetic iteration from a ~500 s solve into a ~40 s redraw. Absence
  of the cache is a HOLE, never a silent fallback;
* element edges and collocation-grid lines are drawn after all quads, so a line on a far element
  is not hidden by a nearer fold (sub-pixel on this surface);
* lighting is one Lambert term over an ambient floor — the reference's ParaView shading is what
  makes its high-|Ψ| slopes read dark brown rather than tan;
* the camera, the 2:1:0.78 scene box, the colour anchors (sampled from the reference render) and
  every annotation offset are style constants, listed together at the top of the figure section
  of the driver.

### Gates

1. **The refinement/order map**, asserted level by level against the sealed AMR metadata rule:
   the criterion flags exactly the two elements containing `x_p` (the `z²` chart puts it on the
   domain boundary), they sit at h-level `L+1+k` in both dimensions with the order **held** at
   `(P, P) = (4, 4)`, every other element descends from exactly one level-`k−1` element and
   carries that ancestor's order **+1**, and every level is 2:1 balanced.
2. **The field**: finite and non-zero at all 130 probes including the whole horizon column
   and the poles; `puncture_chart_position` equal to `(r₀, 0)`; and `|Ψ^ℛ_m| < |Ψ_m|` at every
   inner worldtube-boundary probe (worst ratio 0.0469) with the regularized field 17×
   flatter across the tube than the retarded field.
3. **Determinism**: `--check` after `--regen` at `rtol = 1e-10`, `-t1` and `-t4` identical.
   A sixth gate holds the tracked PNG inside the 200 KB budget and is registered as
   `png_bytes_within_budget`.
4. **VLM** against `ssf_domain.pdf`: **match**. All four structures the manual demands (warped coloured surface with white DG lines, flat far zone, detached red-framed patch, AMR inset) and 9/9 checklist items; the lobe count and placement agree. Five style differences are recorded in `evidence/fig-domain/vlm_verdict.md` — proportions (the published ParaView camera is not recoverable from the PDF), ParaView-vs-Lambert lighting, the added `H = ±1` labels, the drawn cliff face at the tube boundary, and the inset window.

### The mesh sequence is anchored to the all-dimension 2:1 closure

Regenerated at `13ea24c` (`amr — the 2:1 closure in all dims`), which balances the refinement
in **every** dimension rather than only the normal one. The half-domain sequence gains the
tangential closures and becomes

    elements  20, 30, 36, 42, 48, 54, 60          (was 20, 27, 33, 39, 45, 51, 57)
    dofs     320, 632, 1014, 1506, 2120, 2868, 3762

so the per-element level map, `element_id_sha256` and `level_map_sha256` all move with it,
and the field probes follow. What does **not** move: the GMRES counts (`1, 7, 6, 5, 5, 4, 4`),
the sample coordinates, the puncture h-levels and held orders, the order window and every gate
flag. GATE 1b still holds and now exercises the anisotropic case it was written for — two
elements carry the anisotropic order pair `(9, 10)`, i.e. the closure held the order in the
dimension it split and incremented the other, and the level map now contains the `(0, 1)` and
`(1, 0)` h-level pairs of single-dimension balance splits. In the inset the extra splits show
up as a graded band wrapping the cascade, which is closer to the published inset than the
previous ungraded step was.

### The order map is anchored to the fixed `hp_refine`

EllipticCore closes the caller's AMR flags under the 2:1 rule (`balanced_flags`) *before* the
p half, so an element the balancer splits keeps its polynomial order in the split dimension —
spectre `Domain/Amr/UpdateAmrDecision.cpp:82-87` with
`ParallelAlgorithms/Amr/Projectors/Mesh.cpp:30-34`. This golden was generated after that fix
(`6df9f51`) and GATE 1b asserts the closed rule directly: measured, the two balance-split
children of block 1 at level 1 carry `(5, 4)`, not `(5, 5)`. `element_id_sha256`, the drawn
element boundaries and the inset cascade are independent of the closure (the h-refined element
**set** is the same either way); `element_order_x/y`, `level_ndofs`, `level_map_sha256`,
`level_gmres_iterations` and the sampled field values are not, and would move against a tree
without it.

## Fig. 6 — the lattice, the resolved staircase, and what still does not reproduce

### The grid count is a convention mapping, derived from the source

`Amr.Iterations: 8` (yaml:153) counts spectre's AMR **solver iterations**, i.e. *solves*, and the
yaml's own trigger `EveryNIterations {N: 1, Offset: 1}` (yaml:157-160, commented "Run AMR in every
iteration, but not on the initial guess") skips refinement on iteration 0. Tracing spectre's loop
— `src/Elliptic/Executables/Solver.hpp:348-357`
(`RunEventsAndTriggers → StopAmr → ExecutePhaseChange → StepActions(solve) → IncrementIterationId`)
against `src/Elliptic/Amr/Actions.hpp:56` (`HasConverged{num_iterations, 0}`) — the sweep is
**8 solves separated by 7 refinements**. GRTools' `solve_mode_amr(; iterations = n)` counts
`amr_step!`s, so the paper's configuration is `iterations = 7`; the driver derives it as
`cfg.amr.iterations - cfg.amr.trigger_offset` rather than hard-coding a number. The published
figure confirms the reading independently: its top panel carries exactly **8** disconnected
staircases on a 0…32 axis.

`ESSF::amr-driver`'s own default is `iterations = 8`, i.e. **nine** grids — one refinement more
than the paper. The secondary `prod_*` series registered here is exactly that nine-grid sweep on
the yaml's half domain.

### Which lattice the published figure is on — the one deviation from the yaml

The plotted sweep runs on the **full-z** lattice, `equatorial_symmetry = false`: 30 elements and
480 grid points at the initial level, with the puncture contained in **four** elements. The paper's
own input file says the opposite — `ImposeEquatorialSymmetry: True` (yaml:36) — which folds the
domain onto `z² ∈ [0, 1]`. The yaml's `RefinedLevels` table (yaml:78-110) lists ten block entries,
five "middle" rows (`LowerCornerIndex: [·, 0]`) and five "top" rows (`[·, 1]`); with the angular
block bounds `[0., 0.25, 1.]` and `P = 4` that is 20 elements × 16 points = **320** points under
equatorial symmetry, and 30 elements × 16 = **480** without it.

Three measurements say the FIGURE was made on the full lattice, whatever the yaml says:

| quantity | published | full-z (plotted) | half domain (yaml) |
|---|---|---|---|
| ‖b‖₀ = grid 0's plotted initial residual | 1.39401 | **1.39146** (−0.18 %) | 1.7444 (+25.1 %) |
| grid points per AMR grid | 480, 975, 1604, 2427, 3468, 4751, 6300 | **exactly those, all seven** | 320, 632, 1014, 1506, 2120, 2868, 3762 — none of them |
| iterations per grid | 1, 7, 6, 5, 4, 3, 3, 2 (31) | **1, 7, 6, 5, 4, 3, 3, 3 (32)** | 1, 7, 6, 5, 5, 4, 4, 4, 3 (39 over nine grids) |

‖b‖₀ is decisive because the initial guess is zero, so grid 0's plotted initial residual **is** the
source norm. The half-domain sweep is registered too, as the secondary `prod_*` series, so the
yaml-versus-figure discrepancy is data in the golden rather than prose; the figure itself carries a
one-line note saying which lattice it is.

### Why the staircases shorten — the resolved question

This node previously plotted a **flat** `1, 8, 9, 9, …` staircase against the published shortening
one, and the gap was owned as an open finding. It is closed, and the cause was the **multigrid
hierarchy**, not the figure. spectre does not build the hierarchy by re-coarsening the current
grid: it uses the **AMR history** — every grid the sweep has already visited becomes a coarser
multigrid level of the next
(`Elliptic/Executables/SolveScalarSelfForce.hpp:173`;
`ParallelAlgorithms/Amr/Actions/AdjustDomain.hpp:159-163, 203-212`;
`SendDataToChildren.hpp:43-70`; `LinearSolver/Multigrid/ElementActions.hpp:96-135`;
`Multigrid.hpp:56-58`). With `hierarchy = :history` (now the `solve_mode_amr` default,
`fdcd919`) the preconditioner improves with depth exactly as the paper's does.

| grid | iters (paper) | iters (ours) | ‖r₀‖ paper / ours | reduction·iter⁻¹ paper / ours |
|---|---|---|---|---|
| 0 | 1 | 1 | 1.39401 / 1.39146 | — |
| 1 | 7 | 7 | 2.6292e-3 / 2.5725e-3 | 26.0 / 28.3 |
| 2 | 6 | 6 | 7.1549e-4 / 6.6525e-4 | 42.7 / 44.5 |
| 3 | 5 | 5 | 8.8591e-5 / 8.6975e-5 | 48.0 / 51.3 |
| 4 | 4 | 4 | 2.3976e-5 / 2.1682e-5 | 82.3 / 82.6 |
| 5 | 3 | 3 | 4.0581e-6 / 3.6897e-6 | 244.6 / 260.4 |
| 6 | 3 | 3 | 9.7638e-7 / 8.9674e-7 | 191.9 / 192.2 |
| 7 | **2** | **3** | 2.6296e-7 / 2.2148e-7 | 567.0 / 234.3 |
| **total** | **31** | **32** | | |

The paper column is extracted from `linsolv.pdf`'s own **vector paths** (PyMuPDF), not from
pixels. The top panel's y axis is calibrated on the **exact** flipped-y anchors 20.44117 and
103.93375 — the 10⁻¹ and 10⁻¹³ gridlines. A first extraction rounded those to 20.4 / 103.9 and put
the whole published residual series ~1.35 % low (‖b‖₀ 1.375 instead of 1.39401); the iteration
counts and the reduction-per-iteration ratios are invariant under that rescaling, ‖b‖₀ is not, so
the ‖b‖₀ gate band is set from the corrected value. "Reduction per iteration" is
`(‖r₀‖ / ‖r_final‖)^(1/iters)`. An independent re-extraction from the PDF content stream agrees.

### Remaining holes

1. **The last grid takes 3 iterations against the published 2** (32 cumulative against 31). Its
   reduction factor is 234 against the published 567 — the finest grid is where our preconditioner
   is still weakest. Gated as a regression band (`total_iterations ≤ 33`), not papered over.
2. **The walltime panel's amplitude**, because this package builds the Schwarz subdomain matrices
   eagerly — see the next section. The panel's *shape* reproduces; its spike/plateau ratio and its
   y-range do not.
3. ~~The published grid sequence is m-dependent and ours is not.~~ **Closed by `13ea24c`** (the
   2:1 closure applied in all dimensions). The grid sequence is now
   `480, 975, 1604, 2427, 3468, 4751, 6300` — the paper's own convergence-figure abscissa, **all
   seven exact, zero deviation** — and gate 8 asserts it. Our eighth grid, 8139, is the one the
   paper holds as its reference and so has no published counterpart.

### The bottom panel: eager build, attributed by addition

The paper's first-iteration bump is the **lazy** construction of the Schwarz subdomain matrices
inside that iteration (main.tex:702). This package builds them **eagerly** in `rebuild!`, so the
same work is measured in `SolveRecord.build_time` and *not* inside `iteration_time[1]`. The bottom
panel therefore plots `iteration_time[1] + build_time` at each grid's first point: the same
physical work charged to the same iteration, by addition rather than by measurement inside the
timer — a bookkeeping equivalence, not a lazy build. Both columns stay separately registered in
the golden, and the figure itself carries the note *"first point of each level includes the eager
operator build"* under the bottom panel, so the panel cannot be read as a lazy-build measurement.

Measured on the landed full-z golden: the build runs 0.85 s → 24.9 s across the seven refinements,
which is why the panel spans three decades where the published one spans about one and a half, and
why the walltime stamp reads 0:01:04 against the published 0:00:06. Without the attribution the
first/rest ratio is flat, i.e. there is no first-iteration feature here at all.

### Timings are reported, not gated

Wall-clock values are registered under `reported_*` names with a deliberately loose tolerance, so
`--check` cannot turn a busy machine into a red result. What is gated is the structure they carry
(`first_bar_is_slowest_per_level`) and the residual histories, at `rtol = 1e-8` — three orders
looser than the measured `-t1` vs `-t4` spread, which is **bitwise zero** with
`OPENBLAS_NUM_THREADS=1`, and still an order tighter than the 4.2e-12 relative drift an unpinned
BLAS thread count introduces. `residual_log10_flat` is registered beside `residual_flat` because
the harness compares arrays norm-wise: on the raw vector the O(1) initial residual would mask a
factor-2 move at the 1e-13 floor.

## Fig. 4 — the paper's lattice is not the shipped yaml's, and the criteria are not either

`fig_convergence.jl` plots what main.tex:719-724 says it plots: for each of the 21 modes, the
absolute error `|F_r^m − F^m_{r,ref}|` at every AMR level against **that mode's own
highest-resolution level**, over `sqrt(N_points)`, with the cumulative walltime beneath. The
reference level is at zero error by construction and is therefore not drawn, so an 8-grid sweep
gives **7 plotted points** — exactly the vertex count of the published curves. The level count
is the same convention mapping `fig_linsolv.jl` derives (`Amr.Iterations: 8` counts solves and
the trigger's `Offset: 1` skips the initial guess ⇒ 7 refinements, 8 grids).

The driver runs **two lattices from one worker pool** and plots the paper's:

| | `ImposeEquatorialSymmetry` | lattice | level-0 grid | role |
|---|---|---|---|---|
| **plotted** | `false` | full-z 5×3, 30 elements | **480** points | the published one |
| secondary | `true` (yaml:36) | half 5×2, 20 elements | 320 points | the shipped file's own, registered and gated, never drawn |

### What the published figure's own vector paths say

`ssf_convergence.pdf` is a matplotlib vector figure. Its 21 top-panel and 21 bottom-panel
polylines were parsed out of the page content stream and de-projected with the calibration the
file itself carries (x gridlines at 20/40/60/80; y gridlines one decade apart, 1e-12…1e-5 on
top and the "1 s"/"10 s" ticks below). Script and full table: this node's evidence,
`published_figure_extract.py` / `published_figure_vertices.txt`.

* **Seven vertices per curve** — consistent with 8 grids and the last one used as the
  reference. An independent confirmation of the level-count mapping.
* **The gray ramp is `m/23`** (0 at m = 0 to 0.8696 at m = 20), i.e. `gray(linspace(0,1,24))`
  truncated at 21 entries. The driver uses exactly that.
* **The abscissa is the FULL-z lattice, and it now reproduces EXACTLY.** The published
  grid-point sequence is `480, 975, 1604, 2427, 3468, 4751, 6300` (m ≤ 4). Its first two
  entries are the `NumPoints` column of spectre's own
  `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml` regression at the same
  `a = 0.5 M, r₀ = 10 M` — a file whose `ImposeEquatorialSymmetry: False` gives the full-z 5×3
  block lattice, 30 elements at P = 4 = 480 grid points, with `NumContributingElements = 4`
  because the puncture sits on a four-element corner. `anc/ScalarSelfForce.yaml:36` sets it
  `True`, so the published figure was **not** run with the shipped file's own setting. With
  `crp_paper_config(; equatorial_symmetry = false)` and the all-dimension 2:1 closure, our
  plotted sequence is

      480, 975, 1604, 2427, 3468, 4751, 6300      (the reference level is 8139)

  — **the published sequence, level for level, exactly**, which GATE 5 asserts at all seven
  plotted levels. The elements go `30, 47, 59, 71, 83, 95, 107, 119`; grid 1's 47 elements are
  spectre's own `20 × (4,4) + 4 × (4,5) + 23 × (5,5) = 975`. The secondary half-domain sequence
  is `320, 632, 1014, 1506, 2120, 2868, 3762`, registered beside it.
* **The published grids still depend on m, and ours cannot.** From level 1 on the published
  `N_points` splits by mode: 975 for m ≤ 4 and 970 for m ≥ 5, spreading to 6300 (m ≤ 4) down to
  6092 (m ≥ 18) at the last plotted level — a monotone decrease with m. `RefineAtPuncture` and
  `IncreaseResolution`, the only two criteria the shipped yaml enables (yaml:130-132), are both
  m-independent, and our sweep is gated m-independent (`ndofs_identical_across_m = true`,
  GATE 3). The candidate is the yaml's own **commented-out** `TruncationError` block
  (`AbsoluteTarget: 1.e-11`, `RelativeTarget: 1.e-11`, yaml:133-136): it is solution-dependent,
  `Policies.AllowCoarsening: True` is only meaningful beside an error-based criterion, higher m
  has a smaller far-field amplitude and would coarsen more (the direction observed), the
  figure's dashed "Solver tolerance" line sits at **1e-11** (measured off the PDF) while the
  GMRES `AbsoluteResidual` is 1e-12, and main.tex:721 says the error falls "until it reaches the
  tolerances set for the GMRES algorithm **and for AMR**" — a second tolerance the shipped
  configuration does not have. **[HOLE: consistent with every measurement here, not
  demonstrated; no run of ours enables that criterion.]** It is now the ONLY part of the
  abscissa that does not reproduce, and it is a 0.5 % effect (970 against 975) rather than the
  50–78 % lattice mismatch it was hiding behind.

### Which lattice the gates apply to

The published-envelope R² floor (GATE 1a) is **scoped to the plotted lattice**, and the reason
is physical, not cosmetic: `R2_FLOOR` is calibrated on the published curves, which descend ~7
decades on grids reaching 6300 points. At the same 8 AMR steps the half domain only reaches
3762 points and its curves floor at ~2e-11 against the plotted lattice's ~2e-13 — a floor
**49…839× higher**, registered per mode as
`floor_ratio_secondary_over_plotted_per_mode`. On a curve whose tail is already floor, R²
measures what fraction of the fit range is floor rather than how exponentially the method
converges. The lattice-independent half of the claim — every mode's error actually decays
(GATE 1b) — is applied to both, as are the decade floor, the mesh m-independence, and the
per-level GMRES convergence. The secondary series' fits are registered
(`secondary_fit_r2_per_mode`, `secondary_fit_r2_min`, `n_modes_secondary_r2_above_098`) so the
contrast between the two lattices is in the golden rather than in prose.

### The planned R² gate is falsified by the published figure itself

GLOBAL_DAG row 50 planned "21 exponential fits with R² > 0.98 down to a floor ≤ 1e-11". Fitting
`log10(err)` against `sqrt(N_points)` on the **published** curves gives R² ∈ [0.7996, 0.9929],
median 0.9833, with **9 of 21** modes at or above 0.98 (worst: m = 15 at 0.7996, whose level-2
point dips three decades below its neighbours and returns). The landed gates are therefore the
published envelope and a floor on the total drop:

| gate | instrument | landed | measured (21 modes) |
|---|---|---|---|
| 1a (plotted) | least-squares fit of `log10(err)` on `sqrt(N_points)`, 7 points | rate < 0 **and** R² ≥ 0.80 (= the published minimum, rounded down) | rate ∈ [−0.1395, −0.1020] (published [−0.1421, −0.1043]); R² ∈ [0.9143, 0.9932], median 0.9781, 10/21 ≥ 0.98 |
| 1b (both) | the same fit's rate | rate < 0 on every mode of BOTH lattices | true (secondary rates ∈ [−0.1165, −0.0915]) |
| 2 (both) | `log10(err_first/err_last)` | ≥ 2.0 decades (calibration prototype's worst was 2.65) | plotted [6.02, 8.38] (published [6.06, 8.12]); secondary [4.04, 5.60] |
| 3 (both) | `ndofs` per level, all 21 modes | identical | true |
| 4 (both) | `SolveRecord.converged` on every level of every mode | all true | true |
| 5 (plotted) | `ndofs` per level against the published vertices | equal at all **seven** plotted levels, not just the first | true — 480, 975, 1604, 2427, 3468, 4751, 6300 |
| 6 | tracked PNG | ≤ 200 KB | 155 219 B |

Of the two planned numbers, one is now MET and the other is falsified by the paper's own
figure. **The floor ≤ 1e-11** is met on the published lattice: our last plotted error is
**5.6e-14 … 7.6e-13** across the 21 modes, against the published 8.9e-14 … 7.9e-13 — the same
order, on the same 6300-point grids. (On the half-domain lattice, which stops at 3762 points,
it is 1.8e-11 … 5.6e-11, and that gap is the floor ratio GATE 1a is scoped around.)
**R² > 0.98 for all 21** is not met and cannot be: the paper's own figure reaches it for 9 of
21 and ours for 10 of 21, and both fail it for the same structural reason — when a coarse
level's value happens to land near the reference, its "error" is anomalously small and breaks
the straight line. It is landed as the published envelope instead.

### The two-step contraction envelope is reported, not gated

The row-44 contract states the convergence of `F_r^m` as a two-level envelope — each difference
at most half the one two levels back. That is a statement about the **pre-floor** regime, and
at 8 AMR levels the high-m sequences are already on their resolution floor
(`|F_r^20| = 2.1e-9` against a floor of 1.7e-13), where the differences are floor noise.
Measured on the plotted lattice, **1 of the 105** two-step ratios exceeds one half — m = 10, at
**0.7029** — and **none exceeds one**, so the row-44 half-envelope holds for every mode except
m = 10 and the weaker unit contraction holds for all 21. The coarser half domain, whose floor is
49…839× higher, breaks the half-envelope on 7 of 21 modes and reaches 3.09: the same envelope
measured on a sequence that is mostly floor, which is exactly why it is measured and not gated.
Gating either would be gating noise, so
`two_step_ratio_max_per_mode`, `n_two_step_ratios_above_half` and
`modes_with_two_step_above_half` are registered in the golden and neither is asserted.

### Walltime is reported, never gated

The bottom panel is the paper's own cumulative walltime per mode. Published (Apple M2 Pro,
read off the PDF): 0.154 … 0.476 s at the first grid, 4.19 … 5.93 s at the last — "a few
seconds per m-mode" (main.tex:688). Ours, on a shared 256-core Linux box with one BLAS thread
per process and five worker processes competing for it, is 0.50 … 1.39 s at the first grid and
23.4 … 28.7 s at the last — roughly 5× the published figure's, now on grids of exactly the same
size rather than smaller ones, and it swings by 30 % with the machine's load, which is why every
`reported_*` timing row carries a deliberately loose tolerance and none of them is gated.
The driver holds ONE `mode_pool` for the whole figure and runs a throwaway warm-up sweep first,
because a fresh worker charges its entire compilation to its first mode's level 0 — without it
the bottom panel measures the compiler. The warm-up is itself the JIT payment and dominates the
run.


## Fig. 5 — what is plotted, which radii, and the rule that picks the horizon hold

### The figure plots the FIELD, not a mode-sum curve

`fig:modes` is not `|F_r^m|` against `m`. Its caption (main.tex:691-696) is explicit — "Each
plane shows an m-mode of the scalar field, `Ψ_m`, up to m = 20 … Plotted is their complex
magnitude. Within the regularized region, the regular field `Ψ^ℛ_m = Ψ_m − Ψ^𝒫_m` is plotted.
The computational domain extends over `z² ∈ [0, 1]` (equatorial symmetry) and to `R = 10⁴ M`
towards the right, but the plot is clipped at `r = 2 r₀`" — and main.tex:718 adds what the
figure demonstrates: "the decay in power over mode number" and that "regular zero crossings in
the solution of the regularized field within the worldtube display as sharp drops on the chosen
log scale, particularly for the m = 0 mode". So each panel is 21 surfaces over the (r, z²)
domain, which is exactly `sol.u` under the row-42 contract (`Ψ^ℛ` inside the worldtube, `Ψ`
outside). The driver draws the element nodal values themselves — no interpolation — which is
also what the published ParaView render shows, patch by DG patch.

### The panels: two ISCOs and the published prograde radius

| panel | published title | `kerr_isco(a)` | what we plot |
|---|---|---|---|
| a = −0.998 | `r₀ = 9` | 8.9943744548 | **8.9943744548** — the ISCO (0.06 % from the title) |
| a = 0 | `r₀ = 6` | 6 | **6** — the ISCO, and the title |
| a = +0.998 | `r₀ = 2` | 1.2369706552 | **2** — the PUBLISHED title (1.6× the ISCO) |

Each published abscissa ends exactly at twice its title's radius (18, 12, 4 — the caption's
"clipped at `r = 2 r₀`"), so the published panels were run at 9, 6 and 2, and only the middle
one is its spin's ISCO. Panels 1-2 follow the caption, panel 3 follows the published title —
which is what the paper actually drew — and every panel title carries the radius solved.

`r₀ = 2 M` was **not constructible** until `47f0ec0`: `_angular_momentum_circular`
(`constants.jl:13`) divides by `(p − 2) p − a²(x² − 1)`, which on the equator is `p (p − 2)` and
vanishes at `p = 2` with its numerator, and `_carter_circular` (`:17`) inherits it as `0·NaN`,
so **`Lz` and `Q`** were 0/0 while **`E`** was already right; `zc2` (main.tex:453) took the NaN
and the puncture died in `legendre_P3`. It only bit where `r = 2 M` carries a circular geodesic
at all, `a > 1/√2 = 0.70710678` (below that the negative radicand makes the `DomainError`
correct; at `a = 0`, `p = 2` is the horizon, where the pre-fix code returned a silent `E = 0.0`).
Fixed bitwise-neutrally away from the endpoint, with `7082f9a` pinning the orbit and the
puncture at `r₀ = 2 M` and `12b43bd` excluding `p = 0`.

### The worldtube radius is a documented choice

The yaml is templated in `worldtube_radius` and the paper prints one only for its header
example (`r₀ = 10 M`, `r_in = 2 M` ⇒ 4.8). That value is illegal here — the lattice needs
`r₀ − 1.5 w_r > r₊` for a non-empty v block, and 4.8 puts the v→t transition at −1.2 M for
`r₀ = 6`. The driver scales the worldtube with the radial gap the orbit actually has,

    w_r = 0.6 (r₀ − r₊),

which reproduces the paper's own published value EXACTLY at its own point (0.6 (10 − 2) = 4.8),
keeps the v block at a fixed 10 % of `[r₊, r₀]` at every spin, and is the only one-parameter
family in reach that stays legal at all three ISCOs. It is an ASSUMPTION where the paper prints
nothing — and the published figure corroborates it: the grey curtain's cliff, which is the
discontinuity between `Ψ^ℛ` inside the tube and `Ψ` outside, sits at `r₀ + w_r` = 13.75 and
8.4 M in panels 1 and 2 against ≈ 13.5 and ≈ 8.5 M read off the published panels.

### The horizon–axis layer, and the measured rule that decides the hold

**What regularity means here, stated correctly.** `eq:mode_sum` factors `sin^{|m|}θ` **out** of
the field, so `Ψ_m` is the variable that is analytic at the axis and approaches a finite,
generically non-zero limit there — main.tex:227: each mode "approaches a non-zero constant
amplitude … towards the poles". A healthy plane therefore has an **axis:equator ratio of order
one** and a fitted `log|Ψ_m|` vs `log(1 − z²)` slope near **zero**. It must *not* fall like
`(1 − z²)^{m/2}`; that is `Φ`'s m-mode, never `Ψ_m`. (Earlier revisions of this section had it
the wrong way round; the measurement was always a ratio, so nothing measured changed.)

**What goes wrong.** At the horizon the operator degenerates and its second radial branch is
`Ψ_m ∼ (r − r₊)^{iν}` with `ν = 2m|Ω − Ω_H|(r₊² + a²)/(r₊ − r₋)` — unit modulus, infinitely
oscillatory, excluded **only** by the polynomial trial space. Refine the horizon element away
and it appears as a layer on the innermost radial node line. `ν` is 341.6 at m = 20 on panel 1,
5.4 on panel 2 and 139.9 on panel 3, and the layer is present at every spin: it is a question of
**scale**, so the driver never keys on the spin.

**The rule.** Sweep the panel with the package default (`horizon_hold = false`, reference
parity), read `horizon_axis_profile(dom, orbit, mesh, u).layer_ratio` — the near-horizon strip's
maximum over the rest of the plane's — for **every** mode, and re-sweep with
`solve_kwargs = (; horizon_hold = true)` only if the **maximum over m** exceeds `1.05`. The
plotted run's own `layer_ratio` is then gated (**GATE 7**), and both numbers and the setting are
stated on the panel and registered.

**The maximum over m, not the top mode** — measured, and it cost a run to learn: at
`a = +0.998, r₀ = 2` the `layer_ratio` is 0.58 at m = 0, crosses one at m = 6, **peaks at 3.889
at m = 12** and falls back to 3.6e-3 by m = 20, so a top-mode probe reads the panel as healthy
(0.00363) and leaves the hold off.

**What the rule decided here**: panels 1 and 2 pass on the default (`layer_ratio` max 0.1105 and
0.3108, both at m = 0, i.e. no layer at all), panel 3 fails it (3.889) and is re-swept held, which
brings it to **0.5776** — GATE 7 passes on all three. The hold changes the plane by up to 1.4
decades and leaves the particle alone: at this configuration `F_r^m` and `|Ψ^ℛ(x_p)|` agreed
between the two settings to 9.0e-4 and 2.6e-4 relative.

### [HOLE] GMRES stagnates at m = 19 and 20 for `a = +0.998`, `r₀ = 2 M`

`converged = false` on those two modes with a worst final residual of **9.3e-6** against the
yaml's 1e-12 tolerance, stopping far short of the 500-iteration cap; every other mode of every
panel converged. Pre-adoption measurements of the same configuration gave 2.3e-5 (default) and
7.6e-5 (held) at m = 20, with the per-grid staircase **growing** with depth
(`1, 7, 7, 7, 11, 21, 61, 235` at m = 19). The two planes stay in the figure — drawn **dashed and
fainter**, named in the panel note with the plotted run's own residual — and `pK_converged`
carries the verdict per mode; the only break in panel 3's decay claim is the 19 → 20 step.
Handed to `ESSF::fix-horizon-axis-layer`.

### Gates: one hard gate on the rule, the paper's claims registered

**GATE 7 is a real gate**: the plotted run's `layer_ratio` must not exceed 1.05 on any mode. The
driver *chooses* `horizon_hold` to satisfy it, so a violation means the rule failed and the panel
is not drawable as physics. Alongside it, `gatefails` holds only structural faults: a lost mode,
a violated precondition of the extraction (`dom.chart === :z2` and `x_p[2] == 0`, which is what
makes the sealed raw-ingredient `mode_self_force(orbit, m, psi, dpsi_r)` bitwise equal to the
`AmrHistory` method's `:trace` reduction), and an oversized tracked PNG.

The paper's two claims are **measured per panel and registered**, not asserted — the fig_linsolv
precedent — so `--check` keeps comparing numbers and a later change that breaks a claim which
currently holds fails the comparison:

* `pK_decay_monotone_from_m2` — `max|Ψ_m|` strictly decreasing for m ≥ 2 (true, true, false),
  with `pK_decay_first_break_m` (−1, −1, 19) and the particle-local
  `pK_psi_xp_monotone_from_m2` (**true on all three** — with the hold, panel 3's particle-side
  decay is monotone and only its plane's 19 → 20 step breaks, on the two stagnated modes);
* `pK_worldtube_decades` / `pK_outside_decades` with `pK_m0_drops_two_decades` (2.43, 3.22, 3.04
  decades) and `pK_m0_deeper_than_m1_m2` (**true on all three**: m = 0 drops 2.43 against
  0.95/1.23, 3.22 against 0.87/1.21, 3.04 against 1.16/1.85). The honest part stays registered:
  m = 0 is not the deepest plane of all (`pK_deepest_worldtube_mode` = 18, 12, 20);
* the row-58 instrument and the solver's verdict, per mode: `pK_layer_ratio`,
  `pK_axis_equator_ratio`, `pK_axis_exponent`, `pK_branch_index`, `pK_converged`,
  `pK_final_residual`, plus the rule's inputs `pK_horizon_hold` and
  `pK_layer_ratio_default_max`.

### The lattice these numbers belong to

Everything here is the production `ImposeEquatorialSymmetry: True` (yaml:36) `z²` half-domain —
320 dofs at level 0, 4535 at the 8th grid, 63 elements — registered as
`impose_equatorial_symmetry`, `pK_level0_ndofs` and `pK_final_ndofs`. ESSF::fig-convergence
found that the published fig:convergence abscissa implies the paper's own runs used the FULL-z
lattice (`ImposeEquatorialSymmetry: False`, 480 points at level 0); the physics is symmetric, so
the per-m planes are unaffected beyond resolution, but every dof/element count above is
lattice-dependent and should be read as this lattice's.

## Fig. 7 — what the published figure actually samples, what it plots, and the lattice that sets the ordering

### The published figure's own vector paths

`ssf_error_budget.pdf` is a matplotlib PDF, so its eleven data polylines are readable
exactly. Calibrating on the axis tick and gridline STROKES (x: `a = 0` at 144.05 / 370.01 pt
with 97.03 pt per unit spin; y: 29.333 pt per decade in the budget panels and 8.075 pt per
decade in the strips) gives the published curves to about 0.1 %, which is checked
independently: the published `|F_r|` strip at `r₀ = 10 M` reproduces the printed
Warburton–Barack table to **0.10 %** at all seven of its points. The extraction script and
the calibration probe are in this node's evidence
(`published_curve_extraction.py`, `published_axis_calibration.py`, `published_curves.json`).

Three structural facts follow, none of them guessable from the rendered image:

1. **Nine spins are sampled, not five.** The five tick labels (−0.998, −0.5, 0, 0.5, 0.998)
   are not the sample set: markers sit at `a = −0.998, −0.9, −0.7, −0.5, 0, +0.5, +0.7, +0.9,
   +0.998`. The four extra spins are exactly the ones that carry a printed WB cell. This
   driver runs the published nine.
2. **The solid curve stops at |a| = 0.9.** The dashed (truncation) and dotted (resolution)
   curves carry all nine points; the solid "Error to reference" curve and the `|F_r|` strip
   carry only seven, `|a| ≤ 0.9`. The marker count settles it — 2×7 + 2×7 + 2×9 + 2×9 + 3
   legend = **67** filled circles, and the PDF contains exactly 67. So the published figure
   has **no reference point at |a| = 0.998 either**: its reference stops where
   Warburton–Barack's printed range stops. Our blank column there is the paper's own blank,
   not a gap of this reproduction.
3. **The published `F_r` is the BARE partial sum — no fitted tail.** The published
   reference/truncation ratio is 2.81–4.63 at `r₀ = 10` and 3.09–3.61 at the ISCO, i.e. the
   "constant offset" of main.tex:726. That constant is not free: for a `C m^{-p}` tail the
   omit-the-largest-`m` estimator is smaller than the true remainder `Σ_{m>m_max}` by
   `≈ m_max/(p−1)`, which `ESSF::mode-sum-total-force` measured as 3.57–3.65 on this very
   data. If a tail had been added the remainder would have been removed and the offset could
   not be 3.6. Two direct confirmations from that node's sealed table: our **no-tail**
   `|F_r − F_r^WB|` is 3.6e−6 at `a = −0.9` against the published 3.53e−6, and 3.4e−5 at
   `a = +0.5` against the published 3.42e−5 — 2 % and 0.6 %.

The driver therefore plots the paper's own quantity, `total_self_force(set; tail = :none)`,
and registers the Heffernan-tail number beside it (`Fr_heffernan_tail_flat`,
`rel_ref_error_tail_flat`) — the same solved modes, one extra call to the sealed mode sum,
no extra solve. The tail-corrected `F_r` is this reproduction's more accurate number and its
reference error is roughly 5–10× smaller than the plotted one.

### The lattice is what sets the ordering — the row-46 finding, resolved

`ESSF::mode-sum-total-force` measured, on the shipped `ImposeEquatorialSymmetry: True` half
domain, a coherent resolution error 3–6× **above** the truncation estimator — the opposite of
main.tex:723 — and left the cause open. `fig_convergence.jl` had already gated that the
published convergence figure's abscissa is the **full-z 5×3 lattice** (480, 975, 1604, 2427,
3468, 4751, 6300 grid points, level for level), so the published figures were not run with the
shipped file's own setting of `yaml:36`. This node ran the same configurations on both
lattices from one pool and measured the consequence on the TOTAL.

### The worldtube radius

`anc/ScalarSelfForce.yaml:12` templates it and the paper publishes a value only for its header
example (`r₀ = 10 M`, `r_in = 2 M`: `worldtube_radius = 4.8`). At `r₀ = 10` that published
value is legal at every spin (`v→t = 10 − 1.5·4.8 = 2.8 > r₊ ≤ 2`) and is used verbatim, which
also keeps this column directly comparable with row 46's sealed table. At the ISCO it is not:
4.8 puts the `v→t` transition at −1.2 M for `r₀ = 6`. The ISCO column therefore uses
`fig_modes.jl`'s landed one-parameter family, `w_r = 0.6 (r₀ − r₊)`, which reproduces the
paper's own 4.8 exactly at its own point and keeps the v block at a fixed 10 % of `[r₊, r₀]`
at every spin. It is an assumption where the paper prints nothing, it applies to the ISCO
column only, and every ISCO number registered here depends on it.

### The horizon order hold, measured — with one strengthening

The rule is `fig_modes.jl`'s: sweep with the package default (`horizon_hold = false`,
reference parity), read `horizon_axis_profile(dom, orbit, mesh, u).layer_ratio` on **every**
mode, and re-sweep with `solve_kwargs = (; horizon_hold = true)` only when the **maximum over
m** exceeds 1.05. This figure adds one thing that node could assume and this one cannot: when
the default probe crosses the gate the driver re-sweeps held **and keeps whichever setting
measures smaller**. Both maxima are registered (`layer_ratio_default_max_flat`,
`layer_ratio_held_max_flat`, `horizon_hold_flat`).

The hold is not free everywhere: measured in this node's calibration prototype, at
`a = −0.998` at its ISCO the held sweep changes `F_r^{m=0}` by 0.4 % relative
(4.0308e−4 → 4.0144e−4) while its default `layer_ratio` is only 0.110 — which is exactly why
the rule is a measured switch and not a default.

### What reproduces, measured

| claim | ours | published | verdict |
|---|---|---|---|
| relative truncation error, all 18 configurations | see the golden | read off the PDF | **0.3–3.6 %** at 17 of 18 (the exception is the compressed near-extremal ISCO below) |
| relative reference error, all 8 printed WB points | 3.06e−6 … 3.50e−5 | 3.53e−6 … 3.54e−5 | **ratio 0.775–0.988** — slightly more accurate than the paper at every point |
| reference/truncation offset (main.tex:726) | 2.28–3.75 | 2.81–4.63 | analytic `m_max/(p−1) ≈ 3.6`; row 46 measured 3.57–3.65 |
| `\|F_r\|` against the printed WB table, r₀ = 10 | 7 points | 7 points | all 7 to **≤ 1.2e−5 relative** |
| unregularized `F_t` against Table II (9 printed digits) | 8 points | — | all inside 1e−6 relative, the second external gate |
| main.tex:723 ordering, resolution below truncation | **15 of 18** on the paper's lattice | 18 of 18 | **0 of 9** on the shipped half domain at the same 8 grids |

**The row-46 finding is closed as a lattice effect.** `ESSF::mode-sum-total-force` measured, on
the half domain, a resolution error 3–6× ABOVE the truncation estimator and left the cause open.
Running the same nine `r₀ = 10` configurations on both lattices from one pool, at the same eight
AMR grids, gives the controlled variation: the two lattices agree on `F_r` to 1.18–1.26e−10
absolute (the equatorial-symmetry reduction is exact, so they must) while their resolution
errors differ by 88–2370× at `|a| ≤ 0.9`. On the paper's own lattice the ordering of
main.tex:723 holds; on the shipped one it does not, at every single configuration.

### The gates, and the one that was replaced after measurement

GLOBAL_DAG row 53 planned "the reference error within 3× the coherent resolution error". On the
paper's lattice that is unreachable BY CONSTRUCTION, and its failing is the paper's headline
claim rather than a defect: with the resolution floor two to three decades below the mode-sum
remainder, `|F_r − F_r^WB|` IS the remainder. **Measured: 0 of the 8 WB-referenced
configurations satisfy it**, registered as `n_wb_configurations_within_dag_factor`, reported
and never gated. What replaced it is the published figure's own envelope — the same discipline
`fig_convergence.jl` landed for its R² floor:

* **GATE 3** relative truncation within a factor 2 of the published dashed curve (measured
  worst 1.036);
* **GATE 4** relative reference error within a factor 3 of the published solid curve at every
  printed WB point (measured worst 1/0.775 = 1.29);
* **GATE 5** the reference/truncation offset inside [2, 6] (measured 2.28–3.75);
* **GATE 6** the unregularized `F_t` within 1e−6 relative of the printed Table II — no
  regularization and no tail, so it separates a solve/extraction error from a
  regularization/tail one (row 46 §3b);
* **GATE 2** the plotted `layer_ratio` at or below 1.05;
* **GATE 7** the tracked PNG inside the shared 200 KB budget.

GATES 2–6 are scoped to the configurations whose mode sum has actually CONVERGED at `m_max`,
by two measured criteria: the sealed tail estimator refusing to fit, and the last retained
regularized term reaching 1 % of the total. **The second criterion was added after the first
production run** and is recorded as such here and in the driver: the first criterion alone does
not fire at the one configuration that needs it, and its threshold is calibrated on the
PUBLISHED dashed curve (whose largest value anywhere is 2.93e−4, 34× below the bound), never on
ours. Exactly one configuration is excluded, and it is ring-marked on the figure.

### [HOLE] the compressed near-extremal prograde ISCO, and the `|a| = 0.998` resolution error

At `a = +0.998`, `r₀ = r_ISCO = 1.2370 M` the mode sum has not converged at `m = 20`: the last
retained term is 2.4 % of the total and the resolution error 9.5 % of it, against a published
truncation of 1.99e−4. Every mode still reaches the 1e−12 GMRES tolerance and the horizon order
hold fires (default `layer_ratio` 23.94 at `m = 17`; the rule keeps whichever setting measures
smaller). The point stays in the figure, ring-marked.

Separately, at `|a| = 0.998`, `r₀ = 10 M` the relative RESOLUTION error is 1.48e−5 and 3.92e−5
against a published 1.99e−7 and 3.92e−7 — 74× and 100× too large — and those are the only two
other configurations where our ordering inverts. Every other spin at `r₀ = 10` sits between
1.4e−8 and 1.0e−6. The `layer_ratio` instrument reads only 0.098 at both, so whatever costs the
resolution there is NOT the near-horizon layer that instrument measures. Cause not established;
handed to `ESSF::fix-horizon-axis-layer`.

## Not yet reproduced

Figure 3 (`precond`, the preconditioner schematic) is a later node of the
elliptic-ssf DAG. Nothing about it is claimed here.
