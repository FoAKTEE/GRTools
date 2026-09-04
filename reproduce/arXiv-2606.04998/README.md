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
| 4 | `fig:convergence` — the DG convergence of every `F_r^m` up to m = 20 over the hp-AMR sequence, + the cumulative walltime | `fig_convergence.jl` | `expected/fig_convergence.json` | heavy | **preliminary** | **partial** | Two panels for all 21 modes at a = 0.5 M, r₀ = 10 M: `\|F_r^m − F^m_{r,ref}\|` against `sqrt(N_points)` with each mode's own highest-resolution level as its reference, and the cumulative walltime below. Presentation reproduces item by item (2 panels 3:1 sharing x, 21 grayscale lines at the published ramp `m/23` with no markers, log y over ~8 decades, the framed upper-right `m=0`/`m=20` legend, the dashed "Solver tolerance" line with its inside label, "1 s"/"10 s" ticks, light-gray grid, CM serif) and every mode converges exponentially. Two STRUCTURAL deviations, both escalated: the published abscissa is the **full-z 5×3 lattice** (480 → 6300 grid points; its first two values are spectre's own `ScalarCircularOrbit.yaml` regression rows, which set `ImposeEquatorialSymmetry: False`) while ours is the production equatorial-symmetry lattice (320 → 3538), and the published grids **depend on m** while ours cannot, which the shipped criteria (`RefineAtPuncture` ∘ `IncreaseResolution`) cannot produce. See the Fig. 4 section for the evidence and the `TruncationError` hypothesis. |
| 6 | `fig:linsolv` — GMRES residual staircase over AMR levels + walltime per iteration | `fig_linsolv.jl` | `expected/fig_linsolv.json` | heavy | **preliminary** | **partial** | Two panels for the m = 2 mode at a = 0.5 M, r₀ = 10 M: the residual staircase of 8 consecutive AMR levels over a cumulative-iteration axis, and the walltime per iteration. Presentation and structure reproduce item by item (see the checklist below); the per-level ITERATION COUNTS do not — published `1, 7, 6, 5, 4, 3, 3, 2` (31 cumulative) against our `1, 8, 9, 9, 9, 9, 9, 9` (63). All five node gates pass (max 9 < 10 iterations per level; level 0 in 1 iteration; the count bounded across levels; the first iteration slowest on every level; every level below the 1e-12 tolerance), so the paper's stated claims hold — but the shape of the published staircases does not, which is a solver-quality finding escalated to rows 21 and 43, and why this row is `preliminary`. 252.7 s at -t4 in `run_all` with a warm depot (~400 s standalone with precompilation), of which the AMR sweep itself is 12.3 s. Generated on top of `6df9f51` (the `hp_refine` 2:1-closure fix). |
| 5 | `fig:modes` — 21 stacked `log₁₀\|Ψ_m\|` planes over the (r, z²) domain, at the ISCO for three spins | `fig_modes.jl` | `expected/fig_modes.json` | heavy | **preliminary** | **partial** (2 of 3 panels match) | Three 3D panels, `a ∈ {−0.998, 0, +0.998}` at `r₀ = kerr_isco(a) = 8.9944, 6, 1.2370 M`, each carrying 21 semi-transparent wireframe planes of `log₁₀\|Ψ_m\|` over (r, z²) — `Ψ^ℛ_m` inside the worldtube, `Ψ_m` outside, clipped at `r = 2 r₀` — with m = 0…4 in the paper's categorical colours and m ≥ 5 in the published `m/23` grey ramp. NOT a mode-sum curve: the caption plots the FIELD. Panels 1 and 2 **match** item by item, including the grey curtain's cliff at `r₀ + w_r` (13.75 against ≈ 13.5 read off the published panel; 8.4 against ≈ 8.5), the two m = 0 dives per panel, and the two-sided 0…−5 ticks. Panel 3 is **preliminary** and says so on the figure: at the prograde ISCO a **spurious boundary layer on the innermost radial node line `r = r₊`** grows toward the axis and takes over the plane from m ≈ 14 — at m = 17 along `r = r₊` the field RISES from `\|u\| = 1.6e-5` at `z² = 0.003` to `1.05` at `z² = 0.9918`, where axis regularity (`sin^\|m\|θ` is factored out, so the plane must fall like `(1 − z²)^{m/2}`) asks for 0.98 down to 1.8e-18 — a 23-decade violation; one node in (`r > r₊ + 0.02`) the same modes decay physically at −0.087 dec/mode with GMRES converging in 10–11 iterations. It is present at every spin and is a matter of scale (inner-boundary axis:equator ratio at m = 17: 1.3 at a = 0, 4.5e3 at a = −0.998 but at an absolute 3.7e-8, four decades below that panel's maximum, and 6.6e4 at O(1) here), and it is owned by `ESSF::fix-horizon-axis-layer`. Separately, the published prograde panel is at `r₀ = 2 M` — **not** its spin's ISCO — which was unconstructible when this figure landed (`Lz` and `Q` were 0/0 at `p = 2` on the equator) and is constructible since `47f0ec0`; regenerating that panel waits on the layer fix. Both claims are registered per panel, so `--check` compares them instead of aborting. 63 modes (3 × 21) at the paper's 8 AMR grids through ONE `mode_pool` of 5 single-threaded workers: 561 s all told, of which the first panel's 378 s is the pool's JIT (the two warm panels take 86 s each); whole-system load max 0.060 against the 0.90 guard — **heavy**. |

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

## Fig. 6 — how many AMR levels, and what does not reproduce

### The level count is a convention mapping, derived from the source

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

`ESSF::amr-driver` ran `iterations = 8`, i.e. **nine** levels; the ninth `9` of its published
sequence `1, 8, 9, 9, 9, 9, 9, 9, 9` is that extra level. The first eight entries are this
figure's.

### The finding: the staircases do not shorten

| level | iters (paper) | iters (ours) | initial residual paper / ours | reduction per iteration paper / ours |
|---|---|---|---|---|
| 0 | 1 | 1 | 1.375 / 1.744 | — |
| 1 | 7 | 8 | 2.59e-3 / 2.53e-3 | 26.0 / 16.4 |
| 2 | 6 | 9 | 7.06e-4 / 3.26e-3 | 42.7 / 12.2 |
| 3 | 5 | 9 | 8.75e-5 / 3.64e-4 | 48.0 / 9.0 |
| 4 | 4 | 9 | 2.37e-5 / 2.18e-4 | 82.3 / 9.3 |
| 5 | 3 | 9 | 4.01e-6 / 5.12e-5 | 244.6 / 7.6 |
| 6 | 3 | 9 | 9.64e-7 / 1.03e-5 | 191.9 / 6.4 |
| 7 | 2 | 9 | 2.60e-7 / 5.30e-6 | 567.0 / 6.0 |
| **total** | **31** | **63** | | |

The paper column is extracted from `linsolv.pdf`'s own **vector paths** (PyMuPDF; axes calibrated
on the tick-label boxes and the bottom panel's minor-tick strokes, self-consistent to <1 %), not
from pixels. "Reduction per iteration" is `(r₀ / r_final)^(1/iters)`.

Both runs satisfy every claim the caption and main.tex:720 make — `< 10` iterations per level, the
first level converging immediately, and no growth of the count as levels are added. What differs
is that the paper's preconditioner becomes steadily **more** effective with depth while ours
**degrades** (16.4× → 6.0×), and our per-level initial residual is larger at **every** level: at
level 0 it is the source norm `‖b‖` itself, 1.744 against the published 1.375 (a factor 1.27 —
an open [HOLE], see below), and from level 2 on the interpolated warm start is 1.4–13× larger. Owner: `ESSF::solver-stack-core` (the GMRES ∘ multigrid ∘ Schwarz composition)
and `ESSF::amr-driver` (the level policy). Row 43 already measured a nearby A/B: with
`coarse_orders = nothing` (the min rule for the coarse multigrid orders, instead of spectre-parity
pinning at `P`) the same sweep gives `1, 8, 8, 8, 7, 7, 6, 7, 6` — decreasing, though still not
the published shape. That is the first lead to follow.

### The golden is post-`hp_refine`-fix

It is generated on top of `6df9f51` (*fix(ellipticcore)!: hp_refine — hold the order in a
balance-split dimension*), where an element that the 2:1 balance splits no longer also carries a
p increment in the split dimension, matching spectre
(`ParallelAlgorithms/Amr/UpdateAmrDecision.cpp:82-87`, `Amr/Projectors/Mesh.cpp:30-34`).

Measured effect on this figure: the dof sequence moves from
`320, 603, 962, 1425, 2004, 2711, 3558, 4557` to `320, 593, 950, 1411, 1988, 2693, 3538, 4535`,
and the per-level final residuals move in their last digits. The GMRES **iteration counts are
unchanged** at `1, 8, 9, 9, 9, 9, 9, 9`, so the reproduction gap documented above survives the
fix rather than having been caused by it.

### [HOLE] the level-0 source norm is 1.27× the published one

Level 0's initial residual is `‖b‖` itself (the initial guess is zero), and it is **1.744** here
against the published **1.375** — a factor 1.27 that then propagates as an offset into every
later level's starting residual.

**Mass weighting is excluded as the cause.** The candidate explanation was that spectre's Krylov
norm might be a mass-weighted inner product while ours is Euclidean. It is not: spectre's Krylov
inner product is the plain Euclidean 2-norm of the *massive* residual
(`NumericalAlgorithms/LinearSolver/InnerProduct.hpp:23-27,90-96`, with the mass matrix already
applied by the DG operator under `Massive: True`,
`Elliptic/DiscontinuousGalerkin/DgOperator.hpp:621-680`), and ours plots the same quantity —
the Euclidean norm of the massive residual (`Solve/driver.jl:229-232`,
`EllipticCore/src/gmres.jl:144`). Both sides therefore measure the same norm of the same kind of
vector, so the ratio is **not** a norm-convention artefact.

What remains open is the **normalisation of the source itself** — the fixed-source lift built
from the effective source and the worldtube jump data, and any symmetry/mode factor (the
equatorial-symmetry ½, the `e^{imφ}` convention, the ±m folding) that would scale `b` by an
O(1) constant without touching the *shape* of any staircase. A constant rescaling of `b` cannot
change GMRES iteration counts against a *relative* tolerance, but this solve uses an **absolute**
one (`AbsoluteResidual: 1e-12`), so a 1.27× larger `b` costs a fraction of one iteration — far
too little to explain the 63-vs-31 gap, which is why the two findings are separate.

Open, unresolved, and handed to node 56 `linsolv-depth-gap` together with the depth gap.

### The bottom panel: eager build, attributed by addition

The paper's first-iteration bump is the **lazy** construction of the Schwarz subdomain matrices
inside that iteration (main.tex:702). This package builds them **eagerly** in `rebuild!`, so the
same work is measured in `SolveRecord.build_time` and *not* inside `iteration_time[1]`. The bottom
panel therefore plots `iteration_time[1] + build_time` at each level's first point: the same
physical work charged to the same iteration, by addition rather than by measurement inside the
timer — a bookkeeping equivalence, not a lazy build. Both columns stay separately registered in
the golden, and the figure itself carries the note *"first point of each level includes the eager
operator build"* under the bottom panel, so the panel cannot be read as a lazy-build measurement. Measured consequences (envelope over the `-t1` and `-t4` runs): with the attribution the
first/rest ratio is **37…144** (paper: 1.8…2.9); without it, **0.35…3.1**, i.e. scattered about 1
with no trend in level — there is no first-iteration feature here at all without the attribution. Our per-iteration plateau is *faster* than the paper's
(0.034 s vs 0.274 s on the deepest level); only the eager build makes the sweep total larger
(12.3 s against the published stamp's 0:00:06).

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
| 1 | least-squares fit of `log10(err)` on `sqrt(N_points)`, 7 points | rate < 0 **and** R² ≥ 0.80 (= the published minimum, rounded down) | rate ∈ [−0.1372, −0.0807] (published [−0.1421, −0.1043]); R² ∈ [0.8212, 0.9954], median 0.9716, 5/21 ≥ 0.98 |
| 2 | `log10(err_first/err_last)` | ≥ 2.0 decades (calibration prototype's worst was 2.65) | [2.65, 5.82] (published [6.06, 8.12]) |
| 3 | `ndofs` per level, all 21 modes | identical | true |
| 4 | `SolveRecord.converged` on every level of every mode | all true | true |
| 5 | tracked PNG | ≤ 200 KB | see the golden |

Two planned numbers are NOT met and are landed as measured envelopes rather than dropped:
**R² > 0.98 for all 21** (the paper's own figure reaches it for 9 of 21; ours for 5 of 21 —
both fail it, and for the same structural reason: when a coarse level's value happens to land
near the reference, its "error" is anomalously small and breaks the straight line), and the
**floor ≤ 1e-11** (our last plotted error is 4.3e-12 … 3.0e-10 across the 21 modes; the
published one is 8.9e-14 … 7.9e-13, reached on grids of 6300 points against our 3538).

### The two-step contraction envelope is reported, not gated

The row-44 contract states the convergence of `F_r^m` as a two-level envelope — each difference
at most half the one two levels back. That is a statement about the **pre-floor** regime. At 8
AMR levels the high-m sequences are already on their resolution floor from level ~5
(`|F_r^20| = 2.1e-9` against a floor of 1.7e-10), where the differences are floor noise:
measured, 3 of the 105 two-step ratios exceed one half — at m = 10 (0.521), m = 14 (1.057) and
m = 19 (7.687) — and the last two exceed **one**. So the two envelopes have DIFFERENT scopes
here, and the row states both: the row-44 **half**-envelope holds for **m ≤ 9** (m = 10 is the
first mode to break it, at 0.5213), and the weaker **unit** contraction — each difference merely
smaller than the one two levels back — holds for **m ≤ 13** (m = 14 is the first to break it, at
1.0566). Gating either would be gating noise beyond those scopes, so
`two_step_ratio_max_per_mode`, `n_two_step_ratios_above_half` and
`modes_with_two_step_above_half` are registered in the golden and neither is asserted.

### Walltime is reported, never gated

The bottom panel is the paper's own cumulative walltime per mode. Published (Apple M2 Pro,
read off the PDF): 0.154 … 0.476 s at the first grid, 4.19 … 5.93 s at the last — "a few
seconds per m-mode" (main.tex:688). Ours, on a shared 256-core Linux box with one BLAS thread
per process, is roughly twice that at the last plotted grid, on meshes that are 1.78× SMALLER.
The driver holds ONE `mode_pool` for the whole figure and runs a throwaway warm-up sweep first,
because a fresh worker charges its entire compilation to its first mode's level 0 — without it
the bottom panel measures the compiler. The warm-up is itself the JIT payment and dominates the
run.


## Fig. 5 — what is plotted, which radii, and the two reasons panel 3 is preliminary

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

### The radii: the caption says ISCO, the panel titles say 9, 6, 2

| panel | published title | `kerr_isco(a)` | what we ran |
|---|---|---|---|
| a = −0.998 | `r₀ = 9` | 8.9943744548 | 8.9943744548 (the ISCO; 0.06 % from the title) |
| a = 0 | `r₀ = 6` | 6 | 6 (identical) |
| a = +0.998 | `r₀ = 2` | 1.2369706552 | 1.2369706552 (the ISCO; the title is 1.6× it) |

The published titles are corroborated by the published abscissae, each of which ends exactly at
twice its title's radius (18, 12, 4 — the caption's "clipped at `r = 2 r₀`"), so the published
panels were run at `r₀ = 9, 6, 2` and only the middle one is its spin's ISCO. This driver
follows the CAPTION and the node's brief — `r₀ = kerr_isco(a)`, never a typed radius — and every
panel title it draws carries the radius actually solved.

**`r₀ = 2 M` was not constructible when this figure landed — fixed upstream since.**
`CircularOrbit(a, 2.0).Lz` was `NaN`: Spacetime's generic circular closed form
(`capabilities/Spacetime/src/Geodesics/constants.jl:13`) divides by
`(p − 2) p − a²(x² − 1)`, which on the **equator** is exactly `p (p − 2)` and vanishes at
`p = 2` together with its numerator; `_carter_circular` (`constants.jl:17`) inherits it as
`0·NaN`, so **`Lz` and `Q`** were 0/0 while **`E`**, which carries no such factor, was already
right. `zc2 = 4(Lz² + …)` (main.tex:453) took the NaN and the puncture died in `legendre_P3`.
The NaN only bit where `r = 2 M` actually carries a circular geodesic, i.e. outside the
prograde photon orbit, `a > 1/√2 = 0.70710678`; below that the negative radicand makes the
`DomainError` those spins raise physically correct (at `a = 0`, `p = 2` is the horizon, where
the pre-fix code returned a silent `E = 0.0`). Measured: `NaN` at (0.9, 2.0) and (0.998, 2.0),
finite at 1.9, 2.5, 3.0. `ESSF::fix-horizon-axis-layer`'s sibling node removed it at `47f0ec0`
(the Bardeen–Press–Teukolsky equatorial pair behind an exact removable-endpoint trigger,
bitwise-neutral away from `p = 2`), with `7082f9a` pinning the orbit and the puncture at
`r₀ = 2 M` and `12b43bd` excluding `p = 0`; a conditioning band survives around the endpoint
(≈ `log₁₀|p − 2|⁻¹` digits lost, ~7 at `|p − 2| ≈ 1e-8`) and is that row's measured envelope.
So the published third panel IS constructible now — but regenerating it waits on the layer
below, which is a property of the inner boundary and not of `r₀`.

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

### [HOLE] a spurious horizon–axis boundary layer takes over the prograde ISCO panel

At `a = +0.998`, `r₀ = r_ISCO = 1.2370 M` the whole radial gap between horizon and orbit is
0.174 M, so the widest legal worldtube is 0.116 M and both vtu transitions are squeezed into
`[1.063, 1.55] M`. What breaks the m-decay there is **not** a failure to converge — it is a
spurious layer on the innermost radial node line, `r = r₊`, growing toward the axis:

* the mode variable has `sin^|m|θ` factored out, so axis regularity requires each plane to fall
  like `(1 − z²)^{m/2}` as `z² → 1`. MEASURED at m = 17 along `r = r₊` it RISES instead:
  `|u| = 1.6e-5` at `z² = 0.003` to `1.05` at `z² = 0.9918`, against the 0.98 → 1.8e-18 that
  regularity asks for — a **23-decade violation**;
* one node in from the boundary (`r > r₊ + 0.02`) the same modes decay physically, −0.107
  decades per mode over m = 13…20, and GMRES converges in 10–11 iterations. The solver and the
  interior are fine; the boundary row is not;
* the layer exists at every spin and is a question of SCALE. Inner-boundary axis:equator ratio
  at m = 17 / 20: **1.14 / 1.14** at `a = 0` (absent), **4.5e3 / 8.9e3** at `a = −0.998` but at
  an absolute 3.70e-8 — four decades below that panel's maximum, hence invisible — and
  **6.57e4 / 4.08e4** at O(1) for `a = +0.998` (`|u| = 1.05` at m = 17), where it dominates the
  plane from m ≈ 14 and is what `p3_decay_monotone_from_m2 = false` records. Every ratio here is
  re-measured from this figure's own surface cache (`evidence/fig-modes/layer_check_*.log`), and
  the "innermost node line" is the first Legendre-Gauss collocation line inside the
  horizon-adjacent element — `r − r₊` = 1.9e-4 M at this ISCO, 8.6e-3 and 4.4e-3 M on the other
  two panels.

Owned by **`ESSF::fix-horizon-axis-layer`** (the horizon boundary condition and axis
regularity). The paper solves this configuration ("prograde orbits at the ISCO … with spin
magnitude |a|/M = 0.998 … in a few seconds per m-mode"), so it is our limitation, and rows 46
and 53 inherit it — row 53's error budget needs this very configuration.

A SEPARATE, particle-side observation stands and is NOT the layer: at this ISCO
`|Ψ^ℛ(x_p)|` and `F_r^m` moved by 2-3× per refinement step over 8 → 10 → 12 grids, while the
same spin at `r₀ = 2.5 M` (worldtube 8× wider) decayed monotonically in m
(`|Ψ^ℛ(x_p)|` = 1.00e-2 → 2.19e-6 over m = 0 … 20). That control separates a compressed-domain
stall at the particle from the spin-driven layer at the boundary — but it was measured on the
PRE-ADOPTION multigrid hierarchy (at `5700c16`, before `fdcd919`), and its GMRES columns
(27/34/41 and 29/31/79 at the ISCO; 9/35/83/172/284 at `r₀ = 2.5`) are superseded: the landed
panel-3 range is **4–11**. The particle-side numbers have not been re-measured on the adopted
hierarchy.

### Gates are registered claims, not aborts

Both of the paper's claims about this figure are measured per panel and REGISTERED — the
fig_linsolv precedent — so `--check` keeps comparing numbers and a later change that breaks a
claim which currently HOLDS fails the comparison:

* `pK_decay_monotone_from_m2` — `max|Ψ_m|` strictly decreasing for m ≥ 2 (true, true, false),
  with `pK_decay_first_break_m` / `_last_break_m` and the second, particle-local instrument
  `pK_psi_xp_monotone_from_m2` (true, true, false);
* `pK_worldtube_decades` / `pK_outside_decades` — the in- and out-of-worldtube dynamic range of
  every plane, with `pK_m0_drops_two_decades` (2.43, 3.22, 2.65 decades) and
  `pK_m0_deeper_than_m1_m2` (true on all three panels: m = 0 drops 2.43 against 0.95 / 1.23,
  3.22 against 0.87 / 1.21, 2.65 against 1.09 / 0.93). NOTE the honest part: m = 0 is NOT the
  deepest plane of all — high-m planes drop further (`pK_deepest_worldtube_mode` = 18, 12, 19)
  because they oscillate more inside the tube. The paper's "particularly for the m = 0 mode" is
  reproduced in the form it holds — against the neighbouring highlighted modes — and the full
  vector is registered so the reading is checkable.

`gatefails` is reserved for structural faults: a lost mode, a violated precondition of the
extraction (`dom.chart === :z2` and `x_p[2] == 0`, which is what makes the sealed
raw-ingredient `mode_self_force(orbit, m, psi, dpsi_r)` bitwise equal to the `AmrHistory`
method's `:trace` reduction — verified identical at `a = 0, r₀ = 6, m = 2`), and an oversized
tracked PNG.

### The lattice these numbers belong to

Everything here is the production `ImposeEquatorialSymmetry: True` (yaml:36) `z²` half-domain —
320 dofs at level 0, 4535 at the 8th grid, 63 elements — registered as
`impose_equatorial_symmetry`, `pK_level0_ndofs` and `pK_final_ndofs`. ESSF::fig-convergence
found that the published fig:convergence abscissa implies the paper's own runs used the FULL-z
lattice (`ImposeEquatorialSymmetry: False`, 480 points at level 0); the physics is symmetric, so
the per-m planes are unaffected beyond resolution, but every dof/element count above is
lattice-dependent and should be read as this lattice's.

## Not yet reproduced

Figures 3 (`precond`) and 7 (`ssf_error_budget`) are later nodes of the
elliptic-ssf DAG. Nothing about them is claimed here.
