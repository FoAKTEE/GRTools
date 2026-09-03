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
read the manual, render, VLM-check the render alone, VLM-compare against the reference PDF,
decide script-vs-structure, iterate (`.claude/skills/vlm-figure-repro/SKILL.md`).

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
| 2 | `fig:matrix` — DG operator sparsity at the coarsest grid level + element numbering | `fig_matrix.jl` | `expected/fig_matrix.json` | heavy | **solid** | **match** | Two panels: the B/W spy of the 320×320 coarsest-level matrix over a light-gray element-block checkerboard, and the domain sketch numbering the 20 elements that index its rows and columns. Structure reproduced exactly: 20 diagonal blocks + 33 symmetric off-diagonal block pairs, of which the 8 pairs across the four 2:1 interfaces (2↔3, 2↔4, 7↔9, 8↔9, 12↔13, 12↔14, 17↔19, 18↔19 in the paper's numbering) are fully dense 16×16 and the other 25 carry 64 entries. Adjacency graph-isomorphic to the domain topology in BOTH directions (0 matrix-only, 0 topology-only pairs). Reference block-occupancy extracted from `matrix.pdf` at 600 dpi and compared block-by-block: identical. 398.8 s at -t4 — **heavy**, although the numerics are ~10 s: the cost is the first-call specialization of the ComplexF64 `MeshOperator` (104 s) and of the first Plots/GR render. |
| 6 | `fig:linsolv` — GMRES residual staircase over AMR levels + walltime per iteration | `fig_linsolv.jl` | `expected/fig_linsolv.json` | heavy | **preliminary** | **partial** | Two panels for the m = 2 mode at a = 0.5 M, r₀ = 10 M: the residual staircase of 8 consecutive AMR levels over a cumulative-iteration axis, and the walltime per iteration. Presentation and structure reproduce item by item (see the checklist below); the per-level ITERATION COUNTS do not — published `1, 7, 6, 5, 4, 3, 3, 2` (31 cumulative) against our `1, 8, 9, 9, 9, 9, 9, 9` (63). All five node gates pass (max 9 < 10 iterations per level; level 0 in 1 iteration; the count bounded across levels; the first iteration slowest on every level; every level below the 1e-12 tolerance), so the paper's stated claims hold — but the shape of the published staircases does not, which is a solver-quality finding escalated to rows 21 and 43, and why this row is `preliminary`. 252.7 s at -t4 in `run_all` with a warm depot (~400 s standalone with precompilation), of which the AMR sweep itself is 12.3 s. Generated on top of `6df9f51` (the `hp_refine` 2:1-closure fix). |

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
**degrades** (16.4× → 6.0×), and our per-level initial residual (the interpolated warm start) is 1.4–13× larger
from level 2 on. Owner: `ESSF::solver-stack-core` (the GMRES ∘ multigrid ∘ Schwarz composition)
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

### The bottom panel: eager build, attributed by addition

The paper's first-iteration bump is the **lazy** construction of the Schwarz subdomain matrices
inside that iteration (main.tex:702). This package builds them **eagerly** in `rebuild!`, so the
same work is measured in `SolveRecord.build_time` and *not* inside `iteration_time[1]`. The bottom
panel therefore plots `iteration_time[1] + build_time` at each level's first point: the same
physical work charged to the same iteration, by addition rather than by measurement inside the
timer — a bookkeeping equivalence, not a lazy build. Both columns stay separately registered in
the golden. Measured consequences (envelope over the `-t1` and `-t4` runs): with the attribution the
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


## Not yet reproduced

Figures 1 (`ssf_domain`), 3 (`precond`), 4 (`ssf_convergence`), 5 (`ssf_modes`) and
7 (`ssf_error_budget`) are later nodes of the elliptic-ssf DAG; they need the mode-sum
extraction and the m-mode fan-out. Nothing about them is claimed here.
