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
> not pin the value, and each driver states which of its registered numbers depend on it.

> **Why `InitialCoarseLevels: L` means "no coarser levels" at `L = 0`.** Primary source:
> spectre `src/ParallelAlgorithms/LinearSolver/Multigrid/Tags.hpp:30-37` — *"Maximum number of
> coarser levels in the initial multigrid hierarchy. Set to `'0'` to create no coarser levels
> initially. Set to `'Auto'` to coarsen all the way up to single-element blocks."*

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
| 1 | `fig:domain` — the DG computational domain: `\|Ψ_m\|` height- and colour-mapped over the (r, z²) element grid, `\|Ψ^ℛ_m\|` in the worldtube | `fig_domain.jl` | `expected/fig_domain.json` | heavy | **solid** | **match** | A single 3D panel plus an inset: |Ψ_m| height- **and** colour-mapped over the (r, z²) DG element grid with white element edges and the Legendre-Gauss collocation grid, |Ψ^ℛ_m| inside the red-framed worldtube, the v/t/u slicing arc with H = ∓1, and the six-step hp-AMR cascade at x_p magnified bottom-right. Seven AMR levels, 20 → 57 elements and 320 → 3538 dofs, GMRES 1, 8, 9, 8, 9, 9, 9. Structure reproduced: the two |Ψ_m| lobes of the published render — the tan ridge in the polar row over the worldtube's radial span and the sharper wedge just outside the tube in r at the equator — agree in count and position, as do the flat navy near-horizon band, the low navy tube patch and the teal wavezone plateau. `--check` 52/52 at rtol 1e-10 with max|Δ| = 0, `-t1` ≡ `-t4`. ~500 s at -t4, of which the sweep is 223.5 s — **heavy**. |
| 2 | `fig:matrix` — DG operator sparsity at the coarsest grid level + element numbering | `fig_matrix.jl` | `expected/fig_matrix.json` | heavy | **solid** | **match** | Two panels: the B/W spy of the 320×320 coarsest-level matrix over a light-gray element-block checkerboard, and the domain sketch numbering the 20 elements that index its rows and columns. Structure reproduced exactly: 20 diagonal blocks + 33 symmetric off-diagonal block pairs, of which the 8 pairs across the four 2:1 interfaces (2↔3, 2↔4, 7↔9, 8↔9, 12↔13, 12↔14, 17↔19, 18↔19 in the paper's numbering) are fully dense 16×16 and the other 25 carry 64 entries. Adjacency graph-isomorphic to the domain topology in BOTH directions (0 matrix-only, 0 topology-only pairs). Reference block-occupancy extracted from `matrix.pdf` at 600 dpi and compared block-by-block: identical. 398.8 s at -t4 — **heavy**, although the numerics are ~10 s: the cost is the first-call specialization of the ComplexF64 `MeshOperator` (104 s) and of the first Plots/GR render. |

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

    r_c(r) = 97.779 − 1790.6/r   on  [24.4, 10⁴] → [24.4, 97.6],     r_c(r) = r   below r_in.

That is this figure's horizontal coordinate, taken from `cfg.observers.rcc_*` rather than
hard-coded. **Framing:** the published render itself stops at the u-block element boundary
`r = 48.68 M`; measured on `ssf_domain.pdf` the published worldtube patch sits at 4.5–21.5 % of
the drawn ribbon, which the compressed coordinate predicts as 5.6–21.9 % if the ribbon is cut
there (against 3.5–13.5 % for the full `r_c ∈ [1.87, 97.6]`). The driver frames the same way
(`R_FRAME = 52`) — the outermost wavezone element runs off the right edge exactly as it does in
the paper. Nothing is removed from the solve or from the golden.

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
4. **VLM** against `ssf_domain.pdf`: **match**. All four structures the manual demands (warped coloured surface with white DG lines, flat far zone, detached red-framed patch, AMR inset) and 9/9 checklist items; the lobe count and placement agree. Five style differences are recorded in `evidence/fig-domain/vlm_verdict.md` — proportions (the published ParaView camera is not recoverable from the PDF), ParaView-vs-Lambert lighting, the added `H = ±1` labels, the drawn cliff face at the tube boundary, and the inset window.

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

## Not yet reproduced

Figures 3 (`precond`), 4 (`ssf_convergence`), 5 (`ssf_modes`), 6 (`linsolv`) and
7 (`ssf_error_budget`) are later nodes of the elliptic-ssf DAG; they need
the solver stack, the AMR driver and the mode-sum extraction. Nothing about them is claimed
here.
