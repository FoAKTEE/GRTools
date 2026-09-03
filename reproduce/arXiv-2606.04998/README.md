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

## Results

| # | figure / check | driver | golden | tier | status | VLM verdict | notes |
|---|---|---|---|---|---|---|---|
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
off-pattern entries of each of those four diagonal blocks cancel. The structural count is
**9536** for every relative threshold from 1e-8 to 1e-16 — a nine-decade plateau, asserted by
the driver as `nnz_threshold_plateau` — and the figure uses `1e-12·max|A|`. The plateau stops
at 1e-8, not lower: 56 genuine structural entries of the outer wavezone blocks lie between
1e-8 and 1e-6 of the global maximum, since `max|A|` is set by the horizon column while `A^rr`
varies over decades across the compactified domain.

## Not yet reproduced

Figures 1 (`ssf_domain`), 3 (`precond`), 4 (`ssf_convergence`), 5 (`ssf_modes`),
6 (`linsolv`) and 7 (`ssf_error_budget`) are later nodes of the elliptic-ssf DAG; they need
the solver stack, the AMR driver and the mode-sum extraction. Nothing about them is claimed
here.
