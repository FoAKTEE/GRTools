# reproduce/arXiv-2606.04998/fig_matrix.jl — ESSF::fig-matrix.
#
# arXiv:2606.04998 Fig. 2 (`fig:matrix`, main.tex:622-626): "Matrix representation of the DG
# operator at the coarsest grid level (with equatorial symmetry). Bottom: numbering of DG
# elements in the domain corresponding to rows and columns in the matrix. Each element has 4×4
# grid points. Jump conditions are imposed at the worldtube boundary (red rectangle) and
# between the v, t, and u domains (red dotted lines)."
#
# THIN driver: zero physics. Everything is `ScalarSelfForce` /`EllipticCore` API —
# `crp_paper_config` → `block_domain` → `coarsest_operator` (= `domain_coefficients` +
# `jump_hooks` + `boundary_conditions` over `MeshOperator`, walked to the multigrid coarsest
# level) → `build_matrix`, with `vtu_junction_faces` and `worldtube_boundary_faces` supplying
# the two annotated face sets. The driver only arranges, gates and draws.
#
# PARAMETERS (an ASSUMPTION, stated: the caption names none). The production example of the
# paper's own input file is used verbatim — `anc/ScalarSelfForce.yaml:8-12`,
# `-p spin=0 -p orbital_radius=10 -p m_mode=2 -p inner_radius=2 -p outer_radius=1e4
#  -p worldtube_radius=4.8 -p worldtube_dz2=0.25 -p P=4 -p L=0`. The figure itself confirms
# the domain half of that choice: the caption's "each element has 4×4 grid points" is P = 4,
# and the drawn 20-element lattice IS the L = 0 `RefinedLevels` table of yaml:78-110. The
# matrix SPARSITY is m-independent (the mode number enters only the coefficient VALUES), so
# the figure's structure does not depend on the `m_mode = 2` half of the assumption; the
# registered magnitudes (`vtu_matrix_delta_maxabs`, `operator_on_zero_maxabs`) do.
#
# THE COARSEST GRID LEVEL. `LinearSolver.Multigrid.InitialCoarseLevels: L` (yaml:188) counts
# the COARSER levels only, so the hierarchy has `L + 1` levels and at the production `L = 0`
# the coarsest level IS the initial mesh: 20 elements × 16 dofs = 320. That is what the
# reference figure shows (its tick labels run 1…20) and what `coarsest_operator` returns at
# its yaml-twin default. Coarsening further, to the irreducible block structure, would give 10
# elements / 160 dofs — registered here as `nelements_block_structure` so the two readings of
# "coarsest" stay distinguishable.
#
# ELEMENT NUMBERING — a documented CONVENTION difference, gated in both directions. The paper
# numbers elements by block (the yaml flat block order, radial index fastest) and then, WITHIN
# a block, by the RADIAL segment index major and the ANGULAR index minor: its worldtube cell
# reads 3 (lower-left), 4 (upper-left), 5 (lower-right), 6 (upper-right). GRTools' `Mesh2D`
# enumerates the transpose (angular major), so its elements 4 and 5 (and 14, 15) are the
# paper's 5 and 4 (14, 15). Both panels here are drawn in the PAPER's order — a symmetric
# relabeling `P·A·Pᵀ` of the same operator, computed from the mesh's own `ElementId`s by the
# sort key `(block, radial segment, angular segment)`, never from the reference image. Every
# gate below is permutation-invariant and is asserted on the mesh-order data before the
# relabeling; the permutation itself is registered (`element_mesh_index`).
#
# THE SPARSITY THRESHOLD is measured, not chosen. Counting exact nonzeros gives 10112 entries,
# but 576 of them — the 144 off-pattern entries of each of the four COARSE-side diagonal
# blocks of the 2:1 interfaces — are roundoff: the mortar projections satisfy
# `Σ_mortars R·P = I` exactly, so those entries cancel analytically. The count is EXACTLY 9536
# for every relative threshold from 1e-8 down to 1e-16 (a NINE-decade plateau, asserted by the
# driver and registered as `nnz_threshold_plateau`), which is the structural gap; the figure
# and the block table use `1e-12·max|A|`. The plateau does NOT extend to 1e-6: 56 genuine
# structural entries of the outer wavezone blocks fall between 1e-8 and 1e-6 of the global
# maximum, because `max|A|` is set by the horizon column and `A^rr` varies over decades across
# the compactified domain.
#
#   julia --project=reproduce -t4 reproduce/arXiv-2606.04998/fig_matrix.jl --check
#
# Style: `_style/manuals/arXiv-2606.04998.md` §0 + §2 (B/W spy over a light-gray element
# checkerboard, red accents reserved for the worldtube and the vtu seams, serif labels). The
# bottom panel is a SCHEMATIC — the radial axis cannot be to scale (the outer column spans
# 24.4 M → 10⁴ M) — so the per-column drawing widths and the z² row split are presentation
# constants measured off the reference PDF; the cells they contain, their subdivision and
# their numbering all come from the mesh.

include(joinpath(@__DIR__, "..", "_lib", "repro_harness.jl"))
using .ReproHarness

opts = ReproHarness.parse_args()
ReproHarness.setup!()
ReproHarness.require_tier(opts, :fast)

using Printf, LinearAlgebra, SHA
using GrTools
using GrTools.ScalarSelfForce
using GrTools.EllipticCore

# ───────────────────────────── parameters (anc/ScalarSelfForce.yaml:8-12) ────────────────────
const SPIN = 0.0
const R0   = 10.0
const MM   = 2
const RIN  = 2.0            # = r₊ at a = 0 (block_domain requires the horizon)
const ROUT = 1.0e4
const WR   = 4.8
const WDZ2 = 0.25
const PORD = 4
const LREF = 0

for f in (:coarsest_operator, :assemble_operator)
    isdefined(ScalarSelfForce, f) || ReproHarness.hole(
        "ScalarSelfForce.$(f) is not available: Fig. 2 needs the operator-assembly node " *
        "(ESSF::fig-matrix commit 1, capabilities/ScalarSelfForce/src/Solve/assemble.jl).")
end

cfg = crp_paper_config(; spin = SPIN, orbital_radius = R0, m_mode = MM, inner_radius = RIN,
                       outer_radius = ROUT, worldtube_radius = WR, worldtube_dz2 = WDZ2,
                       P = PORD, L = LREF)
orb = CircularOrbit(SPIN, R0)
dom = block_domain(orb, cfg; orders = (PORD, PORD))
wtf = mesh -> worldtube_data(dom, orb, cfg, MM; mesh = mesh)

println("\n[fig_matrix] a = ", SPIN, "M  r₀ = ", R0, "M  m = ", MM, "  P = ", PORD,
        "  L = ", LREF, "   lattice ", dom.nr, "×", dom.nz, " blocks (", dom.chart, ")")
println("[fig_matrix] radial block bounds ", dom.bounds_r)
println("[fig_matrix] InitialCoarseLevels = ", cfg.multigrid.initial_coarse_levels,
        "  ⇒ multigrid levels = ", cfg.multigrid.initial_coarse_levels + 1)

lvl = coarsest_operator(ComplexF64, dom, orb, cfg, MM; worldtube = wtf)
mesh = lvl.mesh
nel = length(mesh.elements)
npe = prod(mesh.orders[1])
all(o -> prod(o) == npe, mesh.orders) ||
    error("fig_matrix: this driver assumes a p-uniform level mesh (orders $(unique(mesh.orders)))")
A = build_matrix(lvl.op)
amax = maximum(abs, A)
@printf("[fig_matrix] coarsest level: %d levels, %d elements × %d dofs = %d;  max|A| = %.6e\n",
        lvl.nlevels, nel, npe, mesh.ndofs, amax)

# ───────────────────────────── the structural sparsity threshold ─────────────────────────────
# measured, not chosen: the count is constant across eight decades (see the header)
const PLATEAU_TOLS = (1e-8, 1e-10, 1e-12, 1e-14, 1e-16)
plateau = [count(x -> abs(x) > t * amax, A) for t in PLATEAU_TOLS]
const TOL = 1e-12 * amax
nnz_struct = count(x -> abs(x) > TOL, A)
nnz_exact = count(!iszero, A)
@printf("[fig_matrix] nnz: exact %d, structural %d; plateau over %s → %s\n",
        nnz_exact, nnz_struct, string(PLATEAU_TOLS), string(plateau))
allequal(plateau) || error("fig_matrix: the sparsity threshold is NOT on a plateau ($plateau) — " *
                           "the structural pattern is not separated from roundoff")

# ───────────────────────────── the paper's element numbering ─────────────────────────────────
_segmid(s) = (s.index + 0.5) / (1 << s.level)
_seglo(s) = s.index / (1 << s.level)
_seghi(s) = (s.index + 1) / (1 << s.level)
# blocks in the yaml flat order, then radial segment MAJOR, angular segment MINOR
paperperm = sortperm(collect(1:nel); by = e -> (mesh.elements[e].block,
                                                _segmid(mesh.elements[e].segments[1]),
                                                _segmid(mesh.elements[e].segments[2])))
papernum = Vector{Int}(undef, nel)              # mesh index -> paper number
for (k, e) in enumerate(paperperm); papernum[e] = k; end

_rng(e) = (mesh.offsets[e] + 1):(mesh.offsets[e] + prod(mesh.orders[e]))
nnz_mesh = [count(x -> abs(x) > TOL, view(A, _rng(a), _rng(b))) for a in 1:nel, b in 1:nel]
nnz_paper = [nnz_mesh[paperperm[a], paperperm[b]] for a in 1:nel, b in 1:nel]

# ───────────────────────────── GATE 1: block adjacency ≡ domain topology ─────────────────────
# every nonzero off-diagonal element block ⇔ a shared face of the mesh, asserted BOTH ways
face_pairs = Set{Tuple{Int,Int}}()
for e in 1:nel, d in 1:2, s in 1:2
    ifc = mesh.faces[e][d][s]
    ifc === nothing && continue
    for nb in ifc.nbs
        push!(face_pairs, (min(e, nb), max(e, nb)))
    end
end
matrix_pairs = Set{Tuple{Int,Int}}()
for a in 1:nel, b in (a + 1):nel
    (nnz_mesh[a, b] > 0 || nnz_mesh[b, a] > 0) && push!(matrix_pairs, (a, b))
end
matrix_only = sort(collect(setdiff(matrix_pairs, face_pairs)))
topo_only = sort(collect(setdiff(face_pairs, matrix_pairs)))
isomorphic = isempty(matrix_only) && isempty(topo_only)
diag_ok = all(a -> nnz_mesh[a, a] > 0, 1:nel)
symmetric_pattern = all((nnz_mesh[a, b] > 0) == (nnz_mesh[b, a] > 0) for a in 1:nel, b in 1:nel)
println("[fig_matrix] GATE 1 sparsity ≅ topology: ", length(matrix_pairs),
        " off-diagonal block pairs vs ", length(face_pairs), " shared faces; matrix-only ",
        matrix_only, ", topology-only ", topo_only, " → ",
        isomorphic ? "ISOMORPHIC" : "MISMATCH")

# ───────────────────────────── the annotated face sets ───────────────────────────────────────
isin(el) = inside_worldtube(dom, el)
wtfaces = worldtube_boundary_faces(mesh, isin)
jfaces = vtu_junction_faces(dom, orb, cfg, MM; mesh = mesh)
wtels = sort(papernum[[mesh.index[el] for el in filter(isin, mesh.elements)]])

_facekey(f) = (papernum[mesh.index[f[1]]], Int(f[2]), Int(f[3]))
_facepair(f) = begin
    e = mesh.index[f[1]]
    nbs = mesh.faces[e][Int(f[2])][Int(f[3])].nbs
    [(min(papernum[e], papernum[nb]), max(papernum[e], papernum[nb])) for nb in nbs]
end
wt_keys = sort(_facekey.(wtfaces))
j_keys = sort(_facekey.(collect(keys(jfaces))))
wt_pairs = sort(unique(vcat([_facepair(f) for f in wtfaces]...)))
j_pairs = sort(unique(vcat([_facepair(f) for f in collect(keys(jfaces))]...)))
j_radii = sort(unique(Float64[j.r for j in values(jfaces)]))
println("[fig_matrix] worldtube: elements ", wtels, ", ", length(wtfaces),
        " boundary faces spanning element pairs ", wt_pairs)
println("[fig_matrix] vtu junctions: ", length(jfaces), " faces at r = ", j_radii,
        " spanning element pairs ", j_pairs)

# every worldtube / junction face is an ORDINARY shared face: neither creates a block that the
# domain topology does not already have (the assertion the caption's red marks must satisfy)
wt_pairs_mesh = Set((min(mesh.index[f[1]], nb), max(mesh.index[f[1]], nb))
                    for f in wtfaces
                    for nb in mesh.faces[mesh.index[f[1]]][Int(f[2])][Int(f[3])].nbs)
j_pairs_mesh = Set((min(mesh.index[f[1]], nb), max(mesh.index[f[1]], nb))
                   for f in keys(jfaces)
                   for nb in mesh.faces[mesh.index[f[1]]][Int(f[2])][Int(f[3])].nbs)
wt_extra = sort(collect(setdiff(wt_pairs_mesh, face_pairs)))
j_extra = sort(collect(setdiff(j_pairs_mesh, face_pairs)))

# ───────────────────────────── GATE 2: the row-40 oracle at this level ───────────────────────
# the vtu junction term IS in the matrix; the worldtube term is BITWISE absent from it and
# lives entirely in the fixed-source lift `operator_on_zero`
fcs = domain_coefficients(dom, orb, cfg)
bcs = boundary_conditions(dom, orb, cfg; mesh = mesh)
_mk(h) = MeshOperator(ComplexF64, mesh, fcs.A, fcs.beta, fcs.gamma;
                      C = cfg.dg.penalty_parameter, hooks = h, bcs = bcs)
_nojd = Dict{Tuple{ElementId,Int,Int},Vector{FaceState{ComplexF64}}}()
_nojf = Dict{Tuple{ElementId,Int,Int},VtuJunction{Float64}}()
A_novtu = build_matrix(_mk(jump_hooks(mesh, worldtube_jump_data(lvl.wt), _nojf)))
A_nowt = build_matrix(_mk(jump_hooks(mesh, _nojd, jfaces)))
d_vtu = maximum(abs, A .- A_novtu)
d_wt = maximum(abs, A .- A_nowt)
oz = maximum(abs, operator_on_zero(lvl.op))
println("[fig_matrix] GATE 2 row-40 oracle: |A−A_novtu|∞ = ", d_vtu,
        " (vtu IS in the matrix), |A−A_nowt|∞ = ", d_wt, " (bitwise ", A == A_nowt,
        "), |operator_on_zero|∞ = ", oz)

# ───────────────────────────── the block-structure reading of "coarsest" ─────────────────────
lvl_blk = coarsest_operator(ComplexF64, dom, orb, cfg, MM; worldtube = wtf, maxlevels = nothing)
@printf("[fig_matrix] fully coarsened (block structure): %d levels, %d elements, %d dofs\n",
        lvl_blk.nlevels, length(lvl_blk.mesh.elements), lvl_blk.mesh.ndofs)

# ───────────────────────────── registered values (the golden) ────────────────────────────────
pattern = UInt8[abs(A[i, j]) > TOL for i in 1:size(A, 1), j in 1:size(A, 2)]
pattern_paper = UInt8[abs(A[mesh.offsets[paperperm[fld(i - 1, npe) + 1]] + mod(i - 1, npe) + 1,
                           mesh.offsets[paperperm[fld(j - 1, npe) + 1]] + mod(j - 1, npe) + 1]) > TOL
                      for i in 1:mesh.ndofs, j in 1:mesh.ndofs]

@repro_result "spin" SPIN
@repro_result "orbital_radius" R0
@repro_result "m_mode" MM
@repro_result "grid_points_per_dim" PORD
@repro_result "initial_refinement_L" LREF
@repro_result "penalty_parameter" cfg.dg.penalty_parameter
@repro_result "block_bounds_r" collect(dom.bounds_r)
@repro_result "block_bounds_z2" collect(dom.bounds_z)
@repro_result "n_blocks" dom.nr * dom.nz
@repro_result "n_multigrid_levels" lvl.nlevels
@repro_result "n_elements" nel
@repro_result "dofs_per_element" npe
@repro_result "ndofs" mesh.ndofs
@repro_result "nelements_block_structure" length(lvl_blk.mesh.elements)
@repro_result "ndofs_block_structure" lvl_blk.mesh.ndofs
@repro_result "element_block" [mesh.elements[paperperm[k]].block for k in 1:nel]
@repro_result "element_mesh_index" paperperm
@repro_result "element_r_level" [mesh.elements[paperperm[k]].segments[1].level for k in 1:nel]
@repro_result "element_r_index" [mesh.elements[paperperm[k]].segments[1].index for k in 1:nel]
@repro_result "element_z_level" [mesh.elements[paperperm[k]].segments[2].level for k in 1:nel]
@repro_result "element_z_index" [mesh.elements[paperperm[k]].segments[2].index for k in 1:nel]
@repro_result "block_nnz" nnz_paper
@repro_result "nnz_exact" nnz_exact
@repro_result "nnz_structural" nnz_struct
@repro_result "nnz_threshold_plateau" plateau
@repro_result "matrix_maxabs" amax

offpairs = sort(collect(matrix_pairs))
@repro_result "n_offdiag_block_pairs" length(offpairs)
@repro_result "offdiag_pair_lo" [min(papernum[p[1]], papernum[p[2]]) for p in offpairs]
@repro_result "offdiag_pair_hi" [max(papernum[p[1]], papernum[p[2]]) for p in offpairs]
densepairs = sort([(min(papernum[a], papernum[b]), max(papernum[a], papernum[b]))
                   for (a, b) in offpairs if nnz_mesh[a, b] == npe^2])
@repro_result "n_dense_block_pairs" length(densepairs)
@repro_result "dense_pair_lo" [p[1] for p in densepairs]
@repro_result "dense_pair_hi" [p[2] for p in densepairs]
@repro_result "graph_isomorphic" isomorphic
@repro_result "n_matrix_only_pairs" length(matrix_only)
@repro_result "n_topology_only_pairs" length(topo_only)
@repro_result "diagonal_blocks_all_nonzero" diag_ok
@repro_result "sparsity_pattern_symmetric" symmetric_pattern

@repro_result "worldtube_elements" wtels
@repro_result "n_worldtube_faces" length(wtfaces)
@repro_result "worldtube_face_element" [k[1] for k in wt_keys]
@repro_result "worldtube_face_dim" [k[2] for k in wt_keys]
@repro_result "worldtube_face_side" [k[3] for k in wt_keys]
@repro_result "worldtube_pair_lo" [p[1] for p in wt_pairs]
@repro_result "worldtube_pair_hi" [p[2] for p in wt_pairs]
@repro_result "n_worldtube_pairs_outside_topology" length(wt_extra)
@repro_result "n_vtu_faces" length(jfaces)
@repro_result "vtu_face_element" [k[1] for k in j_keys]
@repro_result "vtu_face_dim" [k[2] for k in j_keys]
@repro_result "vtu_face_side" [k[3] for k in j_keys]
@repro_result "vtu_pair_lo" [p[1] for p in j_pairs]
@repro_result "vtu_pair_hi" [p[2] for p in j_pairs]
@repro_result "vtu_junction_radii" j_radii
@repro_result "n_vtu_pairs_outside_topology" length(j_extra)
@repro_result "vtu_matrix_delta_maxabs" d_vtu
@repro_result "worldtube_matrix_delta_maxabs" d_wt
@repro_result "operator_on_zero_maxabs" oz
@repro_result "sparsity_sha256_mesh_order" bytes2hex(sha256(collect(vec(pattern))))
@repro_result "sparsity_sha256_paper_order" bytes2hex(sha256(collect(vec(pattern_paper))))

# ───────────────────────────── the figure (presentation only) ────────────────────────────────
using Plots
include(joinpath(@__DIR__, "..", "_style", "emri_style.jl"))
emri_defaults!()
gr()

const RED = RGB(0.90, 0.09, 0.09)          # the paper's reserved worldtube / vtu accent
const CHECK = RGB(0.953, 0.953, 0.953)     # the element-block checkerboard (level 243/255)
# schematic drawing weights measured off the reference PDF (presentation only; the r axis
# cannot be to scale — the outer column spans 24.4 M → 10⁴ M — and the paper's z² axis is
# stretched the same way). The CELLS they contain and their numbering come from the mesh.
const COL_W_5 = (0.0651, 0.0781, 0.2511, 0.1881, 0.4175)
const ROW_H_2 = (0.318, 0.682)
const GAP_X = 0.012
const GAP_Y = 0.005
const SK_Y0, SK_Y1 = -0.62, -0.20          # the sketch band on the shared canvas

"the schematic column-slot edges (cumulative), and the row-band edges"
function _slots(dom)
    w = dom.nr == 5 ? collect(COL_W_5) : fill(1 / dom.nr, dom.nr)
    h = dom.nz == 2 ? collect(ROW_H_2) :
        [(dom.bounds_z[j+1] - dom.bounds_z[j]) / (dom.bounds_z[end] - dom.bounds_z[1])
         for j in 1:dom.nz]
    X = cumsum(vcat(0.0, w ./ sum(w)))
    Y = cumsum(vcat(0.0, h ./ sum(h)))
    return X, Y
end

"one element's schematic rectangle (x0, x1, y0, y1) on the shared canvas"
function _cell_rect(el, dom, X, Y, rcells)
    c = mod(el.block, dom.nr) + 1
    r = fld(el.block, dom.nr) + 1
    # the radial cells of this column, in order; the element's own is located among them
    cells = rcells[c]
    k = findfirst(t -> t[1] == _seglo(el.segments[1]) && t[2] == _seghi(el.segments[1]), cells)
    nc = length(cells)
    x0 = X[c] + (k - 1) * (X[c+1] - X[c]) / nc
    x1 = X[c] + k * (X[c+1] - X[c]) / nc
    yb0 = SK_Y0 + Y[r] * (SK_Y1 - SK_Y0)
    yb1 = SK_Y0 + Y[r+1] * (SK_Y1 - SK_Y0)
    y0 = yb0 + _seglo(el.segments[2]) * (yb1 - yb0)
    y1 = yb0 + _seghi(el.segments[2]) * (yb1 - yb0)
    return (x0 + GAP_X, x1 - GAP_X, y0 + GAP_Y, y1 - GAP_Y)
end

_rect(x0, x1, y0, y1) = Shape([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0])

function draw()
    plt = plot(size = (720, 1180), dpi = 150, legend = false, framestyle = :none,
               grid = false, ticks = false, aspect_ratio = :equal,
               xlims = (-0.115, 1.085), ylims = (-0.78, 1.115),
               left_margin = 0Plots.mm, right_margin = 0Plots.mm,
               top_margin = 0Plots.mm, bottom_margin = 0Plots.mm)

    # ---- TOP: the B/W spy plot over the element-block checkerboard --------------------------
    n = mesh.ndofs
    Z = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        a = fld(i - 1, npe) + 1
        b = fld(j - 1, npe) + 1
        Z[i, j] = (a + b) % 2 == 0 ? 1 : 0                    # checkerboard, gray on (a+b) even
    end
    for i in 1:n, j in 1:n
        pattern_paper[i, j] == 1 && (Z[i, j] = 2)
    end
    xs = ((1:n) .- 0.5) ./ n
    ys = 1 .- ((n:-1:1) .- 0.5) ./ n                          # row 1 at the TOP
    heatmap!(plt, xs, ys, Z[end:-1:1, :];
             c = cgrad([:white, CHECK, :black]; categorical = true), clims = (-0.5, 2.5),
             colorbar = false, linewidth = 0)
    plot!(plt, _rect(0, 1, 0, 1); linecolor = :black, linewidth = 1.1, fillalpha = 0)

    # element tick labels 1…N on the top and left edges, with "Ω_k" at the corner
    ctr(k) = (k - 0.5) / nel
    fs = nel <= 20 ? 9 : 6
    annotate!(plt, -0.052, 1.038, text("\$\\Omega_k\$", 10, :right, :bottom))
    for k in 1:nel
        annotate!(plt, ctr(k), 1.012, text("\$$(k)\$", fs, :center, :bottom))
        annotate!(plt, -0.014, 1 - ctr(k), text("\$$(k)\$", fs, :right, :vcenter))
    end
    # RED rectangles framing the worldtube elements' row band and column band (main.tex:625).
    # The band framing is only meaningful if the worldtube elements are CONTIGUOUS in the
    # paper's numbering (they are: one lattice block, numbered together).
    wlo, whi = minimum(wtels), maximum(wtels)
    wtels == collect(wlo:whi) || error("fig_matrix: the worldtube elements $wtels are not " *
        "contiguous in the paper's numbering — a row/column BAND cannot frame them")
    plot!(plt, _rect((wlo - 1) / nel, whi / nel, 1.0105, 1.0555);
          linecolor = RED, linewidth = 1.4, fillalpha = 0)
    plot!(plt, _rect(-0.0455, -0.0025, 1 - whi / nel, 1 - (wlo - 1) / nel);
          linecolor = RED, linewidth = 1.4, fillalpha = 0)
    # the "16 entries" scale bar: one element block wide
    plot!(plt, [0.0, npe / n], [-0.043, -0.043]; linecolor = :black, linewidth = 2.4)
    annotate!(plt, npe / n + 0.012, -0.043, text("$(npe) entries", 8, :left, :vcenter))

    # ---- BOTTOM: the domain sketch with the element numbering -------------------------------
    X, Y = _slots(dom)
    rcells = [sort(unique([(_seglo(el.segments[1]), _seghi(el.segments[1]))
                           for el in mesh.elements if mod(el.block, dom.nr) + 1 == c]))
              for c in 1:dom.nr]
    for k in 1:nel
        el = mesh.elements[paperperm[k]]
        x0, x1, y0, y1 = _cell_rect(el, dom, X, Y, rcells)
        plot!(plt, _rect(x0, x1, y0, y1); linecolor = :black, linewidth = 0.9, fillalpha = 0)
        annotate!(plt, (x0 + x1) / 2, (y0 + y1) / 2, text("\$$(k)\$", 10, :center, :vcenter))
    end
    # the worldtube: ONE red solid rectangle around the regularized cells (main.tex:625)
    wrs = [_cell_rect(mesh.elements[paperperm[k]], dom, X, Y, rcells) for k in wtels]
    plot!(plt, _rect(minimum(r[1] for r in wrs) - 0.010, maximum(r[2] for r in wrs) + 0.010,
                     minimum(r[3] for r in wrs) - 0.008, maximum(r[4] for r in wrs) + 0.008);
          linecolor = RED, linewidth = 2.0, fillalpha = 0)
    # the v|t and t|u seams: red DOTTED vertical lines at the transition radii
    for r in j_radii
        ib = findfirst(x -> x == r, dom.bounds_r)
        ib === nothing && continue
        plot!(plt, [X[ib], X[ib]], [SK_Y0 - 0.045, SK_Y1 + 0.055];
              linecolor = RED, linewidth = 1.6, linestyle = :dot)
    end
    # the v / t / u slicing spans
    ybr = SK_Y1 + 0.055
    iv = findfirst(x -> x == j_radii[1], dom.bounds_r)
    iu = findfirst(x -> x == j_radii[end], dom.bounds_r)
    plot!(plt, [X[iv], X[iu]], [ybr, ybr]; linecolor = :black, linewidth = 1.1)
    plot!(plt, [X[iv], X[iv] - 0.055], [ybr, ybr + 0.030]; linecolor = :black, linewidth = 1.1,
          arrow = true)
    plot!(plt, [X[iu], X[iu] + 0.055], [ybr, ybr + 0.030]; linecolor = :black, linewidth = 1.1,
          arrow = true)
    annotate!(plt, X[iv] + 0.010, ybr + 0.034, text("\$v\$", 11, :left, :bottom))
    annotate!(plt, (X[iv] + X[iu]) / 2, ybr + 0.012, text("\$t\$", 11, :center, :bottom))
    annotate!(plt, X[iu] + 0.010, ybr + 0.034, text("\$u\$", 11, :left, :bottom))
    # the axes of the sketch
    yax = SK_Y0 - 0.045
    plot!(plt, [-0.005, 1.045], [yax, yax]; linecolor = :black, linewidth = 1.0,
          arrow = true)
    plot!(plt, [-0.005, -0.005], [yax, SK_Y1 + 0.055]; linecolor = :black, linewidth = 1.0,
          arrow = true)
    annotate!(plt, 1.062, yax, text("\$r\$", 10, :left, :vcenter))
    annotate!(plt, -0.034, SK_Y1 + 0.030, text("\$z^2\$", 11, :center, :bottom))
    rticks = [(X[1], "\$r_+\$"), (X[2], "\$r_v\$"),
              ((X[3] + X[4]) / 2, "\$r_p\$"), (X[5], "\$r_u\$"),
              (X[end], "\$R = 10^4\$")]
    for (x, lab) in rticks
        plot!(plt, [x, x], [yax, yax + 0.014]; linecolor = :black, linewidth = 0.9)
        annotate!(plt, x, yax - 0.012, text(lab, 8, :center, :top))
    end
    for (j, lab) in zip(1:(dom.nz + 1), ("\$0\$", "\$0.25\$", "\$1\$"))
        yv = SK_Y0 + Y[j] * (SK_Y1 - SK_Y0)
        plot!(plt, [-0.005, 0.010], [yv, yv]; linecolor = :black, linewidth = 0.9)
        annotate!(plt, -0.016, yv, text(lab, 8, :right, :vcenter))
    end

    plot!(plt; xlims = (-0.115, 1.085), ylims = (-0.78, 1.115))   # heatmap must not widen them
    png = joinpath(@__DIR__, "expected", "fig_matrix.png")
    pdf = joinpath(@__DIR__, "expected", "fig_matrix.pdf")
    mkpath(dirname(png))
    savefig(plt, png)
    savefig(plt, pdf)
    nb = filesize(png)
    println("[fig_matrix] wrote ", relpath(png, ReproHarness.REPRO_ROOT), " (", nb, " bytes) and ",
            relpath(pdf, ReproHarness.REPRO_ROOT), " (", filesize(pdf), " bytes, untracked)")
    nb <= 200 * 1024 || error("fig_matrix.png is $(nb) bytes > 200 KB — reduce the dpi/size")
    return png
end

draw()

# ───────────────────────────── the node's own gates ──────────────────────────────────────────
fails = String[]
isomorphic || push!(fails, "GATE 1: the element-block adjacency graph is NOT isomorphic to " *
                           "the domain topology (matrix-only $matrix_only, topology-only $topo_only)")
diag_ok || push!(fails, "GATE 1: some diagonal element block is empty")
isempty(wt_extra) || push!(fails, "GATE 1: worldtube faces create block(s) the topology does " *
                                  "not have: $wt_extra")
isempty(j_extra) || push!(fails, "GATE 1: vtu junction faces create block(s) the topology does " *
                                 "not have: $j_extra")
d_vtu > 1e-6 || push!(fails, "GATE 2: the vtu junction term is ABSENT from build_matrix " *
                             "(|A−A_novtu|∞ = $d_vtu)")
A == A_nowt || push!(fails, "GATE 2: the worldtube term is PRESENT in build_matrix " *
                            "(|A−A_nowt|∞ = $d_wt), but it must live only in operator_on_zero")
oz > 1e-3 || push!(fails, "GATE 2: operator_on_zero carries no worldtube lift ($oz)")
@repro_result "gates_pass" isempty(fails)
if !isempty(fails)
    println("\n[fig_matrix] FINDING — node gates violated:")
    foreach(f -> println("    ", f), fails)
end

code = ReproHarness.conclude!(opts; paper = "2606.04998", id = "fig_matrix",
                              script = @__FILE__, rtol = 1e-10,
                              tolerances = Dict("vtu_matrix_delta_maxabs" => (rtol = 1e-9, atol = 0.0),
                                                "operator_on_zero_maxabs" => (rtol = 1e-9, atol = 0.0),
                                                "matrix_maxabs" => (rtol = 1e-9, atol = 0.0)))
exit(isempty(fails) ? code : max(code, 1))
