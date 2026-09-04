# test_multigrid.jl — grounded gates for src/multigrid.jl (DAG node: multigrid-hierarchy,
# elliptic-ssf GLOBAL_DAG row 19 + TOPOLOGICAL_SCHEDULE W7-A, CRITICAL PATH). Written
# BEFORE the implementation (TDD red run captured in evidence/multigrid-hierarchy/tests_red.log).
#
# Target: the CRP main.tex:670-676 preconditioner level — "Each GMRES iteration is
# preconditioned with one multigrid V-cycle, where the multigrid hierarchy is built by
# h-coarsening the domain up to the initial level (the irreducible block structure)",
# with an additive Schwarz smoother (2 passes) on every finer level and a direct bottom
# solve on the coarsest. The yaml twin (anc/ScalarSelfForce.yaml:186-194, mirrored by
# ScalarSelfForce `MultigridConfig`): Iterations 1, PreSmoothing True, UseBottomSolver
# True, PostSmoothingAtBottom False, ComplexShift 0.0.
#
# Model problems (coefficient callables only — zero physics; the m-mode integration is
# the solver-stack-ssf node's scope):
#   * DEFINITE complex model  A = −∇·∇ + β,  β = 1 + i/2 on a TWO-BLOCK domain
#     [0,2]×[0,1] (two affine unit blocks joined by a BlockConnection), zero-Dirichlet
#     ghosts, N = 4×4 Gauss points per element, C = 1.5. Uniform h-refinement level k
#     gives exactly k+1 multigrid levels (k = 2, 3, 4 → 3, 4, 5 levels).
#   * MILDLY INDEFINITE complex Helmholtz β = −1 + 0.3i on the same domain.
#   * NONCONFORMING meshes: the 5-element 2:1 + p-jump mesh of test_schwarz.jl and a
#     k = 2 two-block mesh with one element split once more (2:1 h interface) and one
#     column at order 5 (p-jump).
#   * CURVED model: one InverseRadialMap1D(2, 20) block (the CRP eq:compaction wavezone
#     map) — the scope boundary of the mass-weighted transfer identities (gate 4).
#
# Oracle gates (constants calibrated PRE-RED in evidence/multigrid-hierarchy/
# calibrate_multigrid{,2,3}.jl + calibration.log — THROWAWAY prototypes composing the
# sealed W5/W6 APIs, NOT the implementation):
#   (1) COARSENING (REPAIRED — the shipped per-dimension level CAP was falsified on the
#       CRP's own AMR class): one step coarsens every COMPLETE SIBLING GROUP by one level,
#       largest dimension mask first and each group consumed as a whole (letting an element
#       pick its own mask independently overlaps targets: with A, B, C present and D
#       refined, A picks {1} → A|B while C picks {2} → A|C), then the 2:1 balance is
#       restored with the sealed anisotropic `balance_two_to_one`. Gates: the validator's
#       18-element anisotropic counterexample builds a full hierarchy [18, 9, 5, 3, 2, 1]
#       (the cap rule threw "invalid face coverage"); a 300-mesh census over AMR-reachable
#       ANISOTROPIC meshes (base level 2, 3 rounds of h_refine with per-dimension flags —
#       anc/ScalarSelfForce.yaml:145 Isotropy: Anisotropic, Iterations 8) builds 300/300
#       valid, balanced, nested, exactly-tiling hierarchies where the cap rule managed
#       56/300 (and 4/150 with a block seam); uniform meshes shrink geometrically,
#       64 → 16 → 4 → 1 (spectre Hierarchy.cpp:20-46 parity); orders p-coarsen by the MIN
#       rule (spectre Multigrid.hpp "we choose the smaller of the two polynomial degrees").
#   (2) TRANSFER TRIO (`mg_transfer_check`), on the affine, the 2:1+p, and the CURVED mesh:
#       prolongation EXACT on the coarse element's polynomial space (calibrated 8.0e-16 …
#       1.6e-15), R∘P = I (calibrated 8.7e-16 … 1.1e-15), and adjointness
#       ⟨R v, u⟩_c = ⟨v, P u⟩_f in the mass-weighted inner product (calibrated ≤ 4.9e-16) —
#       all gated < 1e-12. At BigFloat(256) the same three are calibrated ≤ 1.1e-76,
#       gated < 1e-70. FINDING (round 2A): "exact on polynomials" must be read in the
#       element's LOGICAL coordinate — a PHYSICAL-coordinate polynomial is not in the
#       coarse space when the block map is not affine (measured 1.16e+03 on the 1/r block).
#   (3) V-CYCLE CONTRACTION < 0.2, INDEPENDENT of the level count: calibrated 0.1715 /
#       0.1718 / 0.1719 at 3 / 4 / 5 levels (spread 4e-4, gated < 0.01) on the definite
#       model, 0.1719 on the mildly indefinite one, 0.1667 on the 2:1+p mesh. Off-gates
#       kept as measured scope: the 5-element mixed mesh contracts at 0.1971 (gated
#       < 0.25 — nearly all smoother, no coarse space to speak of) and the compactified
#       1/r block at 0.3049 (gated < 0.40 — the geometric coarse grid is a poorer
#       approximation on a strongly graded map). The 2 Schwarz passes are load-bearing:
#       1 pass measures 0.3064, 3 passes 0.1049 (round 3D).
#   (4) MASS-WEIGHTED / GALERKIN identities and their SCOPE. On AFFINE blocks the physical
#       mass adjointness M_c·R == Pᵀ·M_f holds to 4.9e-16 and the massive restriction
#       equals M_c ∘ R ∘ M_f⁻¹ to 3.2e-16 (gated < 1e-12). On the CURVED block BOTH fail at
#       O(1) (measured 1.08 and 1.69 relative — the diagonal Gauss mass integrates the
#       polynomial part exactly but not a non-constant Jacobian); the test PINS that
#       failure (> 0.1) so the affine-only scope cannot silently rot. The V-cycle uses the
#       GALERKIN restriction Pᵀ (spectre `operand_is_massive` branch of
#       Spectral::projection_matrix_child_to_parent), which is exact for ANY map.
#   (5) FALSIFIERS: a flipped prolongation sign breaks the trio (measured 3.75 / 2.28) and
#       the contraction (0.3354 > 0.2); dropping the half-interval Jacobian of the
#       restriction breaks R∘P = I and the adjointness (1.32 / 1.00) while leaving the
#       V-cycle untouched — R is NOT on the V-cycle path; and restricting the MASSIVE
#       residual with the nodal L2 restriction instead of Pᵀ degrades the contraction from
#       0.1714 to 0.6756 (gated > 0.5, the Galerkin-consistency falsifier).
#   (6) GMRES ∘ MultigridV vs GMRES ∘ Schwarz on the manufactured discrete solution
#       b := A·u_ex at reltol 1e-10: calibrated 7 vs 32 iterations at 4 levels (k = 3) and
#       7 vs 18 at 3 levels (k = 2) — the multigrid count is LEVEL-INDEPENDENT while the
#       Schwarz-only count nearly doubles; unpreconditioned GMRES needs 213 (k = 2) and
#       does not converge in 400 (k = 3).
#   (7) BOTTOM SEAM: `dense_bottom_solver` solves the coarsest system to a relative
#       residual ≤ 1.3e-15 (gated < 1e-12); a caller-supplied `(z, b)` callable is honored
#       (call counter: exactly one bottom solve per V-cycle) and `bottom = nothing` falls
#       back to the smoother (the yaml UseBottomSolver: False twin) — measured to contract
#       IDENTICALLY (0.1714) on a 2-element coarsest grid, where two overlapping Schwarz
#       passes already solve the bottom system exactly. This is the API seam the
#       multigrid-bottom-solver node (cached sparse UMFPACK LU) plugs into.
#   (8) CACHE + ALLOCATIONS: the per-level Schwarz LU caches build exactly one
#       factorization per subdomain at hierarchy construction (calibrated 170 across the
#       4-level hierarchy) and build ZERO more across V-cycles; the preconditioner apply
#       allocates a CONSTANT number of bytes per call after warm-up (the sealed
#       `smooth!`/`apply!` layers allocate their own workspaces per call — documented, not
#       gated to zero).
#   (9) BigFloat(256) V-cycle: contraction 0.1712 (matching the Float64 0.1715) and
#       GMRES ∘ MultigridV converging to 4.2e-32 in 9 iterations (round 3A).
#  (10) THREAD PARITY: `vcycle!(...; threaded = true)` is BITWISE `threaded = false`
#       (inherited from the bitwise-gated threaded `apply!`/`smooth!`); calibrated as an
#       identical iterate hash at -t1 and -t4 (round 3B).
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting siblings).

using Test
using LinearAlgebra
using REPL          # Julia >=1.11: @doc is a no-op unless REPL is loaded (runtests.jl convention)
using EllipticCore

# ---- calibrated gate constants (evidence/multigrid-hierarchy/calibration.log) ----------
const MG_TRIO_TOL = 1e-12          # calibrated ≤ 1.6e-15 (Float64)
const MG_TRIO_TOL_BIG = 1e-70      # calibrated ≤ 1.1e-76 (BigFloat 256)
const MG_RHO_GATE = 0.2            # calibrated 0.1715 / 0.1718 / 0.1719 at 3 / 4 / 5 levels
const MG_RHO_SPREAD = 0.01         # calibrated spread 4.0e-4 across the three level counts
const MG_RHO_TINY = 0.25           # 5-element mixed mesh, calibrated 0.1971
const MG_RHO_CURVED = 0.40         # compactified 1/r block, calibrated 0.3044
const MG_RHO_NC = 0.25             # 2:1+p two-block mesh, calibrated 0.1983 (repair)
const MG_RHO_CURVED_NC = 0.45      # curved 1/r WITH a 2:1 interface and a p-jump: 0.3707
const MG_RHO_ANISO = 0.2           # anisotropically h-refined mesh, calibrated 0.1679
const MG_CENSUS_N = 300            # anisotropic AMR census size (repair round 1: 300/300)
const MG_COST_RATIO = 1.8          # sum(ndofs)/ndofs_fine, graded mesh at shrink = 0.5: 1.72
const MG_RHO_SIGN = 0.2            # sign falsifier: calibrated 0.3354 (must EXCEED the gate)
const MG_RHO_NONMASSIVE = 0.5      # Galerkin falsifier: calibrated 0.6756 (must EXCEED)
const MG_GMRES_GATE = 12           # calibrated 7 iterations at 3 AND 4 levels
const MG_GMRES_SPREAD = 2          # level-independence of the GMRES count
const MG_BOTTOM_TOL = 1e-12        # calibrated 1.3e-15 relative
const MG_CURVED_SCOPE = 0.1        # the curved-block mass identities must FAIL above this

# ---- tiny deterministic generator (Knuth-MMIX LCG; distinct name/seed per test file) ---
mutable struct MGRNG
    s::UInt64
end
mg_unit!(r::MGRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
mg_cvec!(r::MGRNG, n) = [complex(mg_unit!(r), mg_unit!(r)) for _ in 1:n]

const MG_C = 1.5                                   # CRP main.tex:595 penalty parameter
mg_A11(x, y) = 1.0
mg_A22(x, y) = 1.0
mg_zero2(x, y) = 0.0
const MG_BETA_DEF = 1.0 + 0.5im                    # definite complex model
const MG_BETA_MILD = -1.0 + 0.3im                  # mildly indefinite complex Helmholtz

function mg_zerobcs(mesh, T)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(T, nf))
    end
    return bcs
end

# operator factory: the SAME coefficient closures rediscretized on whatever mesh it is
# handed — exactly what MultigridHierarchy needs for the coarse levels
mg_factory(T, β) =
    mesh -> MeshOperator(T, mesh, (mg_A11, mg_A22), (x, y) -> β, (mg_zero2, mg_zero2);
                         C = MG_C, bcs = mg_zerobcs(mesh, T))

# two affine unit blocks joined along dimension 1: [0,1]×[0,1] and [1,2]×[0,1]
function mg_twoblocks(RT)
    b1 = Block(1, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1)))
    b2 = Block(2, LinearMap1D(RT(1), RT(2)), LinearMap1D(RT(0), RT(1)))
    return [b1, b2], [BlockConnection(1, 1, 2)]
end

function mg_uniform_mesh(RT, k, N)
    blks, conns = mg_twoblocks(RT)
    els = ElementId[]
    for b in (1, 2), iy in 0:(2^k-1), ix in 0:(2^k-1)
        push!(els, ElementId(b, SegmentId(k, ix), SegmentId(k, iy)))
    end
    return Mesh2D(blks, els, (N, N); connections = conns)
end

# the 5-element 2:1 + p-jump mesh (test_schwarz.jl shape, own block id)
function mg_mixed_mesh(RT)
    blk = Block(7, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1)))
    W = ElementId(7, SegmentId(1, 0), SegmentId(0, 0))
    E11 = ElementId(7, SegmentId(2, 2), SegmentId(1, 0))
    E12 = ElementId(7, SegmentId(2, 2), SegmentId(1, 1))
    E21 = ElementId(7, SegmentId(2, 3), SegmentId(1, 0))
    E22 = ElementId(7, SegmentId(2, 3), SegmentId(1, 1))
    els = [W, E11, E12, E21, E22]
    ords = Dict(W => (4, 4), E11 => (4, 4), E12 => (4, 4), E21 => (5, 5), E22 => (5, 5))
    return Mesh2D(blk, els, ords)
end

# k = 2 two-block mesh, one element split once more (2:1 h) and one column at order 5 (p)
function mg_nonconf_mesh(RT)
    blks, conns = mg_twoblocks(RT)
    els = ElementId[]
    ords = Dict{ElementId,NTuple{2,Int}}()
    for b in (1, 2), iy in 0:3, ix in 0:3
        if b == 1 && ix == 1 && iy == 1
            for jy in 0:1, jx in 0:1
                c = ElementId(b, SegmentId(3, 2 + jx), SegmentId(3, 2 + jy))
                push!(els, c)
                ords[c] = (4, 4)
            end
        else
            el = ElementId(b, SegmentId(2, ix), SegmentId(2, iy))
            push!(els, el)
            ords[el] = (b == 2 && ix == 3) ? (5, 5) : (4, 4)
        end
    end
    return Mesh2D(blks, els, ords; connections = conns)
end

# the CRP eq:compaction wavezone map (main.tex:546-549) as a single curved block
function mg_curved_mesh(RT, k)
    blk = Block(9, InverseRadialMap1D(RT(2), RT(20)), LinearMap1D(RT(0), RT(1)))
    els = [ElementId(9, SegmentId(k, ix), SegmentId(k, iy))
           for iy in 0:(2^k-1) for ix in 0:(2^k-1)]
    return Mesh2D(blk, els, (4, 4))
end

mg_erange(mesh, e) = (mesh.offsets[e] + 1):(mesh.offsets[e] + prod(mesh.orders[e]))

# stationary V-cycle iteration on A·u = 0 from a random start: the error IS the iterate,
# so successive norm ratios are the contraction factors
function mg_contraction(mg, n; niter = 8, seed = 0xabc)
    rng = MGRNG(seed)
    x = mg_cvec!(rng, n)
    b = zeros(ComplexF64, n)
    ns = Float64[norm(x)]
    for _ in 1:niter
        vcycle!(x, mg, b)
        push!(ns, norm(x))
    end
    return [ns[i+1] / ns[i] for i in 1:length(ns)-1]
end

# INDEPENDENT reference V-cycle over the public transfer API (no MultigridV): the same
# algorithm written a second way, and the vehicle for the Galerkin falsifier (`massive`).
function mg_ref_vcycle!(x, hier, b, ℓ, passes, bottomfac; massive = true)
    L = length(hier.levels)
    lev = hier.levels[ℓ]
    if ℓ == L
        x .= bottomfac \ Vector(b)
        return x
    end
    smooth!(x, lev.smoother, b; passes = passes)
    y = similar(x)
    apply!(y, lev.op, x)
    r = b .- y
    nc = hier.levels[ℓ+1].mesh.ndofs
    bc = zeros(eltype(x), nc)
    massive ? mg_restrict_residual!(bc, hier, ℓ, r) : mg_restrict_field!(bc, hier, ℓ, r)
    xc = zeros(eltype(x), nc)
    mg_ref_vcycle!(xc, hier, bc, ℓ + 1, passes, bottomfac; massive = massive)
    mg_prolong!(y, hier, ℓ, xc)
    x .+= y
    smooth!(x, lev.smoother, b; passes = passes)
    return x
end

function mg_ref_contraction(hier, passes; niter = 8, seed = 0xabc, massive = true)
    n = hier.levels[1].mesh.ndofs
    fac = lu(build_matrix(hier.levels[end].op))
    rng = MGRNG(seed)
    x = mg_cvec!(rng, n)
    b = zeros(ComplexF64, n)
    ns = Float64[norm(x)]
    for _ in 1:niter
        mg_ref_vcycle!(x, hier, b, 1, passes, fac; massive = massive)
        push!(ns, norm(x))
    end
    return [ns[i+1] / ns[i] for i in 1:length(ns)-1]
end

# per-element physical mass diagonal, rebuilt INDEPENDENTLY through the exported
# element_matrices (not read off the operator) — the oracle for the mass-weighted gates
function mg_massdiag(T, mesh)
    out = Vector{real(float(T))}(undef, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        out[mg_erange(mesh, e)] .= diag(element_matrices(T, blk, el, mesh.orders[e]).M)
    end
    return out
end


# ---- repair-node helpers (H1/H2/H3/M3-M5) ----------------------------------------------
# the validator's 18-element ANISOTROPIC counterexample: four x-columns with x-levels
# 3/3/2/1 and y-levels 1/2/3/2 — two_to_one_balanced, a valid Mesh2D, and fatal to the
# shipped per-dimension level cap (it produced nested, non-tiling ancestors).
function mg_counterexample()
    blk = Block(21, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    els = ElementId[]
    for (sx, ly) in ((SegmentId(3, 0), 1), (SegmentId(3, 1), 2),
                     (SegmentId(2, 1), 3), (SegmentId(1, 1), 2)), j in 0:(2^ly - 1)
        push!(els, ElementId(21, sx, SegmentId(ly, j)))
    end
    return Mesh2D(blk, els, (4, 4))
end

mutable struct MGCRNG
    s::UInt64
end
mg_c01!(r::MGCRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      (r.s >> 33) / 2.0^31)

# an AMR-REACHABLE anisotropic mesh: uniform base level, then `rounds` of h_refine with
# per-dimension NTuple{2,Bool} flags (anc/ScalarSelfForce.yaml:145 Isotropy: Anisotropic)
function mg_aniso_mesh(seed, rounds; base = 2, nblocks = 1)
    r = MGCRNG(UInt64(seed) * 0x9e3779b97f4a7c15 + 1)
    blks = [Block(b, LinearMap1D(Float64(b - 1), Float64(b)), LinearMap1D(0.0, 1.0))
            for b in 1:nblocks]
    conns = nblocks == 2 ? [BlockConnection(1, 1, 2)] : BlockConnection[]
    els = [ElementId(b, SegmentId(base, ix), SegmentId(base, iy))
           for b in 1:nblocks for iy in 0:(2^base-1) for ix in 0:(2^base-1)]
    mesh = Mesh2D(blks, els, (4, 4); connections = conns)
    for _ in 1:rounds
        flags = Dict{ElementId,NTuple{2,Bool}}()
        for el in mesh.elements
            a = mg_c01!(r) < 0.25
            b = mg_c01!(r) < 0.25
            (a || b) && (flags[el] = (a, b))
        end
        isempty(flags) || (mesh = h_refine(mesh, flags))
    end
    return mesh
end

# the CRP puncture cascade: uniform base, then `rounds` of isotropic corner refinement
function mg_graded_mesh(base, rounds)
    blk = Block(23, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    els = [ElementId(23, SegmentId(base, ix), SegmentId(base, iy))
           for iy in 0:(2^base-1) for ix in 0:(2^base-1)]
    mesh = Mesh2D(blk, els, (4, 4))
    for _ in 1:rounds
        corner = [el for el in mesh.elements
                  if el.segments[1].index == 2^el.segments[1].level - 1 &&
                     el.segments[2].index == 2^el.segments[2].level - 1]
        mesh = h_refine(mesh, Set(corner))
    end
    return mesh
end

# a curved (1/r compactified) block carrying BOTH a 2:1 interface and a p-jump
function mg_curved_nonconf()
    blk = Block(9, InverseRadialMap1D(2.0, 20.0), LinearMap1D(0.0, 1.0))
    els = [ElementId(9, SegmentId(2, ix), SegmentId(2, iy)) for iy in 0:3 for ix in 0:3]
    m1 = h_refine(Mesh2D(blk, els, (4, 4)),
                  Set([ElementId(9, SegmentId(2, 1), SegmentId(2, 1))]))
    ords = Dict(el => (el.segments[1].index == 2^el.segments[1].level - 1 ? (5, 5) : (4, 4))
                for el in m1.elements)
    return Mesh2D(m1.blocks, m1.elements, ords; connections = m1.connections)
end

# one full coarsening chain, validated at every level: valid Mesh2D, 2:1 balanced, every
# fine element NESTED in its parent (gap 0 or 1 per dimension) and every coarse element
# EXACTLY tiled by its fine descendants. Returns (element counts, sum(ndofs)/ndofs_fine).
function mg_chain_report(mesh)
    nels = [length(mesh.elements)]
    ndofs = [mesh.ndofs]
    cur = mesh
    while true
        co = coarsen_mesh(cur)
        co === nothing && break
        cm, parent = co
        two_to_one_balanced(cm.elements) || error("unbalanced coarse level")
        area(el) = ldexp(1.0, -el.segments[1].level - el.segments[2].level)
        acc = zeros(Float64, length(cm.elements))
        for (e, el) in enumerate(cur.elements)
            c = cm.elements[parent[e]]
            c.block == el.block || error("parent in another block")
            for d in 1:2
                Δ = el.segments[d].level - c.segments[d].level
                (0 <= Δ <= 1 && (el.segments[d].index >> Δ) == c.segments[d].index) ||
                    error("not nested: $el -> $c")
            end
            acc[parent[e]] += area(el)
        end
        for k in eachindex(cm.elements)
            abs(acc[k] - area(cm.elements[k])) < 1e-12 || error("coarse element $k not tiled")
        end
        cur = cm
        push!(nels, length(cm.elements))
        push!(ndofs, cm.ndofs)
        length(nels) > 40 && error("coarsening did not terminate")
    end
    return nels, sum(ndofs) / ndofs[1]
end

# nonzero fixed Dirichlet data on every boundary face (the operator-on-zero gate)
function mg_nzbcs(mesh, T)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost([T(0.3 + 0.2im) * (1 + 0.1k) for k in 1:nf])
    end
    return bcs
end
mg_nzfactory(T, β) =
    mesh -> MeshOperator(T, mesh, (mg_A11, mg_A22), (x, y) -> β, (mg_zero2, mg_zero2);
                         C = MG_C, bcs = mg_nzbcs(mesh, T))

# a smooth (coarse-grid-mode) starting error: exp(2πi(x+y)) at the element nodes
function mg_lowfreq_start(mesh)
    u = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(ComplexF64, blk, el, mesh.orders[e])
        nx, ny = mesh.orders[e]
        rgn = mg_erange(mesh, e)
        for j in 1:ny, i in 1:nx
            u[rgn[EllipticCore.tensor_lin(i, j, nx)]] = cis(2pi * (em.x[1][i] + em.x[2][j]))
        end
    end
    return u
end

function mg_contraction_from(mg, u0; niter = 8)
    x = copy(u0)
    b = zeros(ComplexF64, length(u0))
    ns = Float64[norm(x)]
    for _ in 1:niter
        vcycle!(x, mg, b)
        push!(ns, norm(x))
    end
    return ns[end] / ns[end-1]
end

@testset "multigrid hierarchy (multigrid-hierarchy)" begin

    # ---------------------------------------------------------------- (1) coarsening ----
    @testset "h-coarsening down to the irreducible block structure" begin
        mesh = mg_uniform_mesh(Float64, 2, 4)
        co = coarsen_mesh(mesh)
        @test co !== nothing
        cmesh, parent = co
        @test length(mesh.elements) == 32
        @test length(cmesh.elements) == 8
        @test length(parent) == length(mesh.elements)
        # every fine element maps to an ancestor, at a per-dimension level gap of 0 or 1
        for (e, el) in enumerate(mesh.elements)
            c = cmesh.elements[parent[e]]
            @test c.block == el.block
            for d in 1:2
                Δ = el.segments[d].level - c.segments[d].level
                @test 0 <= Δ <= 1
                @test (el.segments[d].index >> Δ) == c.segments[d].index
            end
        end
        # ... and the coarse level tiles + stays 2:1 balanced (Mesh2D already validated it)
        @test two_to_one_balanced(cmesh.elements)
        c2 = coarsen_mesh(cmesh)
        @test c2 !== nothing
        @test length(c2[1].elements) == 2                    # ONE element per block
        @test coarsen_mesh(c2[1]) === nothing                # irreducible: fully coarsened

        # p-coarsening by the MIN rule (spectre Multigrid.hpp), on the 2:1 + p mesh
        mm = mg_mixed_mesh(Float64)
        cm, pm = coarsen_mesh(mm)
        @test length(cm.elements) == 2
        for k in eachindex(cm.elements)
            kids = [e for e in eachindex(mm.elements) if pm[e] == k]
            for d in 1:2
                @test cm.orders[k][d] == minimum(mm.orders[e][d] for e in kids)
            end
        end
        @test two_to_one_balanced(cm.elements)
        # the nonconforming two-block mesh coarsens all the way to the block structure
        nc = mg_nonconf_mesh(Float64)
        lev = Any[nc]
        while (c = coarsen_mesh(lev[end])) !== nothing
            push!(lev, c[1])
            @test two_to_one_balanced(c[1].elements)
        end
        # 6 levels, not the 4 of the shipped cap rule: the repaired step merges EVERY
        # complete sibling group (35 -> 10 in one step, where the cap rule managed 35 -> 32).
        @test length(lev) == 6
        @test [length(l.elements) for l in lev] == [35, 10, 6, 4, 3, 2]

        # the hierarchy wraps exactly this: level counts and per-level element counts
        for (k, want) in ((2, [32, 8, 2]), (3, [128, 32, 8, 2]))
            m = mg_uniform_mesh(Float64, k, 4)
            mk = mg_factory(ComplexF64, MG_BETA_DEF)
            hier = MultigridHierarchy(mk(m), mk)
            @test length(hier.levels) == length(want)
            @test [length(l.mesh.elements) for l in hier.levels] == want
            @test hier.levels[1].mesh === m                  # finest level is the caller's
            @test isempty(hier.levels[end].parent)           # coarsest has no parent level
        end
        # maxlevels caps the depth (spectre InitialCoarseLevels)
        m = mg_uniform_mesh(Float64, 3, 4)
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        @test length(MultigridHierarchy(mk(m), mk; maxlevels = 2).levels) == 2
    end

    # ------------------------------------------------------------- (2) transfer trio ----
    @testset "transfer trio: prolongation exact, R∘P = I, adjointness" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        for (name, mesh) in (("affine two-block k=2", mg_uniform_mesh(Float64, 2, 4)),
                             ("2:1 + p nonconforming", mg_nonconf_mesh(Float64)),
                             ("mixed 5-element", mg_mixed_mesh(Float64)),
                             ("curved 1/r block", mg_curved_mesh(Float64, 2)))
            hier = MultigridHierarchy(mk(mesh), mk)
            chk = mg_transfer_check(hier)
            @test chk.prolong_exact < MG_TRIO_TOL
            @test chk.restrict_identity < MG_TRIO_TOL
            @test chk.adjointness < MG_TRIO_TOL
        end

        # explicit R∘P = I on random coarse data through the public transfer calls
        mesh = mg_uniform_mesh(Float64, 2, 4)
        hier = MultigridHierarchy(mk(mesh), mk)
        rng = MGRNG(0x2024)
        for ℓ in 1:(length(hier.levels)-1)
            nf = hier.levels[ℓ].mesh.ndofs
            nc = hier.levels[ℓ+1].mesh.ndofs
            uc = mg_cvec!(rng, nc)
            uf = zeros(ComplexF64, nf)
            mg_prolong!(uf, hier, ℓ, uc)
            back = zeros(ComplexF64, nc)
            mg_restrict_field!(back, hier, ℓ, uf)
            @test maximum(abs.(back .- uc)) < MG_TRIO_TOL
        end
    end

    @testset "transfer trio at BigFloat(256)" begin
        setprecision(BigFloat, 256) do
            RT = BigFloat
            T = Complex{BigFloat}
            blk = Block(11, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1)))
            els = [ElementId(11, SegmentId(1, ix), SegmentId(1, iy)) for iy in 0:1 for ix in 0:1]
            mesh = Mesh2D(blk, els, (4, 4))
            mk = mg_factory(T, T(1) + T(1) / 2 * im)
            hier = MultigridHierarchy(mk(mesh), mk)
            chk = mg_transfer_check(hier)
            @test chk.prolong_exact < MG_TRIO_TOL_BIG
            @test chk.restrict_identity < MG_TRIO_TOL_BIG
            @test chk.adjointness < MG_TRIO_TOL_BIG
            @test chk.prolong_exact isa BigFloat            # precision inherited, not Float64
        end
    end

    # ------------------------------- (4) mass-weighted / Galerkin identities + scope ----
    @testset "mass-weighted adjointness and the massive Galerkin identity" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        function mass_gates(mesh)
            hier = MultigridHierarchy(mk(mesh), mk)
            rng = MGRNG(0x7777)
            adj = 0.0
            mid = 0.0
            for ℓ in 1:(length(hier.levels)-1)
                fm = hier.levels[ℓ].mesh
                cm = hier.levels[ℓ+1].mesh
                Mf = mg_massdiag(ComplexF64, fm)
                Mc = mg_massdiag(ComplexF64, cm)
                v = mg_cvec!(rng, fm.ndofs)
                u = mg_cvec!(rng, cm.ndofs)
                Rv = zeros(ComplexF64, cm.ndofs)
                mg_restrict_field!(Rv, hier, ℓ, v)
                Pu = zeros(ComplexF64, fm.ndofs)
                mg_prolong!(Pu, hier, ℓ, u)
                lhs = dot(u, Mc .* Rv)
                rhs = dot(Pu, Mf .* v)
                adj = max(adj, abs(lhs - rhs) / abs(rhs))
                # massive restriction == M_c ∘ R ∘ M_f⁻¹
                rm = zeros(ComplexF64, cm.ndofs)
                mg_restrict_residual!(rm, hier, ℓ, v)
                rl = zeros(ComplexF64, cm.ndofs)
                mg_restrict_field!(rl, hier, ℓ, v ./ Mf)
                mid = max(mid, maximum(abs.(rm .- Mc .* rl)) / maximum(abs.(rm)))
            end
            return adj, mid
        end
        for mesh in (mg_uniform_mesh(Float64, 2, 4), mg_nonconf_mesh(Float64))
            adj, mid = mass_gates(mesh)                      # AFFINE maps: both exact
            @test adj < MG_TRIO_TOL
            @test mid < MG_TRIO_TOL
        end
        # CURVED map: both identities FAIL at O(1) — the documented scope boundary. The
        # Galerkin restriction Pᵀ (what the V-cycle uses) is exact for any map; the
        # physical-mass forms are not, because the diagonal Gauss mass does not integrate
        # a non-constant Jacobian exactly.
        adj, mid = mass_gates(mg_curved_mesh(Float64, 2))
        @test adj > MG_CURVED_SCOPE
        @test mid > MG_CURVED_SCOPE
    end

    # --------------------------------------------------- (3) V-cycle contraction -------
    @testset "V-cycle contraction < 0.2, independent of the level count" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        ρs = Float64[]
        nlev = Int[]
        for k in (2, 3, 4)
            mesh = mg_uniform_mesh(Float64, k, 4)
            hier = MultigridHierarchy(mk(mesh), mk)
            mg = MultigridV(hier)
            push!(nlev, length(hier.levels))
            r = mg_contraction(mg, mesh.ndofs)
            push!(ρs, r[end])
            @test r[end] < MG_RHO_GATE
        end
        @test nlev == [3, 4, 5]
        @test maximum(ρs) - minimum(ρs) < MG_RHO_SPREAD      # LEVEL INDEPENDENCE
        # mildly indefinite complex Helmholtz on the deepest hierarchy
        mkh = mg_factory(ComplexF64, MG_BETA_MILD)
        mesh = mg_uniform_mesh(Float64, 3, 4)
        mgh = MultigridV(MultigridHierarchy(mkh(mesh), mkh))
        @test mg_contraction(mgh, mesh.ndofs)[end] < MG_RHO_GATE
    end

    @testset "V-cycle on nonconforming (2:1 h and p) and curved meshes" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        m1 = mg_nonconf_mesh(Float64)
        # 0.1983 measured (repair round 2): the repaired coarsening merges every complete
        # sibling group, so the coarse space here is genuinely coarser than the cap rule's
        # 35 -> 32 near-copy (which measured 0.1671). Still a contraction, gate MG_RHO_NC.
        @test mg_contraction(MultigridV(MultigridHierarchy(mk(m1), mk)), m1.ndofs)[end] <
              MG_RHO_NC
        m2 = mg_mixed_mesh(Float64)
        @test mg_contraction(MultigridV(MultigridHierarchy(mk(m2), mk)), m2.ndofs)[end] <
              MG_RHO_TINY
        m3 = mg_curved_mesh(Float64, 2)
        @test mg_contraction(MultigridV(MultigridHierarchy(mk(m3), mk)), m3.ndofs)[end] <
              MG_RHO_CURVED
    end

    # ------------------------------------- (5) independent reference + falsifiers ------
    @testset "independent reference V-cycle and the Galerkin falsifiers" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        mesh = mg_uniform_mesh(Float64, 2, 4)
        hier = MultigridHierarchy(mk(mesh), mk)
        mg = MultigridV(hier)
        # a second, independently written V-cycle over the public transfer API
        rng = MGRNG(0x555)
        b = mg_cvec!(rng, mesh.ndofs)
        z1 = zeros(ComplexF64, mesh.ndofs)
        vcycle!(z1, mg, b)
        z2 = zeros(ComplexF64, mesh.ndofs)
        mg_ref_vcycle!(z2, hier, b, 1, 2, lu(build_matrix(hier.levels[end].op)))
        @test norm(z1 - z2) / norm(z1) < 1e-13
        # GALERKIN falsifier: restricting the massive residual with the NODAL L2
        # restriction instead of Pᵀ badly degrades the contraction
        @test mg_ref_contraction(hier, 2)[end] < MG_RHO_GATE
        @test mg_ref_contraction(hier, 2; massive = false)[end] > MG_RHO_NONMASSIVE
        # SIGN falsifier: flip one element's prolongation in place
        hb = MultigridHierarchy(mk(mesh), mk)
        hb.levels[1].P[1][1] .*= -1
        chk = mg_transfer_check(hb)
        @test chk.prolong_exact > MG_CURVED_SCOPE
        @test chk.restrict_identity > MG_CURVED_SCOPE
        @test mg_contraction(MultigridV(hb), mesh.ndofs)[end] > MG_RHO_SIGN
        # JACOBIAN falsifier: drop the half-interval Jacobian of the restriction. It
        # breaks R∘P = I and the adjointness but NOT the V-cycle — R is not on the
        # V-cycle path (the V-cycle restricts the massive residual with Pᵀ).
        hj = MultigridHierarchy(mk(mesh), mk)
        for e in eachindex(hj.levels[1].R)
            hj.levels[1].R[e][1] .*= 2
        end
        chkj = mg_transfer_check(hj)
        @test chkj.prolong_exact < MG_TRIO_TOL
        @test chkj.restrict_identity > MG_CURVED_SCOPE
        @test chkj.adjointness > MG_CURVED_SCOPE
        @test mg_contraction(MultigridV(hj), mesh.ndofs)[end] < MG_RHO_GATE
    end

    # ------------------------------------------------- (6) GMRES iteration contrast ----
    @testset "GMRES ∘ MultigridV beats GMRES ∘ Schwarz on b := A·u_ex" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        counts = Int[]
        for k in (2, 3)
            mesh = mg_uniform_mesh(Float64, k, 4)
            op = mk(mesh)
            hier = MultigridHierarchy(op, mk)
            mg = MultigridV(hier)
            rng = MGRNG(0xfeed)
            uex = mg_cvec!(rng, mesh.ndofs)
            b = similar(uex)
            apply!(b, op, uex)
            rmg = gmres(op, b; abstol = 0.0, reltol = 1e-10, maxiter = 400, Minv! = mg)
            sw = schwarz_preconditioner(hier.levels[1].smoother; passes = 2)
            rsw = gmres(op, b; abstol = 0.0, reltol = 1e-10, maxiter = 400, Minv! = sw)
            @test rmg.converged
            @test rsw.converged
            @test rmg.residual <= 1e-10 * norm(b)
            @test rmg.iterations <= MG_GMRES_GATE
            @test rmg.iterations < rsw.iterations            # the headline contrast
            push!(counts, rmg.iterations)
        end
        @test abs(counts[1] - counts[2]) <= MG_GMRES_SPREAD  # level-independent count
    end

    # --------------------------------------------------------- (7) bottom-solver seam --
    @testset "bottom-solver seam" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        mesh = mg_uniform_mesh(Float64, 2, 4)
        hier = MultigridHierarchy(mk(mesh), mk)
        cop = hier.levels[end].op
        bot = dense_bottom_solver(cop)
        rng = MGRNG(0x314)
        bc = mg_cvec!(rng, hier.levels[end].mesh.ndofs)
        zc = similar(bc)
        bot(zc, bc)
        yc = similar(zc)
        apply!(yc, cop, zc)
        @test norm(yc - bc) / norm(bc) < MG_BOTTOM_TOL
        # a caller-supplied callable is honored — the seam the multigrid-bottom-solver
        # node (cached sparse UMFPACK LU) plugs into
        calls = Ref(0)
        counting = (z, b) -> (calls[] += 1; bot(z, b))
        mg = MultigridV(hier; bottom = counting)
        u = zeros(ComplexF64, mesh.ndofs)
        vcycle!(u, mg, mg_cvec!(rng, mesh.ndofs))
        @test calls[] == 1                                   # once per V-cycle
        # bottom = nothing ⇒ the smoother runs at the bottom (UseBottomSolver: False twin).
        # FINDING (round 4): on THIS coarsest grid — one element per block, 2 elements —
        # two overlapping Schwarz passes already solve the bottom system essentially
        # exactly, so the two settings contract identically (0.1714 vs 0.1714). The
        # bottom solver earns its keep on the CRP's 5×3-block coarsest grid, not here.
        mgn = MultigridV(hier; bottom = nothing)
        r = mg_contraction(mgn, mesh.ndofs)
        @test r[end] < MG_RHO_GATE                           # still a contraction
        @test abs(r[end] - mg_contraction(MultigridV(hier), mesh.ndofs)[end]) < 0.01
    end

    # ------------------------------------------- (8) cache reuse + allocation stability -
    @testset "Schwarz LU cache reuse across V-cycles, stable allocations" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        mesh = mg_uniform_mesh(Float64, 3, 4)
        hier = MultigridHierarchy(mk(mesh), mk)
        nsub = sum(length(l.smoother.subdomains) for l in hier.levels)
        builds0 = sum(l.smoother.cache.builds for l in hier.levels)
        @test builds0 == nsub                                # exactly one LU per subdomain
        @test nsub == 170                                    # calibrated 4-level total
        mg = MultigridV(hier)
        rng = MGRNG(0x99)
        v = mg_cvec!(rng, mesh.ndofs)
        z = zeros(ComplexF64, mesh.ndofs)
        mg(z, v)                                             # warm-up
        a1 = @allocated mg(z, v)
        a2 = @allocated mg(z, v)
        a3 = @allocated mg(z, v)
        @test a1 == a2 == a3                                 # constant per apply
        @test sum(l.smoother.cache.builds for l in hier.levels) == builds0   # zero rebuilds
        # the Minv! adapter is the pure linear map v ↦ M⁻¹v (independent of z's contents)
        z2 = mg_cvec!(rng, mesh.ndofs)
        mg(z2, v)
        @test z2 == z
    end

    # ------------------------------------------------------------- (9) BigFloat V-cycle -
    @testset "BigFloat(256) V-cycle and GMRES" begin
        setprecision(BigFloat, 256) do
            RT = BigFloat
            T = Complex{BigFloat}
            blk = Block(11, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1)))
            els = [ElementId(11, SegmentId(1, ix), SegmentId(1, iy)) for iy in 0:1 for ix in 0:1]
            mesh = Mesh2D(blk, els, (4, 4))
            mk = mg_factory(T, T(1) + T(1) / 2 * im)
            op = mk(mesh)
            hier = MultigridHierarchy(op, mk)
            mg = MultigridV(hier)
            rng = MGRNG(0x1111)
            x = convert(Vector{T}, mg_cvec!(rng, mesh.ndofs))
            b = zeros(T, mesh.ndofs)
            ns = [Float64(norm(x))]
            for _ in 1:5
                vcycle!(x, mg, b)
                push!(ns, Float64(norm(x)))
            end
            @test ns[end] / ns[end-1] < MG_RHO_GATE
            uex = convert(Vector{T}, mg_cvec!(rng, mesh.ndofs))
            bb = similar(uex)
            apply!(bb, op, uex)
            r = gmres(op, bb; abstol = BigFloat(1e-30), reltol = 0, maxiter = 100, Minv! = mg)
            @test r.converged
            @test r.residual < BigFloat(1e-30)
            @test norm(r.x - uex) / norm(uex) < BigFloat(1e-28)
        end
    end


    # ------------------------------------- (H1) anisotropic coarsening: the repair ------
    @testset "anisotropic h-coarsening: counterexample, census, geometric shrink" begin
        # (a) the validator's 18-element counterexample — fatal to the shipped cap rule
        m = mg_counterexample()
        @test two_to_one_balanced(m.elements)
        nels, _ = mg_chain_report(m)
        @test nels == [18, 9, 5, 3, 2, 1]
        # (b) census over AMR-reachable ANISOTROPIC meshes (the CRP's Isotropy: Anisotropic
        #     policy class, anc/ScalarSelfForce.yaml:145): every level of every hierarchy
        #     must be a valid, balanced, nested, exactly-tiling Mesh2D. The shipped cap rule
        #     failed 244/300 of these (104 balance violations, 140 face-coverage errors).
        ok = 0
        for seed in 1:MG_CENSUS_N
            try
                mg_chain_report(mg_aniso_mesh(seed, 3))
                ok += 1
            catch
            end
        end
        @test ok == MG_CENSUS_N
        ok2 = 0
        for seed in 501:560
            try
                mg_chain_report(mg_aniso_mesh(seed, 3; nblocks = 2))
                ok2 += 1
            catch
            end
        end
        @test ok2 == 60                                   # cross-block seams included
        # (c) uniform meshes shrink geometrically (spectre Hierarchy.cpp:20-46 parity)
        blk = Block(31, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
        u8 = Mesh2D(blk, [ElementId(31, SegmentId(3, ix), SegmentId(3, iy))
                          for iy in 0:7 for ix in 0:7], (4, 4))
        nels8, ratio8 = mg_chain_report(u8)
        @test nels8 == [64, 16, 4, 1]
        @test ratio8 < 1.34                               # calibrated 1.328
        # (d) the graded puncture cascade and an anisotropic mesh both build cleanly
        @test mg_chain_report(mg_graded_mesh(3, 3))[1][1] == 73
        @test length(mg_chain_report(mg_aniso_mesh(7, 3))[1]) > 1
    end

    @testset "contraction on graded, anisotropic and curved-nonconforming meshes" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        for (mesh, gate) in ((mg_graded_mesh(2, 3), MG_RHO_GATE),      # calibrated 0.1712
                             (mg_graded_mesh(3, 3), MG_RHO_GATE),      # calibrated 0.1714
                             (mg_aniso_mesh(7, 3), MG_RHO_ANISO),      # calibrated 0.1679
                             (mg_curved_nonconf(), MG_RHO_CURVED_NC))  # calibrated 0.3707
            hier = MultigridHierarchy(mk(mesh), mk)
            @test mg_contraction(MultigridV(hier), mesh.ndofs)[end] < gate
        end
        # the curved + 2:1 + p hierarchy also satisfies the transfer trio
        hc = MultigridHierarchy(mk(mg_curved_nonconf()), mk)
        chk = mg_transfer_check(hc)
        @test chk.prolong_exact < MG_TRIO_TOL
        @test chk.restrict_identity < MG_TRIO_TOL
        @test chk.adjointness < MG_TRIO_TOL
    end

    # ------------------------------------------------- (M3) the level-shrink policy ----
    @testset "shrink policy: fewer levels at unchanged contraction on graded meshes" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        mesh = mg_graded_mesh(2, 3)
        hfull = MultigridHierarchy(mk(mesh), mk)                       # default shrink = 1
        hskip = MultigridHierarchy(mk(mesh), mk; shrink = 0.5)
        @test length(hskip.levels) < length(hfull.levels)              # 5 vs 10
        cost(h) = sum(l.mesh.ndofs for l in h.levels) / h.levels[1].mesh.ndofs
        @test cost(hfull) > MG_COST_RATIO                              # calibrated 2.840
        @test cost(hskip) <= MG_COST_RATIO                             # calibrated 1.720
        @test mg_contraction(MultigridV(hskip), mesh.ndofs)[end] < MG_RHO_GATE  # 0.1712
        # composed transfers stay exact across the skipped step
        chk = mg_transfer_check(hskip)
        @test chk.prolong_exact < MG_TRIO_TOL
        @test chk.restrict_identity < MG_TRIO_TOL
        @test chk.adjointness < MG_TRIO_TOL
        @test_throws ArgumentError MultigridHierarchy(mk(mesh), mk; shrink = 0.0)
    end

    # ---------------------------- (H2) nonzero fixed Dirichlet data / operator-on-zero --
    @testset "nonzero fixed boundary data: the operator-on-zero contract" begin
        mesh = mg_uniform_mesh(Float64, 2, 4)
        mk = mg_nzfactory(ComplexF64, MG_BETA_DEF)
        op = mk(mesh)
        hier = MultigridHierarchy(op, mk)
        mg = MultigridV(hier)
        # the gate is not vacuous: fixed data are present on the fine AND the coarse levels
        @test norm(operator_on_zero(op)) > 1.0
        @test norm(operator_on_zero(hier.levels[2].op)) > 1.0
        rng = MGRNG(0x5150)
        f = mg_cvec!(rng, mesh.ndofs)
        blin = f .- operator_on_zero(op)                  # the LINEARIZED right-hand side
        r = gmres(op, blin; abstol = 0.0, reltol = 1e-12, maxiter = 200, Minv! = mg)
        @test r.converged
        @test norm(r.x - direct_solve(op, blin)) / norm(direct_solve(op, blin)) < 1e-10
        # the V-cycle is a LINEAR map: the coarse levels solve for a CORRECTION with a zero
        # guess and a restricted-residual source, so no fixed data may leak in. A regression
        # that fed the unadjusted b (or used apply_full!) to any level fails here.
        z1 = zeros(ComplexF64, mesh.ndofs)
        z2 = zeros(ComplexF64, mesh.ndofs)
        mg(z1, f)
        mg(z2, 2 .* f)
        @test norm(z2 .- 2 .* z1) / norm(2 .* z1) < 1e-13   # calibrated 0.0 (bitwise)
    end

    # ------------------------------------------------------ (H3) cross-AMR cache reuse --
    @testset "shared subdomain caches across rebuilds and AMR steps" begin
        mesh = mg_uniform_mesh(Float64, 3, 4)
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        cache = SubdomainCache()
        h1 = MultigridHierarchy(mk(mesh), mk; caches = cache)
        nsub = sum(length(l.smoother.subdomains) for l in h1.levels)
        @test cache.builds == nsub == 170
        @test cache.hits == 0
        h2 = MultigridHierarchy(mk(mesh), mk; caches = cache)   # identical meshes
        @test cache.builds == nsub                              # ZERO new factorizations
        @test cache.hits == nsub                                # every subdomain reused
        # after refining the FINEST level only, the untouched coarser levels still hit
        b0, h0 = cache.builds, cache.hits
        mref = h_refine(mesh, el -> true)
        h3 = MultigridHierarchy(mk(mref), mk; caches = cache)
        @test length(h3.levels) == length(h1.levels) + 1
        @test cache.builds - b0 == 512                          # only the new finest level
        @test cache.hits - h0 == nsub                           # all 170 coarser subdomains
        # per-level caches remain the default
        hdef = MultigridHierarchy(mk(mesh), mk)
        @test all(l.smoother.cache !== cache for l in hdef.levels)
    end

    # ------------------------------------- (M4) rediscretization, not a Galerkin product -
    @testset "coarse operators are rediscretized, not P'AP" begin
        mesh = mg_uniform_mesh(Float64, 2, 4)
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        op = mk(mesh)
        hier = MultigridHierarchy(op, mk)
        Ac = build_matrix(hier.levels[2].op)
        @test Ac == build_matrix(mk(hier.levels[2].mesh))     # the factory's own operator
        nc = hier.levels[2].mesh.ndofs
        P = zeros(ComplexF64, mesh.ndofs, nc)
        ec = zeros(ComplexF64, nc)
        col = zeros(ComplexF64, mesh.ndofs)
        for j in 1:nc
            fill!(ec, 0)
            ec[j] = 1
            mg_prolong!(col, hier, 1, ec)
            P[:, j] = col
        end
        G = transpose(P) * build_matrix(op) * P
        @test opnorm(Ac - G, Inf) / opnorm(Ac, Inf) > 0.1     # calibrated 1.3497
    end

    # --------------------------------------------- (M5) contraction over several starts --
    @testset "contraction over independent and low-frequency starts" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        for mesh in (mg_uniform_mesh(Float64, 3, 4), mg_graded_mesh(3, 3))
            mg = MultigridV(MultigridHierarchy(mk(mesh), mk))
            worst = maximum(mg_contraction(mg, mesh.ndofs; seed = s)[end]
                            for s in (0xabc, 0x1234, 0xfeed, 0x77))
            worst = max(worst, mg_contraction_from(mg, mg_lowfreq_start(mesh)))
            @test worst < MG_RHO_GATE                          # calibrated 0.1720 / 0.1727
        end
    end

    # -------------------------------------- (10) thread parity + argument validation ----
    @testset "thread parity, argument validation, docstrings" begin
        mk = mg_factory(ComplexF64, MG_BETA_DEF)
        mesh = mg_uniform_mesh(Float64, 2, 4)
        hier = MultigridHierarchy(mk(mesh), mk)
        mg = MultigridV(hier)
        rng = MGRNG(0x2222)
        b = mg_cvec!(rng, mesh.ndofs)
        zs = zeros(ComplexF64, mesh.ndofs)
        zt = zeros(ComplexF64, mesh.ndofs)
        vcycle!(zs, mg, b; threaded = false)
        vcycle!(zt, mg, b; threaded = true)
        @test zs == zt                                       # BITWISE, any thread count
        # length validation
        @test_throws ArgumentError vcycle!(zeros(ComplexF64, 3), mg, b)
        @test_throws ArgumentError vcycle!(zs, mg, zeros(ComplexF64, 3))
        @test_throws ArgumentError MultigridHierarchy(mk(mesh), mk; maxlevels = 0)
        @test_throws ArgumentError MultigridV(hier; passes = 0)
        # every export carries a docstring (Docs.Binding: @doc does not evaluate a loop
        # variable, so the binding is looked up explicitly)
        for nm in (:MultigridHierarchy, :MultigridV, :coarsen_mesh, :mg_prolong!,
                   :mg_restrict_residual!, :mg_restrict_field!, :mg_transfer_check,
                   :dense_bottom_solver, :vcycle!)
            doc = string(Base.Docs.doc(Base.Docs.Binding(EllipticCore, nm)))
            @test !occursin("No documentation found", doc)
        end
    end
end

# ============================================================================================
# The EXPLICIT COARSE-MESH CHAIN form of MultigridHierarchy (DAG row 56, linsolv-depth-gap).
#
# Written BEFORE the implementation; the red run is captured in
# emri-disk results/elliptic-ssf/paper_ELLIPTIC-SSF/evidence/linsolv-depth-gap/tests_red.log.
#
# WHY THE FORM EXISTS. The automatic constructor derives every coarse level by h-coarsening the
# CURRENT mesh, which can only ever express an h-hierarchy at a fixed (or MIN-rule) order. The
# reference's own hierarchy is not of that shape: `SolveScalarSelfForce.hpp:173` sets
# `amr::keep_coarse_grids = true`, `Amr/Actions/AdjustDomain.hpp:157-212` gives EVERY element a
# child at `grid_index + 1` on every AMR step (split children where the flag is `Split`, one
# same-shape child carrying the p increment otherwise) and
# `Amr/Actions/SendDataToChildren.hpp:43-70` keeps the parent alive as a coarser grid, so the
# grids the V-cycle walks (`Multigrid/ElementActions.hpp:96-135` sets the very
# `amr::Tags::{ParentId, ChildIds, ParentMesh}` triple `Multigrid.hpp:160-186` traverses) ARE the
# AMR history — consecutive levels differing by one h split AND one p increment. That chain
# cannot be recovered from the finest mesh alone; it has to be supplied.
#
# GATES: (1) handed the automatic hierarchy's own meshes the explicit form reproduces it
# BITWISE, level for level, transfer for transfer, and its V-cycle is bitwise the same map;
# (2) a p-ONLY link (identical elements, lower order) — the link h-coarsening cannot express —
# builds and passes the sealed transfer trio; (3) an hp chain of `hp_refine` steps (the
# reference's own shape) builds, passes the trio and contracts; (4) every way a chain can fail
# to be a coarsening is REFUSED with an ArgumentError rather than silently mis-transferred.
#
# GATE CONSTANTS calibrated PRE-RED on the throwaway prototype
# `evidence/linsolv-depth-gap/calibrate_chain.jl` (log `calibrate_chain.log`), which composes
# the sealed pieces the implementation will compose without being the implementation.
const MGC_TRIO_TOL = 1e-12      # calibrated 6.7e-16 … 1.7e-15 over the p-only and hp chains
const MGC_RHO_HP = 0.25         # calibrated 0.1708, 3-level hp chain (cf. MG_RHO_NC = 0.25)
const MGC_RHO_PONLY = 0.05      # calibrated 0.0045, 2-level p-only chain
@testset "explicit coarse-mesh chain (linsolv-depth-gap)" begin
    mk = mg_factory(ComplexF64, MG_BETA_DEF)

    # ---- (1) the automatic hierarchy's own chain reproduces it bitwise ------------------
    @testset "the automatic chain, supplied explicitly, is bitwise the same hierarchy" begin
        mesh = mg_uniform_mesh(Float64, 3, 4)
        op = mk(mesh)
        auto = MultigridHierarchy(op, mk)
        @test length(auto.levels) >= 3
        chain = [l.mesh for l in auto.levels[2:end]]
        expl = MultigridHierarchy(op, mk, chain)
        @test length(expl.levels) == length(auto.levels)
        for ℓ in eachindex(auto.levels)
            @test expl.levels[ℓ].mesh.elements == auto.levels[ℓ].mesh.elements
            @test expl.levels[ℓ].mesh.orders == auto.levels[ℓ].mesh.orders
            @test expl.levels[ℓ].parent == auto.levels[ℓ].parent
            @test expl.levels[ℓ].P == auto.levels[ℓ].P
            @test expl.levels[ℓ].R == auto.levels[ℓ].R
        end
        rng = MGRNG(0x5601)
        b = mg_cvec!(rng, mesh.ndofs)
        za = zeros(ComplexF64, mesh.ndofs)
        ze = zeros(ComplexF64, mesh.ndofs)
        vcycle!(za, MultigridV(auto), b)
        vcycle!(ze, MultigridV(expl), b)
        @test za == ze                                        # BITWISE the same V-cycle
        # an EMPTY chain is the single-level hierarchy (the bottom solve is the whole problem)
        one = MultigridHierarchy(op, mk, typeof(mesh)[])
        @test length(one.levels) == 1
        @test isempty(one.levels[1].parent)
    end

    # ---- (2) a p-only link: the link h-coarsening cannot express ------------------------
    @testset "p-only link (same elements, lower order) transfers correctly" begin
        fine = mg_uniform_mesh(Float64, 2, 5)
        coarse = mg_uniform_mesh(Float64, 2, 4)
        @test fine.elements == coarse.elements                # same h, different p
        @test fine.orders != coarse.orders
        hier = MultigridHierarchy(mk(fine), mk, [coarse])
        @test length(hier.levels) == 2
        @test hier.levels[1].parent == collect(1:length(fine.elements))
        chk = mg_transfer_check(hier)
        @test chk.prolong_exact < MGC_TRIO_TOL
        @test chk.restrict_identity < MGC_TRIO_TOL
        @test chk.adjointness < MGC_TRIO_TOL
        @test mg_contraction(MultigridV(hier), fine.ndofs)[end] < MGC_RHO_PONLY
    end

    # ---- (3) an hp chain of hp_refine steps: the reference's own shape ------------------
    @testset "hp chain of hp_refine steps builds, transfers and contracts" begin
        m0 = mg_uniform_mesh(Float64, 2, 4)
        # one corner element split, every other element p-incremented — the CRP criterion
        flags1 = Dict(m0.elements[1] => (true, true))
        m1 = hp_refine(m0, flags1)
        flags2 = Dict(m1.elements[1] => (true, true))
        m2 = hp_refine(m1, flags2)
        @test m2.ndofs > m1.ndofs > m0.ndofs
        hier = MultigridHierarchy(mk(m2), mk, [m1, m0])
        @test length(hier.levels) == 3
        @test [l.mesh.ndofs for l in hier.levels] == [m2.ndofs, m1.ndofs, m0.ndofs]
        chk = mg_transfer_check(hier)
        @test chk.prolong_exact < MGC_TRIO_TOL
        @test chk.restrict_identity < MGC_TRIO_TOL
        @test chk.adjointness < MGC_TRIO_TOL
        mg = MultigridV(hier)
        rho = mg_contraction(mg, m2.ndofs)
        @test rho[end] < MGC_RHO_HP
    end

    # ---- (4) every non-coarsening chain is refused --------------------------------------
    @testset "a chain that is not a coarsening is refused" begin
        m0 = mg_uniform_mesh(Float64, 2, 4)
        m1 = mg_uniform_mesh(Float64, 3, 4)                   # FINER than m0
        m2 = mg_uniform_mesh(Float64, 1, 4)                   # two levels coarser than m1
        # a "coarse" mesh that is finer: no ancestor of a fine element is present
        @test_throws ArgumentError MultigridHierarchy(mk(m0), mk, [m1])
        # a level gap of 2 per dimension: the transfers admit 0 or 1
        @test_throws ArgumentError MultigridHierarchy(mk(m1), mk, [m2])
        # a chain mesh from a different block lattice
        other = mg_curved_mesh(Float64, 2)
        @test_throws ArgumentError MultigridHierarchy(mk(m0), mk, [other])
    end
end
