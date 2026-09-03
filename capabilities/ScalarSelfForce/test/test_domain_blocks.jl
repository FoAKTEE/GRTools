# test_domain_blocks.jl — ESSF node `domain-blocks-map` (schedule W7-D, GLOBAL_DAG row 39),
# written BEFORE the implementation (TDD). The node delivers `Domain/blocks.jl`: the PURE
# DATA mapping from the sealed yaml-twin `RunConfig` (Domain/config.jl) to an EllipticCore
# `Mesh2D` — the CRP computational domain of arXiv:2606.04998 §III.A.1 (main.tex:541-549):
# 5 radial × 2 angular irreducible blocks under equatorial symmetry (the (r, z²) chart of the
# production yaml) or 5 × 3 without it (the (r, z) chart of the spectre non-equatorial twin),
# `LinearMap1D` everywhere except the outer wavezone column, which carries the 1/r
# compactification eq:compaction (main.tex:546-549) — plus the per-block METADATA the
# downstream nodes bind to: the vtu slicing region of each block, the worldtube
# (regularized-region) block, and the classification of every domain-boundary face.
#
# ORACLE (the DAG row's gate column) — the READ-ONLY mirrors, transcribed by hand below:
#   * `ref-paper/arxiv-2606.04998/src/anc/ScalarSelfForce.yaml:39-127` — the production
#     `DomainCreator.AlignedLattice`: `BlockBounds` {inner, r0−1.5w, r0−w, r0+w, r0+3w, R} ×
#     {0, dz², 1}, `Distributions` [Linear ×4, Inverse] × [Linear ×2], `RefinedLevels` (the
#     10 single-block patches; source column L+1), and `NullSlicingBlocks`
#     ["Block0", "Block5"] (v) / ["Block4", "Block9"] (u) — which pins BOTH the flat block
#     numbering (radial index fastest: Block(i,j) ↦ i + 5j) and the region of each column.
#   * `ref-code/spectre/tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:29-120` — the
#     NON-equatorial twin (`ImposeEquatorialSymmetry: False`): angular `BlockBounds`
#     [−1, −0.5, 0.5, 1] for the same worldtube half-width dz² = 0.25 (i.e. the z rows are the
#     SIGNED SQUARE ROOTS of the z² rows — the same physical domain), `Distributions`
#     [Linear ×3] angular, `NullSlicingBlocks` Block(0,j)/Block(4,j) for every row j, and 15
#     `RefinedLevels` patches whose refinement is the SAME function of the radial column as in
#     the 5×2 production table. Both facts are transcribed as literals here.
#   * `main.tex:541-543` (the 5×3 / 5×2 split and the v/t/t/t/u radial region order),
#     `main.tex:546-549` (eq:compaction), `main.tex:552-554` (h refinement + 2:1 balance),
#     `main.tex:529` (fig:domain caption: horizon r₊ → R = 1e4 M, z² = 0 equator → z² = 1
#     poles, regularized region = the red rectangle around the puncture).
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here):
#   * runtime-config-fixture (`Domain/config.jl`) — `RunConfig`/`DomainConfig`/`OrbitConfig`/
#     `BlockRefinement` and every derived number (v→t = r0 − 1.5w, t→u = r0 + 3w, the block
#     bounds, the L+1 source pre-refinement, the m = 0 branch). The bounds are asserted BOTH
#     against the config's own fields and against the yaml's Jinja2 EXPRESSIONS (in IEEE
#     doubles 10.0 − 4.8·1.5 = 2.8000000000000007, not the decimal-pretty 2.8).
#   * mesh2d-blocks / matrix-free-operator (EllipticCore) — `Block`, `LinearMap1D`,
#     `InverseRadialMap1D`, `ElementId`/`SegmentId`, `two_to_one_balanced`, `Mesh2D`,
#     `BlockConnection`, `boundary_faces`, `MeshOperator`, `apply_full!`, `DirichletGhost`.
#   * vtu-slicing-jumps — `vtu_slicing(orbit, config)` + `vtu_domain` (the region classifier);
#     first-order-flux-form / equatorial-symmetry-z2 — `fluxform_callables(orbit, config,
#     region)` (the RunConfig-driven chart selector). This node NEVER re-derives H, r_v, r_u
#     or any coefficient: it only decides WHICH sealed bundle each block gets.
#
# THIS node's gates (calibrated pre-red in the evidence dir, calibrate_blocks{,2,3}.jl):
#  (1) bounds/maps — radial bounds exactly the config's and the yaml expressions'; z rows
#      exact in both charts; `InverseRadialMap1D` ONLY in the u column with endpoint-exact
#      `map_forward(±1)`; `LinearMap1D` everywhere else.
#  (2) layout — 5×2 / 5×3 block and connection counts, the flat yaml block numbering against
#      the `NullSlicingBlocks` NAMES, per-block refinement from `RefinedLevels`, element
#      counts for L ∈ {0,1,2}, `two_to_one_balanced`, every interior block face connected
#      exactly once (reciprocal interface topology), and the boundary-face classification
#      (horizon / outer / equator / pole) against an independent enumeration.
#  (3) operator wiring — a `MeshOperator` built on a block's ROOT element from the sealed
#      per-region `fluxform_callables` with exact-trace `DirichletGhost`s (every boundary term
#      of `strong_residual` then vanishes identically, so `apply_full!` == M·(first-order
#      form)) reproduces an INDEPENDENT transcription of △_m Ψ (eq:triangle main.tex:247-251 +
#      the replacement rules main.tex:306-308, region-constant H from the printed vtu table
#      main.tex:316-322, chain-ruled to the z² chart) to 1e-10 — calibrated worst 6.2e-13 over
#      the whole parameter sweep, while the SAME check with the callables of a WRONG region
#      lands at 2.1e-2 … 7.5e+1 (the falsifier: this is what proves the region tags).
#      Plus (3d), the PRODUCTION-order statement, machine-checked rather than asserted in
#      prose: at P = 4 (element space = degree ≤ 3 per dimension) a probe that IS in the space
#      and whose angular flux product stays in it (deg_w ≤ N−2, since A^{z²} = 4w(1−w)/(r²+a²)
#      is degree 2 in w) reproduces △_m to ~1e-14; what survives at representable degree is
#      ANGULAR, not radial — deg (0, N−1) leaves 5.4e-2 while deg (N−1, 0) stays at 1.8e-14,
#      the radial sector cancelling identically because both sides write it as the collocation
#      of the flux; and the operator is CONSISTENT at P = 4: at fixed N = 4 the deviation
#      h-converges at rate 2.29-3.45 on the source column.
#  (3e) the :horizon precondition — `block_domain` REQUIRES the inner radial bound to be the
#      outer horizon r₊ (main.tex:529), because that placement is the only reason the inner
#      faces need no boundary data: A^rr = Δ/(r²+a²) is exactly zero AT r₊ but measures
#      5.88e-2 at r_in = 2 with a = 0.5 and −4.0e-2 inside the horizon (a = 0.9, r_in = 1.3).
#      Every spin variation below therefore moves the inner bound to that spin's r₊, and the
#      throw is gated.
#  (4) refinement/variation/precision — the source (worldtube) column carries L+1 in both
#      dimensions, element counts per block match 2^(nx+ny) for L ∈ {0,1,2}; parameter
#      variation (r0, w, dz², spin incl. a < 0, m = 0, m = 20) keeps every invariant;
#      BigFloat(256) bounds and maps inherit precision; the m = 0 config builds an IDENTICAL
#      lattice (its only difference is the absent outer BC — node 41's business).
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore
using ScalarSelfForce.EllipticCore: refelement, element_matrices, element_coefficients

# ---- independent transcriptions of the two mirrors (parse-free literals) -----------------
# anc/ScalarSelfForce.yaml:55-63 rendered at the header production example (yaml:9-12), in the
# yaml's OWN Jinja2 expressions so the IEEE doubles match what spectre parses
const DB_R = [2.0, 10.0 - 4.8 * 1.5, 10.0 - 4.8, 10.0 + 4.8, 10.0 + 4.8 * 3, 1.0e4]
const DB_Z2 = [0.0, 0.25, 1.0]                       # yaml:63
const DB_Z = [-1.0, -0.5, 0.5, 1.0]                  # spectre ScalarCircularOrbit.yaml:59
# yaml:78-110 / spectre :71-120 — the refinement as a function of the RADIAL COLUMN, L = 0
const DB_COLREF = [(0, 0), (0, 0), (1, 1), (0, 1), (1, 0)]
const DB_NR = 5
const DB_VCOL, DB_UCOL = 0, 4                        # main.tex:541-543 + NullSlicingBlocks

# per-column element counts and the flat lattice id (radial index fastest, yaml Block0..Block9)
db_colref(L) = [(L + a, L + b) for (a, b) in DB_COLREF]
db_bid(i, j) = i + DB_NR * j
db_nel_col(L, i) = 1 << sum(db_colref(L)[i+1])

# a RunConfig twin with ImposeEquatorialSymmetry = false. `crp_paper_config` hard-codes the
# yaml's `True` (yaml:36) and exposes no keyword for it, so the 5×3 layout — described in
# main.tex:541-543 and realized by the spectre non-equatorial input file — is exercised
# through a field-wise reconstruction of the sealed structs.
function db_full_z(cfg::RunConfig{T}) where {T}
    OC = typeof(cfg.orbit)
    oc = OC((f === :impose_equatorial_symmetry ? false : getfield(cfg.orbit, f)
             for f in fieldnames(OC))...)
    return RunConfig{T}((f === :orbit ? oc : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

db_orbit(cfg) = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)

# `block_domain` requires the inner radial bound to BE the outer horizon (main.tex:529), and
# the yaml's inner_radius = 2 is r₊ only at the production spin a = 0. Any spin variation must
# therefore move the inner bound with the horizon; `db_cfg` is `crp_paper_config` with that
# one substitution, so a = 0 reproduces the production config exactly.
db_cfg(; spin = 0.0, orbital_radius = 10.0, kwargs...) =
    crp_paper_config(; spin = spin, orbital_radius = orbital_radius,
                     inner_radius = CircularOrbit(spin, orbital_radius).rp, kwargs...)

# ---- independent expected topology, enumerated from the transcribed tables ----------------
# elements of the lattice, in the layout order the node contracts: block-major (radial index
# fastest, then angular row), ix-fastest within a block
function db_expected_elements(L, nz)
    els = ElementId[]
    for j in 0:nz-1, i in 0:DB_NR-1
        nx, ny = db_colref(L)[i+1]
        for iy in 0:(1<<ny)-1, ix in 0:(1<<nx)-1
            push!(els, ElementId(db_bid(i, j), SegmentId(nx, ix), SegmentId(ny, iy)))
        end
    end
    return els
end

# every domain-boundary face with its class, from the lattice perimeter
function db_expected_boundary(L, nz, chart)
    out = Dict{Tuple{ElementId,Int,Int},Symbol}()
    for j in 0:nz-1, i in 0:DB_NR-1
        nx, ny = db_colref(L)[i+1]
        for iy in 0:(1<<ny)-1, ix in 0:(1<<nx)-1
            el = ElementId(db_bid(i, j), SegmentId(nx, ix), SegmentId(ny, iy))
            i == 0 && ix == 0 && (out[(el, 1, 1)] = :horizon)
            i == DB_NR - 1 && ix == (1 << nx) - 1 && (out[(el, 1, 2)] = :outer)
            j == 0 && iy == 0 &&
                (out[(el, 2, 1)] = chart === :z2 ? :equator : :pole)
            j == nz - 1 && iy == (1 << ny) - 1 && (out[(el, 2, 2)] = :pole)
        end
    end
    return out
end

# ---- the spectral instrument for gate (3) ------------------------------------------------
# deterministic complex bivariate polynomial in the BLOCK-logical coordinates (no Random dep)
mutable struct DbLcg
    s::UInt64
end
_db_u!(g::DbLcg) = ((g.s = 0x5851f42d4c957f2d * g.s + 0x14057b7ef767814f;
                     (g.s >> 11) * 2.0^-52 - 1.0))
_db_c!(g::DbLcg) = complex(_db_u!(g), _db_u!(g))
const DB_DX, DB_DY = 6, 4
const DB_COEF = let g = DbLcg(0x2606049980000039)
    [_db_c!(g) for _ in 0:DB_DX, _ in 0:DB_DY]
end
db_poly(z1, z2) = sum(DB_COEF[i+1, j+1] * z1^i * z2^j for i in 0:DB_DX, j in 0:DB_DY)

# INDEPENDENT reference: eq:triangle (main.tex:247-251) + the hyperboloidal replacement rules
# (main.tex:306-308) with a region-CONSTANT boost H (H′ ≡ 0, main.tex:347) and the radial
# chart binding ∂_{r_*} = Δ/(r²+a²)∂_r (main.tex:356); raw Δ = r² − 2Mr + a². In the z² chart
# the printed angular operator is chain-ruled with ∂_z = 2z∂_w (w = z²):
#   ∂_z((1−z²)∂_zΨ) − 2m z ∂_zΨ = 4w(1−w)Ψ_ww + (2 − 6w − 4mw)Ψ_w .
function db_triangle(orb, m, H, chart, r, y, Dr, Dy, U)
    a, M, Om = orb.a, orb.M, orb.Omega
    yt = transpose(y)
    Dl = @. r^2 - 2M * r + a^2
    wr = @. r^2 + a^2
    lam = Dl ./ wr
    dUr = Dr * U
    dUy = U * transpose(Dy)
    D1 = lam .* dUr
    D2 = lam .* (Dr * (lam .* dUr))
    D1h = D1 .+ (im * m * Om * H) .* U
    D2h = D2 .+ (2im * m * Om * H) .* D1 .- (m^2 * Om^2 * H^2) .* U
    if chart === :z
        z2 = yt .^ 2
        Zang = ((1 .- z2) .* dUy) * transpose(Dy) .- 2m .* (yt .* dUy) .- m * (m + 1) .* U
        Sig = @. r^2 + a^2 * z2
    else
        Zang = (4 .* yt .* (1 .- yt)) .* ((U * transpose(Dy)) * transpose(Dy)) .+
               (2 .- 6 .* yt .- 4m .* yt) .* dUy .- m * (m + 1) .* U
        Sig = @. r^2 + a^2 * yt
    end
    @. -(wr / Dl) * D2h + (2a^2 / (r * wr) - 2im * m * a / Dl) * D1h -
       (Zang + m^2 * Om * (2M * r * (Om * wr - 2a) + Om * Dl * Sig) / Dl * U +
        2 * (a^2 - im * m * a * r - M * r) / r^2 * U) / wr
end

db_H_table(region) = region === :v ? -1 : region === :t ? 0 : 1   # main.tex:316-322

# Build a one-element mesh on the block's ROOT element (a complete tiling of the block, so
# every face is a domain-boundary face), apply the assembled operator built from the sealed
# `region` bundle to an exactly-representable polynomial with exact-trace Dirichlet ghosts,
# and compare against M·db_triangle with boost `Href` (defaults to the region's table value).
function db_identity(cfg, region, bid, N; Href = nothing, k = 0, deg = (DB_DX, DB_DY))
    orb = db_orbit(cfg)
    dom = block_domain(orb, cfg)
    blk = dom.mesh.blocks[findfirst(b -> b.id == bid, dom.mesh.blocks)]
    els = ElementId[]
    for iy in 0:(1<<k)-1, ix in 0:(1<<k)-1          # ix fastest, the layout convention
        push!(els, ElementId(bid, SegmentId(k, ix), SegmentId(k, iy)))
    end
    mesh1 = Mesh2D(blk, els, (N, N))
    fc = fluxform_callables(orb, cfg, region)
    dx, dy = deg
    P(z1, z2) = sum(DB_COEF[i+1, j+1] * z1^i * z2^j for i in 0:dx, j in 0:dy)
    xi = refelement(Float64, N).nodes
    # block-logical coordinates of one element's nodes along `dim`
    function zeta(el, dim)
        lo, hi = logical_interval(Float64, el, dim)
        return [(lo * (1 - x) + hi * (1 + x)) / 2 for x in xi]
    end
    u = zeros(ComplexF64, mesh1.ndofs)
    for (e, el) in enumerate(els)
        z1, z2 = zeta(el, 1), zeta(el, 2)
        off = mesh1.offsets[e]
        u[off+1:off+N*N] .= vec([P(p, q) for p in z1, q in z2])
    end
    bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}()
    for (el, d, sidx) in boundary_faces(mesh1)
        lo, hi = logical_interval(Float64, el, d)
        zf = sidx == 1 ? lo : hi
        zt = zeta(el, 3 - d)
        tr = d == 1 ? [P(zf, q) for q in zt] : [P(p, zf) for p in zt]
        bcs[(el, d, sidx)] = DirichletGhost(ComplexF64.(tr))
    end
    op = MeshOperator(ComplexF64, mesh1, fc.A, fc.beta, fc.gamma;
                      C = cfg.dg.penalty_parameter, bcs = bcs)
    y = apply_full!(zeros(ComplexF64, mesh1.ndofs), op, u)
    ref = zeros(ComplexF64, mesh1.ndofs)
    D = refelement(Float64, N).D
    for (e, el) in enumerate(els)
        em = element_matrices(ComplexF64, blk, el, (N, N))
        z1, z2 = zeta(el, 1), zeta(el, 2)
        U = [P(p, q) for p in z1, q in z2]
        R = db_triangle(orb, cfg.orbit.m_mode_number,
                        Href === nothing ? db_H_table(region) : Href, dom.chart,
                        em.x[1], em.x[2], (1 ./ em.jac[1]) .* D, (1 ./ em.jac[2]) .* D, U)
        off = mesh1.offsets[e]
        ref[off+1:off+N*N] .= em.M * vec(R)
    end
    return maximum(abs, y .- ref) / maximum(abs, ref)
end

@testset "Domain block lattice → EllipticCore Mesh2D (domain-blocks-map)" begin

    @testset "exports" begin
        for s in (:BlockDomain, :block_domain, :block_lattice_bounds, :block_region,
                  :block_lattice_index, :worldtube_blocks, :boundary_face_classes,
                  :boundary_class, :block_coefficients, :domain_coefficients)
            @test Base.isexported(ScalarSelfForce, s)
        end
    end

    @testset "(1) block bounds, z rows, and the per-column coordinate maps" begin
        cfg = crp_paper_config()
        orb = db_orbit(cfg)
        rb, zb, chart = block_lattice_bounds(cfg)
        # equal to the SEALED config's own field and to the yaml Jinja2 expressions
        @test rb == cfg.domain.block_bounds_r
        @test rb == DB_R
        @test chart === :z2
        @test zb == cfg.domain.block_bounds_z2
        @test zb == DB_Z2
        # the decimal-pretty reading is only correct to roundoff — the doubles are the yaml's
        @test rb ≈ [2.0, 2.8, 5.2, 14.8, 24.4, 1.0e4] atol = 1e-14
        @test rb[2] != 2.8                       # 10.0 − 4.8·1.5 = 2.8000000000000007
        @test issorted(rb) && allunique(rb)
        # the yaml's inner_radius IS the outer horizon at the production spin (main.tex:529):
        # 2 = r₊(a = 0) exactly — the precondition `block_domain` enforces
        @test rb[1] == orb.rp == 2.0
        # the vtu transitions ARE block bounds (main.tex:323) — the invariant that makes the
        # region a per-block property
        @test cfg.orbit.hyperboloidal_slicing_transitions[1] == rb[2]
        @test cfg.orbit.hyperboloidal_slicing_transitions[3] == rb[5]

        # the non-equatorial twin: z rows are the signed square roots of the z² rows
        cfz = db_full_z(cfg)
        rb2, zb2, chart2 = block_lattice_bounds(cfz)
        @test rb2 == rb
        @test chart2 === :z
        @test zb2 == DB_Z
        @test zb2[3]^2 == cfg.domain.block_bounds_z2[2]      # 0.5² == 0.25 exactly
        @test zb2 == -reverse(zb2)

        dom = block_domain(orb, cfg)
        @test dom isa BlockDomain{Float64}
        @test (dom.nr, dom.nz) == (5, 2)
        @test dom.bounds_r == rb && dom.bounds_z == zb && dom.chart === :z2
        # Distributions[r] = [Linear ×4, Inverse] (yaml:66); Distributions[z²] = Linear (yaml:68)
        for j in 0:dom.nz-1, i in 0:DB_NR-1
            blk = dom.mesh.blocks[findfirst(b -> b.id == db_bid(i, j), dom.mesh.blocks)]
            mx, my = blockmap(blk, 1), blockmap(blk, 2)
            @test my isa LinearMap1D
            @test map_forward(my, -1.0) == zb[j+1] && map_forward(my, 1.0) == zb[j+2]
            if i == DB_NR - 1
                @test mx isa InverseRadialMap1D          # eq:compaction, main.tex:546-549
                @test mx.r_u == rb[i+1] && mx.R == rb[i+2]
                @test map_forward(mx, -1.0) == rb[i+1]   # endpoint-EXACT by contract
                @test map_forward(mx, 1.0) == rb[i+2]
                # 1/r is affine in ξ: the midpoint is the harmonic mean of the endpoints
                @test map_forward(mx, 0.0) ≈ 2 * rb[i+1] * rb[i+2] / (rb[i+1] + rb[i+2])
            else
                @test mx isa LinearMap1D
                @test map_forward(mx, -1.0) == rb[i+1] && map_forward(mx, 1.0) == rb[i+2]
            end
        end
    end

    @testset "(2) 5×2 / 5×3 layout, block numbering, connections, refinement, elements" begin
        for (cfg, nz, chart) in ((crp_paper_config(), 2, :z2),
                                 (db_full_z(crp_paper_config()), 3, :z))
            orb = db_orbit(cfg)
            dom = block_domain(orb, cfg)
            @test (dom.nr, dom.nz) == (DB_NR, nz)
            @test dom.chart === chart
            @test length(dom.mesh.blocks) == DB_NR * nz
            @test sort([b.id for b in dom.mesh.blocks]) == collect(0:DB_NR*nz-1)
            # flat yaml numbering: radial index fastest
            for j in 0:nz-1, i in 0:DB_NR-1
                @test dom.block_ids[i+1, j+1] == db_bid(i, j)
                @test block_lattice_index(dom, db_bid(i, j)) == (i, j)
            end
            # one connection per interior block face, each exactly once
            @test length(dom.mesh.connections) == (DB_NR - 1) * nz + DB_NR * (nz - 1)
            @test allunique((c.block, c.dim) for c in dom.mesh.connections)
            for c in dom.mesh.connections
                i, j = block_lattice_index(dom, c.block)
                ni, nj = block_lattice_index(dom, c.neighbor)
                @test (c.dim == 1 ? (ni, nj) == (i + 1, j) : (ni, nj) == (i, j + 1))
            end
            # RefinedLevels (yaml:78-110) → per-block initial refinement, elements, 2:1
            L = cfg.domain.initial_levels[1]
            for j in 0:nz-1, i in 0:DB_NR-1
                @test dom.refinements[db_bid(i, j)] == db_colref(L)[i+1]
                @test count(e -> e.block == db_bid(i, j), dom.mesh.elements) ==
                      db_nel_col(L, i)
            end
            @test dom.mesh.elements == db_expected_elements(L, nz)
            @test two_to_one_balanced(dom.mesh.elements)
            @test dom.mesh.ndofs ==
                  length(dom.mesh.elements) * prod(cfg.domain.initial_grid_points)
            @test all(o -> o == cfg.domain.initial_grid_points, dom.mesh.orders)

            # every interior element face carries an interface, reciprocally
            bset = Set(boundary_faces(dom.mesh))
            nint = 0
            for (e, el) in enumerate(dom.mesh.elements), d in 1:2, sidx in 1:2
                iface = dom.mesh.faces[e][d][sidx]
                if (el, d, sidx) in bset
                    @test iface === nothing
                else
                    @test iface !== nothing
                    nint += 1
                    for k in eachindex(iface.nbs)
                        nb = iface.nbs[k]
                        back = dom.mesh.faces[nb][d][3-sidx]
                        @test back !== nothing
                        @test back.nbs[iface.nbslot[k]] == e
                    end
                end
            end
            @test nint == 4 * length(dom.mesh.elements) - length(bset)
        end
        # counts pinned as literals at the production point (L = 0)
        let d2 = block_domain(db_orbit(crp_paper_config()), crp_paper_config()),
            d3 = block_domain(db_orbit(crp_paper_config()), db_full_z(crp_paper_config()))
            @test (length(d2.mesh.blocks), length(d2.mesh.elements),
                   length(d2.mesh.connections), d2.mesh.ndofs) == (10, 20, 13, 320)
            @test (length(d3.mesh.blocks), length(d3.mesh.elements),
                   length(d3.mesh.connections), d3.mesh.ndofs) == (15, 30, 22, 480)
        end
    end

    @testset "(2b) boundary-face classification: horizon / outer / equator / pole" begin
        for (cfg, nz, chart, tot) in ((crp_paper_config(), 2, :z2, 18),
                                      (db_full_z(crp_paper_config()), 3, :z, 20))
            dom = block_domain(db_orbit(cfg), cfg)
            cls = boundary_face_classes(dom)
            exp = db_expected_boundary(cfg.domain.initial_levels[1], nz, chart)
            @test cls == exp
            @test Set(keys(cls)) == Set(boundary_faces(dom.mesh))
            @test length(cls) == tot
            for f in keys(cls)
                @test boundary_class(dom, f) === cls[f]
            end
            counts = Dict(c => count(==(c), values(cls)) for c in Set(values(cls)))
            if chart === :z2
                @test counts == Dict(:horizon => 2, :outer => 2, :equator => 7, :pole => 7)
            else
                @test counts == Dict(:horizon => 3, :outer => 3, :pole => 14)
                @test !(:equator in keys(counts))     # the equator is INTERIOR without symmetry
            end
            # the physical coordinate of every classified face
            rb, zb, _ = block_lattice_bounds(cfg)
            for ((el, d, sidx), c) in cls
                blk = dom.mesh.blocks[findfirst(b -> b.id == el.block, dom.mesh.blocks)]
                lo, hi = logical_interval(Float64, el, d)
                x = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
                if c === :horizon
                    @test x == rb[1]
                elseif c === :outer
                    @test x == rb[end]
                elseif c === :equator
                    @test x == zb[1] == 0
                else
                    @test c === :pole
                    @test (x == zb[1] || x == zb[end])
                end
            end
        end
    end

    @testset "(3a) region tags, worldtube block, and the sealed-bundle wiring" begin
        for cfg in (crp_paper_config(), db_cfg(spin = 0.5),
                    db_full_z(db_cfg(spin = -0.7)))
            orb = db_orbit(cfg)
            dom = block_domain(orb, cfg)
            v = vtu_slicing(orb, cfg)
            for j in 0:dom.nz-1, i in 0:DB_NR-1
                bid = db_bid(i, j)
                # the tag equals the SEALED classifier at the block's radial midpoint
                mid = (dom.bounds_r[i+1] + dom.bounds_r[i+2]) / 2
                @test block_region(dom, bid) === vtu_domain(v, mid)
                @test block_region(dom, bid) ===
                      (i == DB_VCOL ? :v : i == DB_UCOL ? :u : :t)
                # and at EVERY volume collocation node of every element of the block
                blk = dom.mesh.blocks[findfirst(b -> b.id == bid, dom.mesh.blocks)]
                for el in filter(e -> e.block == bid, dom.mesh.elements)
                    em = element_matrices(ComplexF64, blk, el, (4, 4))
                    @test all(r -> vtu_domain(v, r) === block_region(dom, bid), em.x[1])
                end
            end
            # the worldtube (regularized-region) block holds the particle at (r0, equator)
            wt = worldtube_blocks(dom)
            @test length(wt) == 1
            i, j = block_lattice_index(dom, only(wt))
            @test dom.bounds_r[i+1] < orb.r0 < dom.bounds_r[i+2]
            @test dom.bounds_z[j+1] <= 0 <= dom.bounds_z[j+2]
            @test block_region(dom, only(wt)) === :t
            @test only(wt) == (dom.chart === :z2 ? 2 : 7)
            @test dom.refinements[only(wt)] ==
                  (cfg.domain.initial_levels[1] + 1, cfg.domain.initial_levels[2] + 1)
            # per-block coefficients ARE the sealed per-region bundle (no new algebra)
            for bid in 0:DB_NR*dom.nz-1
                bc = block_coefficients(dom, orb, cfg, bid)
                fc = fluxform_callables(orb, cfg, block_region(dom, bid))
                r, y = (dom.bounds_r[1] + dom.bounds_r[2]) / 2, 0.17
                bid % DB_NR == 0 || (r = (dom.bounds_r[3] + dom.bounds_r[4]) / 2)
                @test bc.beta(r, y) == fc.beta(r, y)
                @test bc.gamma[1](r, y) == fc.gamma[1](r, y)
                @test bc.gamma[2](r, y) == fc.gamma[2](r, y)
                @test bc.A[1](r, y) == fc.A[1](r, y) && bc.A[2](r, y) == fc.A[2](r, y)
            end
        end
    end

    @testset "(3b) the global bundle == the per-block bundles, BITWISE, on every element" begin
        # the driver (solver-stack-ssf) needs ONE coefficient bundle for the whole mesh:
        # `MeshOperator` evaluates β and γ at volume nodes only and 𝔄 at faces too — and 𝔄 is
        # H-INDEPENDENT — so the radius-dispatching bundle is exact everywhere it is sampled
        for cfg in (db_cfg(spin = 0.5), db_full_z(crp_paper_config()))
            orb = db_orbit(cfg)
            dom = block_domain(orb, cfg)
            gc = domain_coefficients(dom, orb, cfg)
            for (e, el) in enumerate(dom.mesh.elements)
                blk = dom.mesh.blocks[findfirst(b -> b.id == el.block, dom.mesh.blocks)]
                bc = block_coefficients(dom, orb, cfg, el.block)
                em = element_matrices(ComplexF64, blk, el, dom.mesh.orders[e])
                cb = element_coefficients(em, blk, el, bc.A, bc.beta, bc.gamma)
                cg = element_coefficients(em, blk, el, gc.A, gc.beta, gc.gamma)
                @test cb.beta == cg.beta
                @test cb.gamma == cg.gamma
                @test cb.A == cg.A
                @test cb.Aface == cg.Aface
            end
        end
    end

    @testset "(3c) assembled operator == △_m per region (1e-10) + wrong-region falsifier" begin
        cfg = crp_paper_config()
        worst = 0.0
        for j in 0:1, i in 0:DB_NR-1
            bid = db_bid(i, j)
            dev = db_identity(cfg, i == 0 ? :v : i == DB_NR - 1 ? :u : :t, bid, 12)
            worst = max(worst, dev)
            @test dev <= 1e-10
        end
        println("domain-blocks-map gate (3c) worst region-identity rel dev: ", worst)
        # the falsifier: the reference keeps the printed-table boost of the block's TRUE
        # region while the operator is built from another region's callables
        for (bid, good) in ((0, :v), (4, :u), (2, :t))
            H = db_H_table(good)
            for bad in (:v, :t, :u)
                dev = db_identity(cfg, bad, bid, 12; Href = H)
                bad === good ? (@test dev <= 1e-10) : (@test dev > 1e-3)
            end
        end
        # production order P = 4: the operator BUILDS and applies finitely on every block
        # with zero-data ghosts. The QUANTITATIVE production-order statement is gate (3d);
        # the identity above is gated at N = 12, where the degree-(6,4) probe is inside the
        # element space
        orb = db_orbit(cfg)
        dom = block_domain(orb, cfg)
        for bid in 0:9
            blk = dom.mesh.blocks[findfirst(b -> b.id == bid, dom.mesh.blocks)]
            els = filter(e -> e.block == bid, dom.mesh.elements)
            m1 = Mesh2D(blk, els, cfg.domain.initial_grid_points)
            bc = block_coefficients(dom, orb, cfg, bid)
            bcs = Dict{Tuple{ElementId,Int,Int},GhostBoundary}(
                k => DirichletGhost(zeros(ComplexF64, cfg.domain.initial_grid_points[1]))
                for k in boundary_faces(m1))
            op = MeshOperator(ComplexF64, m1, bc.A, bc.beta, bc.gamma;
                              C = cfg.dg.penalty_parameter, bcs = bcs)
            y = apply!(zeros(ComplexF64, m1.ndofs), op, ones(ComplexF64, m1.ndofs))
            @test all(isfinite, y)
            @test eltype(y) == ComplexF64
        end
    end

    @testset "(3d) P = 4: representable-degree exactness, the angular residual, h convergence" begin
        # The production-order statement, machine-checked. At N = 4 the element space is
        # degree ≤ 3 per dimension, so the N = 12 gate's degree-(6,4) probe is NOT in it and
        # the shipped instrument reports O(1) there (4.29 / 5.54 / 3.14 on blocks 0 / 2 / 4) —
        # that number measures the INSTRUMENT's exact-trace Dirichlet data (the ghost carries
        # the exact polynomial's face value while the element carries its degree-3
        # interpolant, and the difference rides the penalty term), not a defect of the
        # mapping. The three assertions below separate the channels; the calibration that
        # isolates them is evidence/domain-blocks-map/calibrate_followup2.log.
        cfg = crp_paper_config()
        N = cfg.domain.initial_grid_points[1]                       # 4
        for (region, bid) in ((:v, 0), (:t, 2), (:u, 4))
            # (i) a probe inside the space whose ANGULAR flux product also stays inside it
            #     (deg_w ≤ N−2, because A^{z²} = 4w(1−w)/(r²+a²) is degree 2 in w) is exact
            @test db_identity(cfg, region, bid, N; deg = (N - 1, N - 2)) <= 1e-12
            @test db_identity(cfg, region, bid, N; deg = (2, 2)) <= 1e-12
            # (ii) the RADIAL sector cancels identically — both sides write it as the
            #      collocation of the flux — so radial degree N−1 is exact too
            @test db_identity(cfg, region, bid, N; deg = (N - 1, 0)) <= 1e-12
            # (iii) what survives at representable degree is ANGULAR: one more degree in w
            #       pushes A^{z²}·∂_wΨ out of the N-point space
            @test db_identity(cfg, region, bid, N; deg = (0, N - 1)) > 1e-3
            @test db_identity(cfg, region, bid, N; deg = (0, N - 2)) <= 1e-12
        end
        # (iv) and the operator is CONSISTENT at the production order: at FIXED N = 4 the
        #      deviation h-converges under uniform refinement of the source column
        devs = [db_identity(cfg, :t, 2, N; k = kk) for kk in 0:3]
        rates = [log2(devs[i] / devs[i+1]) for i in 1:3]
        @test issorted(devs; rev = true)
        @test minimum(rates) >= 2
        @test devs[1] / devs[end] >= 100
        println("domain-blocks-map gate (3d) P = 4 h convergence, source column: devs = ",
                devs, "  rates = ", rates)
    end

    @testset "(3e) the :horizon precondition is ENFORCED, not assumed" begin
        # `:horizon` needs no boundary data ONLY because the inner bound is r₊, where the
        # sealed A^rr is exactly zero (main.tex:385 + 529). Away from r₊ it is not.
        for (a, rin) in ((0.5, 2.0), (-0.7, 2.0), (0.9, 1.3), (-0.9, 1.3))
            orb = CircularOrbit(a, 10.0)
            @test flux_coeff_r(orb, orb.rp) == 0                    # AT the horizon: exact
            @test abs(flux_coeff_r(orb, rin)) > 1e-2                # away from it: not
            bad = crp_paper_config(spin = a, inner_radius = rin)
            @test bad.domain.block_bounds_r[1] == rin
            @test_throws ArgumentError block_domain(orb, bad)       # refused at construction
        end
        # inside the horizon the coefficient even changes sign (ellipticity is lost)
        @test flux_coeff_r(CircularOrbit(0.9, 10.0), 1.3) < 0
        # the horizon-anchored twin of each of those configs builds
        for a in (0.5, -0.7, 0.9, -0.9)
            cfg = db_cfg(spin = a)
            orb = db_orbit(cfg)
            @test cfg.domain.block_bounds_r[1] == orb.rp
            dom = block_domain(orb, cfg)
            @test flux_coeff_r(orb, dom.bounds_r[1]) == 0
            # every :horizon face sits at r₊ with an identically vanishing radial flux
            for ((el, d, sidx), c) in boundary_face_classes(dom)
                c === :horizon || continue
                blk = dom.mesh.blocks[findfirst(b -> b.id == el.block, dom.mesh.blocks)]
                lo, hi = logical_interval(Float64, el, d)
                @test flux_coeff_r(orb, map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)) == 0
            end
        end
        # a tiny perturbation is still refused: the bound must BE the horizon
        orb = CircularOrbit(0.5, 10.0)
        @test_throws ArgumentError block_domain(orb,
            crp_paper_config(spin = 0.5, inner_radius = orb.rp * (1 + 1e-8)))
    end

    @testset "(4a) refinement levels L ∈ {0,1,2} and the L+1 source pre-refinement" begin
        for L in 0:2, sym in (true, false)
            cfg = crp_paper_config(L = L)
            sym || (cfg = db_full_z(cfg))
            dom = block_domain(db_orbit(cfg), cfg)
            nz = sym ? 2 : 3
            @test length(dom.mesh.elements) == nz * sum(db_nel_col(L, i) for i in 0:DB_NR-1)
            @test length(dom.mesh.elements) == nz * 10 * 4^L
            @test two_to_one_balanced(dom.mesh.elements)
            @test dom.mesh.elements == db_expected_elements(L, nz)
            # source column: L+1 in BOTH dimensions ⇒ 4^(L+1) elements per block
            wt = only(worldtube_blocks(dom))
            @test dom.refinements[wt] == (L + 1, L + 1)
            @test count(e -> e.block == wt, dom.mesh.elements) == 4^(L + 1)
            @test boundary_face_classes(dom) == db_expected_boundary(L, nz, dom.chart)
        end
    end

    @testset "(4b) parameter variation keeps every invariant" begin
        cases = (crp_paper_config(),
                 db_cfg(spin = 0.5),
                 db_cfg(spin = -0.7),
                 db_cfg(spin = 0.9),
                 db_cfg(spin = -0.9),
                 crp_paper_config(m_mode = 0),
                 crp_paper_config(m_mode = 20),
                 db_cfg(orbital_radius = 8.0, worldtube_radius = 2.0,
                        worldtube_dz2 = 0.125, outer_radius = 500.0, P = 6, L = 1),
                 db_full_z(db_cfg(spin = 0.5)),
                 db_full_z(db_cfg(spin = -0.7, orbital_radius = 8.0,
                                  worldtube_radius = 2.0, worldtube_dz2 = 0.36)))
        worst = 0.0
        for cfg in cases
            orb = db_orbit(cfg)
            dom = block_domain(orb, cfg)
            nz = cfg.orbit.impose_equatorial_symmetry ? 2 : 3
            rb, zb, _ = block_lattice_bounds(cfg)
            @test rb == cfg.domain.block_bounds_r
            @test issorted(rb) && allunique(rb)
            @test issorted(zb) && allunique(zb)
            @test length(dom.mesh.blocks) == DB_NR * nz
            @test two_to_one_balanced(dom.mesh.elements)
            @test Set(keys(boundary_face_classes(dom))) == Set(boundary_faces(dom.mesh))
            v = vtu_slicing(orb, cfg)
            for bid in 0:DB_NR*nz-1
                i, _ = block_lattice_index(dom, bid)
                @test block_region(dom, bid) ===
                      vtu_domain(v, (rb[i+1] + rb[i+2]) / 2)
            end
            @test length(worldtube_blocks(dom)) == 1
            # the identity survives the parameter change in every region
            off = cfg.orbit.impose_equatorial_symmetry ? 0 : DB_NR
            for (region, i) in ((:v, 0), (:t, 1), (:t, 2), (:t, 3), (:u, 4))
                dev = db_identity(cfg, region, db_bid(i, 0) + off, 12)
                worst = max(worst, dev)
                @test dev <= 1e-10
            end
        end
        println("domain-blocks-map gate (4b) worst identity rel dev over the sweep: ", worst)
    end

    @testset "(4c) m = 0 builds an identical lattice" begin
        c0, c2 = crp_paper_config(m_mode = 0), crp_paper_config(m_mode = 2)
        d0, d2 = block_domain(db_orbit(c0), c0), block_domain(db_orbit(c2), c2)
        @test d0.mesh.elements == d2.mesh.elements
        @test d0.mesh.ndofs == d2.mesh.ndofs
        @test d0.bounds_r == d2.bounds_r && d0.bounds_z == d2.bounds_z
        @test d0.regions == d2.regions
        @test worldtube_blocks(d0) == worldtube_blocks(d2)
        @test boundary_face_classes(d0) == boundary_face_classes(d2)
        @test c0.domain.outer_radial_bc === nothing      # node 41's business, not the mesh's
        @test c2.domain.outer_radial_bc isa SommerfeldBC
    end

    @testset "(4d) BigFloat(256) precision inheritance" begin
        setprecision(BigFloat, 256) do
            rp = CircularOrbit(big"0.5", big"10.0").rp
            cfg = crp_paper_config(spin = big"0.5", orbital_radius = big"10.0",
                                   inner_radius = rp, outer_radius = big"1.0e4",
                                   worldtube_radius = big"4.8", worldtube_dz2 = big"0.25")
            @test cfg isa RunConfig{BigFloat}
            orb = db_orbit(cfg)
            dom = block_domain(orb, cfg)
            @test dom isa BlockDomain{BigFloat}
            @test dom.mesh isa Mesh2D{BigFloat}
            @test eltype(dom.bounds_r) == BigFloat && eltype(dom.bounds_z) == BigFloat
            @test all(b -> precision(b) == 256, dom.bounds_r)
            @test dom.bounds_r == cfg.domain.block_bounds_r
            ublk = dom.mesh.blocks[findfirst(b -> b.id == 4, dom.mesh.blocks)]
            mx = blockmap(ublk, 1)
            @test mx isa InverseRadialMap1D{BigFloat}
            @test map_forward(mx, big"-1.0") == dom.bounds_r[5]     # endpoint-exact
            @test map_forward(mx, big"1.0") == dom.bounds_r[6]
            @test precision(map_forward(mx, big"0.3")) == 256
            # the non-equatorial chart's z rows inherit precision through the sqrt
            domz = block_domain(orb, db_full_z(cfg))
            @test eltype(domz.bounds_z) == BigFloat
            @test domz.bounds_z[3]^2 == cfg.domain.block_bounds_z2[2]
            @test precision(domz.bounds_z[3]) == 256
            # structure is precision-independent
            @test domz.mesh.elements ==
                  block_domain(db_orbit(crp_paper_config()),
                               db_full_z(crp_paper_config())).mesh.elements
        end
    end

    @testset "(5) construction validation refuses inconsistent configurations" begin
        cfg = crp_paper_config()
        orb = db_orbit(cfg)
        # orbit ≠ config point (the sealed vtu_slicing contract)
        @test_throws ArgumentError block_domain(CircularOrbit(0.5, 10.0), cfg)
        @test_throws ArgumentError block_domain(CircularOrbit(0.0, 8.0), cfg)
        # slicing transitions that are NOT block bounds would put two regions in one block
        bad = RunConfig{Float64}((f === :orbit ?
                typeof(cfg.orbit)((g === :hyperboloidal_slicing_transitions ?
                    (3.0, 3.0, 24.4, 24.4) : getfield(cfg.orbit, g)
                    for g in fieldnames(typeof(cfg.orbit)))...) : getfield(cfg, f)
                for f in fieldnames(RunConfig{Float64}))...)
        @test_throws ArgumentError block_domain(orb, bad)
        # a radial Distributions column that is not [Linear ×(n−1), Inverse]
        baddist = RunConfig{Float64}((f === :domain ?
                typeof(cfg.domain)((g === :distributions_r ?
                    [:Linear, :Linear, :Inverse, :Linear, :Inverse] : getfield(cfg.domain, g)
                    for g in fieldnames(typeof(cfg.domain)))...) : getfield(cfg, f)
                for f in fieldnames(RunConfig{Float64}))...)
        @test_throws ArgumentError block_domain(orb, baddist)
        # per-element p refinement (RefinedGridPoints) is not this node's scope
        badp = RunConfig{Float64}((f === :domain ?
                typeof(cfg.domain)((g === :refined_grid_points ? [(6, 6)] :
                    getfield(cfg.domain, g)
                    for g in fieldnames(typeof(cfg.domain)))...) : getfield(cfg, f)
                for f in fieldnames(RunConfig{Float64}))...)
        @test_throws ArgumentError block_domain(orb, badp)
        # unknown block id
        # a RefinedLevels patch outside the lattice in EITHER dimension (spectre PARSE_ERRORs
        # on both: AlignedLattice.cpp:123-131 tests upper_corner_index[i] >= block_bounds[i]
        # size for every i). The radial overflow was already refused; the angular one is the
        # symmetric check.
        for badpatch in (BlockRefinement((0, 1), (1, 3), (0, 0)),      # row 2 of a 2-row chart
                         BlockRefinement((4, 0), (6, 1), (0, 0)))      # column 5 of 5
            badref = RunConfig{Float64}((f === :domain ?
                    typeof(cfg.domain)((g === :refined_levels ?
                        vcat(cfg.domain.refined_levels, [badpatch]) : getfield(cfg.domain, g)
                        for g in fieldnames(typeof(cfg.domain)))...) : getfield(cfg, f)
                    for f in fieldnames(RunConfig{Float64}))...)
            @test_throws ArgumentError block_domain(orb, badref)
        end
        # overlapping RefinedLevels regions that disagree are an error here, NOT spectre's
        # silent last-write-wins (deliberate divergence, documented in blocks.jl)
        clash = RunConfig{Float64}((f === :domain ?
                typeof(cfg.domain)((g === :refined_levels ?
                    vcat(cfg.domain.refined_levels,
                         [BlockRefinement((2, 0), (3, 1), (7, 7))]) : getfield(cfg.domain, g)
                    for g in fieldnames(typeof(cfg.domain)))...) : getfield(cfg, f)
                for f in fieldnames(RunConfig{Float64}))...)
        @test_throws ArgumentError block_domain(orb, clash)
        dom = block_domain(orb, cfg)
        @test_throws KeyError block_region(dom, 99)
        @test_throws ArgumentError block_lattice_index(dom, 99)
    end
end
