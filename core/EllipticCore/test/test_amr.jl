# test_amr.jl — grounded gates for src/amr.jl (DAG nodes: amr-h-refinement GLOBAL_DAG row 22
# and amr-p-refinement row 23; elliptic-ssf TOPOLOGICAL_SCHEDULE W7-B). Written BEFORE the
# implementation (TDD red runs captured in
# evidence/amr-{h,p}-refinement/tests_red.log).
#
# Target: the CRP hp-AMR mesh transformations (arXiv:2606.04998 main.tex:552-560, 714):
#   h — "split all elements containing the puncture in half in both dimensions in each
#       refinement step … when splitting elements in half, we retain a two-to-one balance at
#       element interfaces … in cases where the two-to-one balance would be violated, the
#       offending elements are split in half as well" (main.tex:552-554);
#   p — "begin with 4 x 4 grid points per element initially, and then increase the polynomial
#       order in all elements that do not contain the puncture by one in both dimensions in
#       each refinement step" (main.tex:558-560);
#   and the Lagrange carry-over of the current solution onto the new mesh.
# The refinement CRITERION and the solve/AMR driver loop are NOT here — criteria are
# caller-supplied predicates (this file uses synthetic ones) and the loop is the
# `amr-driver` node (row 43).
#
# Oracle gates (worker node contract; constants calibrated PRE-RED in the evidence dir:
# calibrate_amr.jl, calibrate_balance.jl, calibrate_mms.jl — throwaway prototypes, NOT the
# implementation):
#   (a) flag forms (Set / Dict / predicate / per-dimension tuple) agree; children inherit the
#       parent order; an empty flag set is the identity transformation.
#   (b) the paper's six-step centred-puncture cascade reproduces a CLOSED-FORM expected element
#       set computed independently here (1, 4, 16, 28, 40, 52, 64 elements = 12k - 8 for k >= 2,
#       with 12 elements per level 2..k-1 and 16 at level k), and a six-step OFF-CENTRE cascade
#       on a two-block layout reproduces a NAIVE independent re-implementation of split +
#       balance (Float64 interval arithmetic, O(n^2) fixpoint) at every step — that cascade
#       fires the balancer on every step, including ACROSS the block connection.
#   (c) exhaustive random gate: 150 random block layouts x 4 rounds of random per-dimension
#       flag sets — after enforcement `two_to_one_balanced` holds, the cross-block transverse
#       balance holds, the balancer is IDEMPOTENT (a second call returns the identical vector),
#       the result is a leaf tiling of the input (exact measure conservation), and the mesh is
#       constructible through the PUBLIC Mesh2D constructor. Falsifier: the UNBALANCED split
#       output violates 2:1 in a large fraction of the trials (the balancer is load-bearing).
#   (d) carry-over: EXACT (<= 1e-13, calibrated 4.4e-16) on polynomials of the element degree
#       for every (Delta-level, child-index, order) combination reached; on a smooth function
#       the error is <= the Gauss interpolation-theory bound H^N 2^N N!/(2N)! max|f^(N)| in x
#       plus the Lebesgue-weighted y term (calibrated worst err/bound 0.179), and decays by
#       2^-N per level (calibrated within [0.9, 1.2] x 2^-N); the one-level same-order transfer
#       is BITWISE the sealed mortar prolongation `projection_prolong(T, N, N, MORTAR_LOWER/
#       UPPER)` (no duplicated interpolation semantics); ComplexF64 linearity; coarsening is a
#       refused ArgumentError.
#   (e) the transformed mesh satisfies every Mesh2D invariant and carries a MeshOperator built
#       with the SAME coefficient closures: apply! is finite and linear, and the manufactured
#       Poisson L2 error drops by 0.072 (gate < 0.15) from the coarse to the h-refined mesh.
#   (f) purity: the source Mesh2D and the source solution vector are untouched, and the new
#       mesh shares no mutable state with the old one.
#   (g) Float64 and BigFloat(256) carry-over (BigFloat polynomial residual calibrated 3.5e-77).
#
# Node 23 (p) gates:
#   (h) metadata: after k = 1..8 p steps with a FIXED flag set, N = 4 + k in every non-flagged
#       element and N = 4 in every flagged one (main.tex:558-560); per-dimension flags hold only
#       the flagged dimension; `increment` is honoured.
#   (i) p transfer is exact on polynomials of the SOURCE element degree (<= 1e-13) and is the
#       identity when the order does not change (bitwise).
#   (j) the combined h+p step `hp_refine` == p_refine-then-h_refine bookkeeping: flagged
#       elements h-split and keep N = 4, everything else keeps its element and gains one order;
#       balance-induced children inherit their parent's (already incremented) order.
#   (k) the p-transformed mesh is valid: a MeshOperator on the p-nonconforming result builds and
#       applies (sealed mortars), and the manufactured Poisson L2 error drops by 0.071
#       (gate < 0.15) for a uniform p step.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists (each
# @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as the sibling test
# ---- files, distinct name/seed) ---------------------------------------------------------
mutable struct AMRRNG
    s::UInt64
end
amr_next!(r::AMRRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f; r.s >> 11)
amr_u01!(r::AMRRNG) = Float64(amr_next!(r)) / 2.0^53
amr_pick!(r::AMRRNG, n) = Int(amr_next!(r) % n) + 1
amr_cvec!(r::AMRRNG, n) = [complex(2 * amr_u01!(r) - 1, 2 * amr_u01!(r) - 1) for _ in 1:n]

const AMR_C = 1.5                                  # CRP main.tex:595 penalty parameter

# ---- INDEPENDENT (naive, obviously-correct) refinement-tree arithmetic, written only from
# ---- the dyadic definition of SegmentId — no `src/amr.jl` and no `src/mesh2d.jl` predicate
amr_ivl(s) = (s.index / 2.0^s.level, (s.index + 1) / 2.0^s.level)
amr_lv(id, d) = id.segments[d].level
function amr_touch(a, b)
    x = amr_ivl(a); y = amr_ivl(b)
    return x[2] == y[1] || y[2] == x[1]
end
function amr_ovl(a, b)
    x = amr_ivl(a); y = amr_ivl(b)
    return max(x[1], y[1]) < min(x[2], y[2])
end
amr_adjacent(a, b) =
    a.block == b.block &&
    ((amr_touch(a.segments[1], b.segments[1]) && amr_ovl(a.segments[2], b.segments[2])) ||
     (amr_touch(a.segments[2], b.segments[2]) && amr_ovl(a.segments[1], b.segments[1])))

# naive split of one element along the flagged dimensions (child order: x fastest)
function amr_split1(id, dims)
    out = [id]
    for d in 1:2
        dims[d] || continue
        out = vcat([[child_id(x, d, 0), child_id(x, d, 1)] for x in out]...)
    end
    return out
end

# naive O(n^2) 2:1 balance fixpoint (independent re-implementation of the gate (c) oracle)
function amr_naive_balance(ids, conns)
    cur = collect(ids)
    sweeps = 0
    while true
        req = Dict{ElementId,Vector{Bool}}()
        ask!(id, d) = (v = get!(req, id, [false, false]); v[d] = true)
        n = length(cur)
        for i in 1:n-1, j in i+1:n
            a, b = cur[i], cur[j]
            amr_adjacent(a, b) || continue
            for d in 1:2
                δ = amr_lv(a, d) - amr_lv(b, d)
                δ >= 2 && ask!(b, d)
                δ <= -2 && ask!(a, d)
            end
        end
        for c in conns                                    # cross-block: transverse only
            d = c.dim; t = 3 - d
            for a in cur, b in cur
                (a.block == c.block && b.block == c.neighbor) || continue
                a.segments[d].index == (1 << a.segments[d].level) - 1 || continue
                b.segments[d].index == 0 || continue
                amr_ovl(a.segments[t], b.segments[t]) || continue
                δ = amr_lv(a, t) - amr_lv(b, t)
                δ >= 2 && ask!(b, t)
                δ <= -2 && ask!(a, t)
            end
        end
        isempty(req) && return cur, sweeps
        sweeps += 1
        sweeps > 200 && error("naive balance did not terminate")
        nxt = ElementId[]
        for id in cur
            haskey(req, id) ? append!(nxt, amr_split1(id, Tuple(req[id]))) : push!(nxt, id)
        end
        cur = nxt
    end
end

# cross-block transverse 2:1, checked independently
function amr_crossblock_ok(ids, conns)
    for c in conns
        d = c.dim; t = 3 - d
        for a in ids, b in ids
            (a.block == c.block && b.block == c.neighbor) || continue
            a.segments[d].index == (1 << a.segments[d].level) - 1 || continue
            b.segments[d].index == 0 || continue
            amr_ovl(a.segments[t], b.segments[t]) || continue
            abs(amr_lv(a, t) - amr_lv(b, t)) <= 1 || return false
        end
    end
    return true
end

# is `anc` an ancestor-or-equal of `id`?
function amr_isanc(anc, id)
    anc.block == id.block || return false
    for d in 1:2
        sa = anc.segments[d]; si = id.segments[d]
        si.level >= sa.level || return false
        (si.index >> (si.level - sa.level)) == sa.index || return false
    end
    return true
end
amr_area(id) = 1 / 2.0^(amr_lv(id, 1) + amr_lv(id, 2))

# closed logical containment of a point (px, py) in [0,1]^2 block-logical units
amr_contains(id, px, py) = begin
    a = amr_ivl(id.segments[1]); b = amr_ivl(id.segments[2])
    a[1] <= px <= a[2] && b[1] <= py <= b[2]
end

# ---- CLOSED-FORM expected centred-puncture cascade (independent of any implementation) ----
# C(l) = the four level-l elements whose closed rectangle contains the block centre.
amr_centre(blk, l) = l == 0 ? [ElementId(blk, SegmentId(0, 0), SegmentId(0, 0))] :
    [ElementId(blk, SegmentId(l, ix), SegmentId(l, iy))
     for iy in (2^(l-1) - 1, 2^(l-1)), ix in (2^(l-1) - 1, 2^(l-1))]
amr_kids(ids) = vcat([amr_split1(id, (true, true)) for id in ids]...)
function amr_expected_cascade(blk, k)
    k == 0 && return [ElementId(blk, SegmentId(0, 0), SegmentId(0, 0))]
    out = ElementId[]
    for l in 2:k-1                              # 12 non-centre children per intermediate level
        append!(out, setdiff(amr_kids(amr_centre(blk, l - 1)), amr_centre(blk, l)))
    end
    append!(out, amr_kids(amr_centre(blk, k - 1)))          # 16 at the finest level
    return out
end

# ---- fixtures ------------------------------------------------------------------------------
amr_blk1() = Block(21, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
amr_root1() = ElementId(21, SegmentId(0, 0), SegmentId(0, 0))
amr_blocks2() = [Block(1, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0)),
                 Block(2, LinearMap1D(1.0, 2.0), LinearMap1D(0.0, 1.0))]
amr_conns2() = [BlockConnection(1, 1, 2)]
amr_roots2() = [ElementId(1, SegmentId(0, 0), SegmentId(0, 0)),
                ElementId(2, SegmentId(0, 0), SegmentId(0, 0))]

amr_one2(x, y) = 1.0
amr_zero2(x, y) = 0.0
amr_uex(x, y) = sinpi(x) * sinpi(y)                   # zero Dirichlet on [0,1]^2
amr_src(x, y) = 2 * pi^2 * amr_uex(x, y)              # -div grad u = S
amr_A11(x, y) = 1 + 0.3x + 0.1y^2
amr_A22(x, y) = 1 + 0.2y
amr_beta(x, y) = 0.3 + 0.2im

function amr_zerobcs(mesh)
    b = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        b[(el, d, sidx)] = DirichletGhost(zeros(ComplexF64, mesh.orders[mesh.index[el]][3-d]))
    end
    return b
end

# nodal sample of f on a mesh, in the mesh's global layout
function amr_sample(mesh, blk, f, ::Type{T} = Float64) where {T}
    v = Vector{T}(undef, mesh.ndofs)
    for (e, id) in enumerate(mesh.elements)
        em = element_matrices(real(T), blk, id, mesh.orders[e])
        Nx, Ny = em.N
        v[mesh.offsets[e]+1:mesh.offsets[e]+Nx*Ny] .=
            vec(T[f(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
    end
    return v
end

# mass-weighted L2 error of a nodal vector against an exact function
function amr_l2err(mesh, blk, u, f)
    s = 0.0
    for (e, id) in enumerate(mesh.elements)
        em = element_matrices(Float64, blk, id, mesh.orders[e])
        Nx, Ny = em.N
        d = u[mesh.offsets[e]+1:mesh.offsets[e]+Nx*Ny] .-
            vec([f(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
        s += real(dot(d, em.M * d))
    end
    return sqrt(s)
end

# manufactured Poisson solve on a mesh (direct bypass) -> L2 error
function amr_poisson_l2(mesh, blk)
    op = MeshOperator(ComplexF64, mesh, (amr_one2, amr_one2), amr_zero2,
                      (amr_zero2, amr_zero2); C = AMR_C, bcs = amr_zerobcs(mesh))
    b = zeros(ComplexF64, mesh.ndofs)
    for (e, id) in enumerate(mesh.elements)
        em = element_matrices(Float64, blk, id, mesh.orders[e])
        Nx, Ny = em.N
        b[mesh.offsets[e]+1:mesh.offsets[e]+Nx*Ny] .=
            em.M * vec([amr_src(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
    end
    x = direct_solve(op, b .- operator_on_zero(op))
    return amr_l2err(mesh, blk, x, amr_uex)
end

# Gauss interpolation-error bound on a half-width-H interval with N nodes:
#   |f - I f| <= H^N 2^N N!/(2N)! max|f^(N)|
amr_gauss_bound(H, N, dmax) =
    Float64((big(H)^N) * (big(2.0)^N) * factorial(big(N)) / factorial(big(2N)) * dmax)

# ============================================================================================
@testset "AMR h/p mesh transformations (amr-h-refinement, amr-p-refinement)" begin

# ---- (a) flag forms, order inheritance, identity -------------------------------------------
@testset "(a) flag forms, child order inheritance, identity transformation" begin
    blk = amr_blk1()
    ids = [ElementId(21, SegmentId(1, i), SegmentId(1, j)) for j in 0:1 for i in 0:1]
    ords = Dict(ids[1] => (4, 4), ids[2] => (5, 5), ids[3] => (4, 6), ids[4] => (7, 4))
    mesh = Mesh2D(blk, ids, ords)
    target = ids[1]

    m_set  = h_refine(mesh, Set([target]))
    m_dict = h_refine(mesh, Dict(target => true))
    m_pred = h_refine(mesh, id -> id == target)
    m_tup  = h_refine(mesh, Dict(target => (true, true)))
    m_vec  = h_refine(mesh, [target])
    @test m_set.elements == m_dict.elements == m_pred.elements == m_tup.elements == m_vec.elements
    @test m_set.orders == m_dict.orders == m_pred.orders == m_tup.orders == m_vec.orders

    # one flagged element in both dimensions -> 4 children, 3 untouched siblings
    @test length(m_set.elements) == 7
    @test count(id -> amr_isanc(target, id) && id != target, m_set.elements) == 4
    # children inherit the parent's order; untouched elements keep theirs
    for (e, id) in enumerate(m_set.elements)
        anc = amr_isanc(target, id) ? target : id
        @test m_set.orders[e] == ords[anc]
    end

    # per-dimension flags: split dimension 1 only -> 2 children at (level+1, level)
    m_x = h_refine(mesh, Dict(target => (true, false)))
    @test length(m_x.elements) == 5
    kids = filter(id -> amr_isanc(target, id) && id != target, m_x.elements)
    @test length(kids) == 2
    @test all(amr_lv(k, 1) == 2 && amr_lv(k, 2) == 1 for k in kids)
    m_y = h_refine(mesh, Dict(target => (false, true)))
    kids_y = filter(id -> amr_isanc(target, id) && id != target, m_y.elements)
    @test length(kids_y) == 2
    @test all(amr_lv(k, 1) == 1 && amr_lv(k, 2) == 2 for k in kids_y)

    # empty flag set is the identity transformation (elements, orders, layout)
    m_id = h_refine(mesh, Set{ElementId}())
    @test m_id.elements == mesh.elements
    @test m_id.orders == mesh.orders
    @test m_id.offsets == mesh.offsets && m_id.ndofs == mesh.ndofs

    # h_split is the pure address-level transformation (no balance, no mesh)
    @test h_split(ids, Set([target])) == m_set.elements
    @test h_split(ids, Set{ElementId}()) == ids
    # deterministic and reproducible
    @test h_refine(mesh, Set([target])).elements == m_set.elements
end

# ---- (b1) the paper's centred-puncture cascade, six steps -----------------------------------
@testset "(b1) six-step centred-puncture cascade == closed-form expected cascade" begin
    blk = amr_blk1()
    mesh = Mesh2D(blk, [amr_root1()], (4, 4))
    punct(id) = amr_contains(id, 0.5, 0.5)          # the synthetic flag rule (main.tex:552)
    expected_counts = [4, 16, 28, 40, 52, 64]       # 12k - 8 for k >= 2
    for k in 1:6
        mesh = h_refine(mesh, punct)
        exp_ids = amr_expected_cascade(21, k)
        @test length(mesh.elements) == expected_counts[k]
        @test Set(mesh.elements) == Set(exp_ids)
        # the LEVEL MAP (refinement level per ElementId) matches element by element
        got = Dict(id => (amr_lv(id, 1), amr_lv(id, 2)) for id in mesh.elements)
        want = Dict(id => (amr_lv(id, 1), amr_lv(id, 2)) for id in exp_ids)
        @test got == want
        @test all(v -> v[1] == v[2], values(got))                   # isotropic cascade
        # per-level histogram: 12 at each level 2..k-1, 16 at level k (k >= 2)
        if k >= 2
            hist = Dict{Int,Int}()
            for id in mesh.elements
                hist[amr_lv(id, 1)] = get(hist, amr_lv(id, 1), 0) + 1
            end
            @test hist[k] == 16
            @test all(hist[l] == 12 for l in 2:k-1)
        end
        @test two_to_one_balanced(mesh.elements)
        @test isapprox(sum(amr_area, mesh.elements), 1.0; atol = 1e-14)   # exact tiling
        @test all(o -> o == (4, 4), mesh.orders)                          # h keeps the order
    end
end

# ---- (b2) off-centre cascade on a two-block layout vs the naive independent oracle ---------
@testset "(b2) six-step off-centre two-block cascade == naive independent balance" begin
    blocks = amr_blocks2(); conns = amr_conns2()
    mesh = Mesh2D(blocks, amr_roots2(), (4, 4); connections = conns)
    naive = copy(amr_roots2())
    punct(id) = id.block == 1 && amr_contains(id, 0.87, 0.31)
    balance_fired = 0
    crossblock_splits = 0
    for k in 1:6
        mesh = h_refine(mesh, punct)
        raw = ElementId[]
        for id in naive
            append!(raw, punct(id) ? amr_split1(id, (true, true)) : [id])
        end
        nb, _ = amr_naive_balance(raw, conns)
        balance_fired += (length(nb) - length(raw)) > 0
        crossblock_splits += count(id -> id.block == 2, nb) - 1
        @test Set(mesh.elements) == Set(nb)
        @test two_to_one_balanced(mesh.elements)
        @test amr_crossblock_ok(mesh.elements, conns)
        naive = nb
    end
    # step 1 only halves the block-1 root (level 1 vs the block-2 root's level 0 — already
    # balanced); the balancer fires on each of the remaining five steps
    @test balance_fired == 5
    @test crossblock_splits > 0         # ... including ACROSS the block connection
end

# ---- (c) exhaustive random gate: balance, idempotence, tiling, constructibility -------------
@testset "(c) random layouts x random flag sets: 2:1 enforcement is exhaustive" begin
    rng = AMRRNG(0x2606049823456789)
    bad_balance = String[]; bad_cross = String[]; bad_idem = String[]
    bad_tiling = String[]; bad_mesh = String[]
    unbalanced_raw = 0; balance_fired = 0; rounds = 0; maxels = 0
    for trial in 1:150
        nb = amr_pick!(rng, 3)
        blocks = [Block(k, LinearMap1D(Float64(k - 1), Float64(k)), LinearMap1D(0.0, 1.0))
                  for k in 1:nb]
        conns = [BlockConnection(k, 1, k + 1) for k in 1:nb-1]
        mesh = Mesh2D(blocks, [ElementId(k, SegmentId(0, 0), SegmentId(0, 0)) for k in 1:nb],
                      (4, 4); connections = conns)
        for round in 1:4
            flags = Dict{ElementId,NTuple{2,Bool}}()
            for id in mesh.elements
                amr_u01!(rng) < 0.35 &&
                    (flags[id] = (amr_u01!(rng) < 0.8, amr_u01!(rng) < 0.8))
            end
            old = copy(mesh.elements)
            raw = h_split(old, flags)
            two_to_one_balanced(raw) || (unbalanced_raw += 1)
            local newmesh
            try
                newmesh = h_refine(mesh, flags)
            catch err
                push!(bad_mesh, "trial $trial round $round: $err")
                break
            end
            rounds += 1
            maxels = max(maxels, length(newmesh.elements))
            balance_fired += length(newmesh.elements) > length(raw)
            two_to_one_balanced(newmesh.elements) ||
                push!(bad_balance, "trial $trial round $round")
            amr_crossblock_ok(newmesh.elements, conns) ||
                push!(bad_cross, "trial $trial round $round")
            # idempotence: a second enforcement returns the IDENTICAL vector (value + order)
            balance_two_to_one(newmesh.elements; connections = conns) == newmesh.elements ||
                push!(bad_idem, "trial $trial round $round")
            # leaf tiling of the input: unique ancestor per element + exact measure
            ok = all(count(a -> amr_isanc(a, id), old) == 1 for id in newmesh.elements) &&
                 abs(sum(amr_area, newmesh.elements) - sum(amr_area, old)) < 1e-12
            ok || push!(bad_tiling, "trial $trial round $round")
            mesh = newmesh
            length(mesh.elements) > 400 && break
        end
    end
    @test isempty(bad_balance)
    @test isempty(bad_cross)
    @test isempty(bad_idem)
    @test isempty(bad_tiling)
    @test isempty(bad_mesh)
    @test rounds >= 500                       # the gate is not vacuous
    @test balance_fired > 100                 # enforcement really fires
    @test unbalanced_raw > 100                # FALSIFIER: the raw split IS 2:1-violating
end

# ---- (d) solution carry-over -----------------------------------------------------------------
@testset "(d1) carry-over is EXACT on polynomials of the element degree" begin
    blk = amr_blk1()
    worst = 0.0
    for N in (4, 6)
        # start from a mixed-level mesh so that Delta = 1 and Delta = 2 both occur
        mesh = Mesh2D(blk, [amr_root1()], (N, N))
        mesh = h_refine(mesh, id -> true)                       # uniform -> level 1
        mesh = h_refine(mesh, id -> amr_contains(id, 0.5, 0.5)) # centre -> level 2
        # a polynomial of degree N-1 in each variable is represented EXACTLY on every element
        poly(x, y) = sum(((0.7)^k) * x^k for k in 0:N-1) * sum(((-0.4)^k) * y^k for k in 0:N-1)
        u = amr_sample(mesh, blk, poly)
        for flagrule in (id -> true,
                         id -> amr_contains(id, 0.5, 0.5),
                         id -> amr_lv(id, 1) == 1)
            new = h_refine(mesh, flagrule)
            v = transfer_solution(mesh, new, u)
            @test length(v) == new.ndofs
            worst = max(worst, maximum(abs.(v .- amr_sample(new, blk, poly))))
        end
    end
    @test worst < 1e-13                       # calibrated 4.4e-16
end

@testset "(d2) smooth carry-over <= the Gauss interpolation-theory bound" begin
    blk = amr_blk1()
    a, b, ph = 0.7, 1.3, 0.4
    fsm(x, y) = exp(a * x) * cos(b * y + ph)
    ratios = Float64[]
    for N in (4, 5, 6), lev in (0, 1)
        parent = ElementId(21, SegmentId(lev, 0), SegmentId(lev, 0))
        mesh = Mesh2D(blk, [ElementId(21, SegmentId(lev, i), SegmentId(lev, j))
                            for j in 0:(2^lev - 1) for i in 0:(2^lev - 1)], (N, N))
        u = amr_sample(mesh, blk, fsm)
        new = h_refine(mesh, id -> id == parent)
        v = transfer_solution(mesh, new, u)
        err = maximum(abs.(v .- amr_sample(new, blk, fsm)))
        H = 0.5 / 2.0^lev                                   # parent physical half-width
        dxmax = a^N * exp(a * 2H); dymax = b^N * exp(a * 2H)
        Ex = amr_gauss_bound(H, N, dxmax); Ey = amr_gauss_bound(H, N, dymax)
        Λ = maximum(sum(abs, projection_prolong(Float64, N, N, MORTAR_LOWER); dims = 2))
        bound = Ex + Λ * Ey
        @test err <= bound                                  # interpolation theory
        push!(ratios, err / bound)
    end
    @test maximum(ratios) <= 0.5                            # calibrated worst 0.179
    @test minimum(ratios) >= 1e-3                           # not vacuously tiny

    # error decays by 2^-N per refinement level (calibrated within [0.9, 1.2] x 2^-N)
    for N in (4, 6)
        prev = NaN
        for lev in 0:3
            mesh = Mesh2D(blk, [ElementId(21, SegmentId(lev, i), SegmentId(lev, j))
                                for j in 0:(2^lev - 1) for i in 0:(2^lev - 1)], (N, N))
            u = amr_sample(mesh, blk, fsm)
            new = h_refine(mesh, id -> true)
            err = maximum(abs.(transfer_solution(mesh, new, u) .- amr_sample(new, blk, fsm)))
            lev >= 2 && @test 0.9 <= (err / prev) / 2.0^-N <= 1.2
            prev = err
        end
    end
end

@testset "(d3) one-level same-order carry-over IS the sealed mortar prolongation" begin
    blk = amr_blk1()
    for N in (4, 5, 6, 8)
        mesh = Mesh2D(blk, [amr_root1()], (N, N))
        new = h_refine(mesh, id -> true)
        @test length(new.elements) == 4
        Plo = projection_prolong(Float64, N, N, MORTAR_LOWER)
        Phi = projection_prolong(Float64, N, N, MORTAR_UPPER)
        u = [Float64(k) for k in 1:N*N]                    # arbitrary nodal data
        v = transfer_solution(mesh, new, u)
        U = reshape(u, N, N)
        for (e, id) in enumerate(new.elements)
            Px = id.segments[1].index == 0 ? Plo : Phi
            Py = id.segments[2].index == 0 ? Plo : Phi
            want = vec(Px * U * transpose(Py))
            got = v[new.offsets[e]+1:new.offsets[e]+N*N]
            @test got == want                              # BITWISE
        end
    end
end

@testset "(d4) complex linearity, refused coarsening, layout validation" begin
    blk = amr_blk1()
    mesh = Mesh2D(blk, [amr_root1()], (5, 5))
    new = h_refine(mesh, id -> true)
    rng = AMRRNG(0xa11cea5e)
    u = amr_cvec!(rng, mesh.ndofs); w = amr_cvec!(rng, mesh.ndofs)
    α = 0.3 - 1.7im
    lhs = transfer_solution(mesh, new, α .* u .+ w)
    rhs = α .* transfer_solution(mesh, new, u) .+ transfer_solution(mesh, new, w)
    @test maximum(abs.(lhs .- rhs)) <= 1e-13 * maximum(abs.(rhs))
    @test eltype(transfer_solution(mesh, new, u)) == ComplexF64
    # identity transfer (same mesh) is bitwise the input
    @test transfer_solution(mesh, mesh, u) == u
    # coarsening / join is out of scope: refused, not silently wrong
    @test_throws ArgumentError transfer_solution(new, mesh, transfer_solution(mesh, new, u))
    @test_throws ArgumentError transfer_solution(mesh, new, u[1:end-1])
end

# ---- (e) Mesh2D invariants + MeshOperator on the transformed mesh ---------------------------
@testset "(e) transformed mesh is a valid Mesh2D carrying a MeshOperator" begin
    blocks = amr_blocks2(); conns = amr_conns2()
    mesh = Mesh2D(blocks, amr_roots2(), (4, 4); connections = conns)
    for k in 1:3
        mesh = h_refine(mesh, id -> id.block == 1 && amr_contains(id, 0.87, 0.31))
    end
    # Mesh2D invariants, re-derived here from the public fields
    @test length(mesh.elements) == length(mesh.orders) == length(mesh.offsets)
    @test mesh.index == Dict(id => e for (e, id) in enumerate(mesh.elements))
    @test mesh.offsets[1] == 0
    @test all(mesh.offsets[e] == mesh.offsets[e-1] + prod(mesh.orders[e-1])
              for e in 2:length(mesh.elements))
    @test mesh.ndofs == mesh.offsets[end] + prod(mesh.orders[end])
    @test length(mesh.faces) == length(mesh.elements)
    @test mesh.connections == conns
    @test !isempty(boundary_faces(mesh))
    # ... and it is re-constructible through the PUBLIC constructor
    twin = Mesh2D(blocks, mesh.elements,
                  Dict(id => mesh.orders[e] for (e, id) in enumerate(mesh.elements));
                  connections = conns)
    @test twin.elements == mesh.elements && twin.ndofs == mesh.ndofs

    # a MeshOperator with the SAME coefficient closures builds and applies
    op = MeshOperator(ComplexF64, mesh, (amr_A11, amr_A22), amr_beta,
                      (amr_zero2, amr_zero2); C = AMR_C, bcs = amr_zerobcs(mesh))
    rng = AMRRNG(0x0ffee123)
    u = amr_cvec!(rng, mesh.ndofs); w = amr_cvec!(rng, mesh.ndofs)
    y = similar(u)
    apply!(y, op, u)
    @test all(isfinite, y)
    y2 = similar(u); y3 = similar(u)
    apply!(y2, op, w); apply!(y3, op, u .+ w)
    @test maximum(abs.(y3 .- y .- y2)) <= 1e-12 * maximum(abs.(y3))

    # manufactured Poisson: h refinement reduces the L2 error (calibrated ratio 0.072)
    blk = amr_blk1()
    coarse = Mesh2D(blk, [ElementId(21, SegmentId(1, i), SegmentId(1, j))
                          for j in 0:1 for i in 0:1], (4, 4))
    fine = h_refine(coarse, id -> true)
    ec = amr_poisson_l2(coarse, blk); ef = amr_poisson_l2(fine, blk)
    @test ef < 0.15 * ec
end

# ---- (f) purity ------------------------------------------------------------------------------
@testset "(f) the transformations are pure functions (source mesh untouched)" begin
    blocks = amr_blocks2(); conns = amr_conns2()
    mesh = Mesh2D(blocks, amr_roots2(), (4, 4); connections = conns)
    mesh = h_refine(mesh, id -> true)
    els0 = copy(mesh.elements); ords0 = copy(mesh.orders); offs0 = copy(mesh.offsets)
    idx0 = copy(mesh.index); nd0 = mesh.ndofs; conn0 = copy(mesh.connections)
    faces0 = copy(mesh.faces)
    u = collect(1.0:mesh.ndofs); u0 = copy(u)

    new = h_refine(mesh, id -> id.block == 1)
    v = transfer_solution(mesh, new, u)

    @test mesh.elements == els0 && mesh.orders == ords0 && mesh.offsets == offs0
    @test mesh.index == idx0 && mesh.ndofs == nd0 && mesh.connections == conn0
    @test mesh.faces == faces0
    @test u == u0
    # no shared mutable state with the transformed mesh
    @test new.elements !== mesh.elements && new.orders !== mesh.orders
    @test new.offsets !== mesh.offsets && new.index !== mesh.index
    @test v !== u
    # the source stays usable afterwards
    @test two_to_one_balanced(mesh.elements)
end

# ---- (g) BigFloat(256) ------------------------------------------------------------------------
@testset "(g) BigFloat(256) carry-over" begin
    setprecision(BigFloat, 256) do
        blk = Block(31, LinearMap1D(BigFloat(0), BigFloat(1)),
                        LinearMap1D(BigFloat(0), BigFloat(1)))
        root = ElementId(31, SegmentId(0, 0), SegmentId(0, 0))
        N = 5
        mesh = Mesh2D(blk, [root], (N, N))
        mesh = h_refine(mesh, id -> true)
        poly(x, y) = sum((big(0.7)^k) * x^k for k in 0:N-1) *
                     sum((big(-0.4)^k) * y^k for k in 0:N-1)
        u = amr_sample(mesh, blk, poly, BigFloat)
        @test eltype(u) == BigFloat
        new = h_refine(mesh, id -> amr_contains(id, 0.25, 0.25))
        v = transfer_solution(mesh, new, u)
        @test eltype(v) == BigFloat
        @test all(precision(x) == 256 for x in v)
        @test maximum(abs.(v .- amr_sample(new, blk, poly, BigFloat))) < 1e-60
        # Float64 twin of the same transformation agrees to Float64 accuracy
        blkf = amr_blk1()
        meshf = h_refine(Mesh2D(blkf, [amr_root1()], (N, N)), id -> true)
        polyf(x, y) = sum((0.7^k) * x^k for k in 0:N-1) * sum(((-0.4)^k) * y^k for k in 0:N-1)
        uf = amr_sample(meshf, blkf, polyf)
        newf = h_refine(meshf, id -> amr_contains(id, 0.25, 0.25))
        vf = transfer_solution(meshf, newf, uf)
        @test maximum(abs.(Float64.(v) .- vf)) < 1e-12 * maximum(abs.(vf))
    end
end

end # testset
