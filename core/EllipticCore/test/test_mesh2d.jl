# test_mesh2d.jl — grounded gates for src/mesh2d.jl (DAG node: mesh2d-blocks,
# elliptic-ssf GLOBAL_DAG row 9 + TOPOLOGICAL_SCHEDULE W2-B). Written BEFORE the
# implementation (TDD red run captured in evidence/mesh2d-blocks/pkg_test_red.log).
#
# Oracle gates (GLOBAL_DAG row 9 + node contract):
#   * eq:compaction wavezone map r(ξ) = 2 r_u R / (r_u + R + (r_u − R)ξ)  [CRP main.tex:546-549]
#     with EXACT endpoints r(−1) == r_u, r(+1) == R (exact equality, NOT a tolerance), and the
#     1/r-affine grounding identity 1/r(ξ) = (1+ξ)/(2R) + (1−ξ)/(2r_u) on a ξ grid
#   * linear block map with EXACT endpoints x(−1) == a, x(+1) == b
#   * round-trip x(ξ(x)) ≤ 1e-13 relative across BOTH map families (calibrated: worst case
#     5.1e-14 at the paper-realistic contrast (r_u, R) = (25.3, 1e4); linear ~1e-16)
#   * analytic Jacobians: positivity on [−1,1]; closed-form r-map identity
#     dr/dξ = (R − r_u) r² / (2 r_u R); central-FD agreement 1e-10 — Float64 FD for the linear
#     map (truncation-free), BigFloat(256) FD (h = 1e-30) for the compaction map, whose Float64
#     FD noise floor eps·(f/f′)/h + truncation h²((R−r_u)/den)² bottoms out ABOVE 1e-10 at the
#     steep outer end; a Float64-vs-BigFloat analytic bridge at 1e-13 closes the chain
#   * ElementId refinement-tree arithmetic (block id + per-dimension level/index): parent/child
#     inverses, neighbor symmetry, block-boundary neighbors === nothing, dyadic logical
#     intervals EXACT (children bisect the parent interval with shared-face equality)
#   * 2:1-balance predicate over an element set (face adjacency only — "two-to-one balance at
#     element interfaces", main.tex:554): grounded hand cases (normal-dim violation,
#     tangential-dim violation, diff-exactly-1 legality boundary, corner-contact tolerance),
#     then a randomized battery: free random trees vs an INDEPENDENT Rational-arithmetic
#     reference checker; ripple-balanced constructive trees (CRP split-offenders strategy)
#     must all pass; a double-split breaker applied to each must fail
#   * BigFloat(256) pass-through on both maps (endpoints exact, round-trip ≤ 1e-70,
#     Float64-map + BigFloat-ξ promotion)
#
# Test helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using EllipticCore

# ---- tiny deterministic generator (keeps the stdlib Random OUT of the test env; same
# ---- Knuth-MMIX LCG convention as test_refelement.jl, distinct name/seed) ----------------
mutable struct MeshRNG
    s::UInt64
end
next_u!(r::MeshRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f; r.s)
# uniform integer in 1:n (top bits; modulo bias irrelevant at test scale)
rand_1n!(r::MeshRNG, n::Int) = Int(next_u!(r) >> 33) % n + 1

# ---- INDEPENDENT reference for the 2:1 predicate: exact Rational-arithmetic intervals ------
# (the implementation uses integer dyadic shifts; this path uses Rational{Int} — a genuinely
#  different arithmetic route to the same face-adjacency semantics)
seg_rat(s) = (2 * (s.index // (1 << s.level)) - 1,
              2 * ((s.index + 1) // (1 << s.level)) - 1)

function rat_face_adjacent(a, b)
    a.block == b.block || return false
    a1lo, a1hi = seg_rat(a.segments[1]); b1lo, b1hi = seg_rat(b.segments[1])
    a2lo, a2hi = seg_rat(a.segments[2]); b2lo, b2hi = seg_rat(b.segments[2])
    touch1 = a1hi == b1lo || b1hi == a1lo
    touch2 = a2hi == b2lo || b2hi == a2lo
    over1 = max(a1lo, b1lo) < min(a1hi, b1hi)
    over2 = max(a2lo, b2lo) < min(a2hi, b2hi)
    return (touch1 && over2) || (touch2 && over1)
end

function ref_two_to_one(ids)
    for i in 1:length(ids)-1, j in i+1:length(ids)
        a, b = ids[i], ids[j]
        rat_face_adjacent(a, b) || continue
        for d in 1:2
            abs(a.segments[d].level - b.segments[d].level) <= 1 || return false
        end
    end
    return true
end

# ---- refinement-tree generators (leaf sets always tile the block root exactly) -------------
# split leaves[i] once: mode 1 = dim 1 only, 2 = dim 2 only, 3 = both dims (CRP h-split)
function split_at!(leaves, i, mode)
    l = popat!(leaves, i)
    if mode == 1 || mode == 2
        push!(leaves, child_id(l, mode, 0), child_id(l, mode, 1))
    else
        c0, c1 = child_id(l, 1, 0), child_id(l, 1, 1)
        push!(leaves, child_id(c0, 2, 0), child_id(c0, 2, 1),
                      child_id(c1, 2, 0), child_id(c1, 2, 1))
    end
    return leaves
end

function random_tree(rng, nsplits; block=0)
    leaves = [ElementId(block, SegmentId(0, 0), SegmentId(0, 0))]
    for _ in 1:nsplits
        split_at!(leaves, rand_1n!(rng, length(leaves)), rand_1n!(rng, 3))
    end
    return leaves
end

# ripple generator: split l in both dims, first recursively splitting every coarser
# face-adjacent leaf (the CRP "split the offending elements as well" strategy,
# main.tex:554) — the output is 2:1-balanced BY CONSTRUCTION, independent of the
# predicate under test (adjacency via the Rational reference path above)
function split_balanced!(leaves, l)
    while true
        off = findfirst(m -> m != l && rat_face_adjacent(m, l) &&
                             (m.segments[1].level < l.segments[1].level ||
                              m.segments[2].level < l.segments[2].level), leaves)
        off === nothing && break
        split_balanced!(leaves, leaves[off])
    end
    return split_at!(leaves, findfirst(==(l), leaves), 3)
end

# breaker: double-split the globally max-level leaf toward its sibling face — that face is
# block-interior (a sibling exists at level ≥ 1) and every leaf across it has level ≤ L in the
# split dimension, so the new level-(L+2) leaf violates 2:1 across a face. Returns nothing for
# an unrefined tree (never happens in the batteries below — every tree gets ≥ 1 split).
function break_tree(leaves)
    best, bestd, bestlev = nothing, 0, 0
    for m in leaves, d in 1:2
        lev = m.segments[d].level
        if lev > bestlev
            best, bestd, bestlev = m, d, lev
        end
    end
    bestlev == 0 && return nothing
    l, d = best, bestd
    w_side = iseven(l.segments[d].index) ? 1 : 0     # sibling side: even index → hi face
    c_side = child_id(l, d, w_side)
    c_other = child_id(l, d, 1 - w_side)
    g_side = child_id(c_side, d, w_side)
    g_other = child_id(c_side, d, 1 - w_side)
    out = [m for m in leaves if m != l]
    append!(out, (c_other, g_other, g_side))
    return out
end

@testset "mesh2d (mesh2d-blocks)" begin

@testset "linear map: endpoints EXACT, monotone, round-trip 1e-13" begin
    for (a, b) in ((0.0, 1.0), (-1.0, 1.0), (1.7, 123.4), (-5.3, -2.9), (2.1, 1.0e4))
        m = LinearMap1D(a, b)
        @test map_forward(m, -1.0) == a          # exact equality — no tolerance
        @test map_forward(m, 1.0) == b
        @test map_forward(m, -1) == a            # Int logical coordinate promotes
        @test map_forward(m, 1) == b
        @test map_forward(m, 0.0) ≈ (a + b) / 2 rtol = 1e-15
        xs = [map_forward(m, ξ) for ξ in range(-1.0, 1.0, length=101)]
        @test all(diff(xs) .> 0)                 # strictly increasing (orientation-preserving)
        scale = max(abs(a), abs(b))
        for x in range(a, b, length=41)
            @test abs(map_forward(m, map_inverse(m, x)) - x) <= 1e-13 * scale
        end
    end
    @test_throws ArgumentError LinearMap1D(1.0, 1.0)   # degenerate
    @test_throws ArgumentError LinearMap1D(3.0, -1.0)  # orientation-reversing
end

@testset "eq:compaction map: endpoints EXACT, 1/r-affine, round-trip 1e-13" begin
    # r(ξ) = 2 r_u R / (r_u + R + (r_u − R)ξ), main.tex:546-549; includes the paper's own
    # contrast class (R = 1e4 M, r_u = tens of M) and the extreme r_u = 10.
    for (ru, R) in ((25.3, 1.0e4), (100.0, 1.0e4), (10.0, 1.0e4), (2.0, 32.0))
        m = InverseRadialMap1D(ru, R)
        @test map_forward(m, -1.0) == ru         # eq:compaction endpoint — EXACT
        @test map_forward(m, 1.0) == R
        @test map_forward(m, -1) == ru
        @test map_forward(m, 1) == R
        rs = [map_forward(m, ξ) for ξ in range(-1.0, 1.0, length=201)]
        @test all(diff(rs) .> 0)                 # strictly increasing toward the wavezone
        # grounding: 1/r is the affine interpolant between 1/r_u and 1/R (eq:compaction
        # rearranged) — catches any wrong rearrangement of the implementation formula
        for ξ in range(-1.0, 1.0, length=21)
            @test 1 / map_forward(m, ξ) ≈ (1 + ξ) / (2R) + (1 - ξ) / (2ru) rtol = 1e-14
        end
        rsamp = collect(exp.(range(log(ru), log(R), length=61)))
        push!(rsamp, R * (1 - 1e-12))            # hardest point: map is steepest at the outer end
        for r in rsamp
            @test abs(map_forward(m, map_inverse(m, r)) - r) <= 1e-13 * r
        end
    end
    @test_throws ArgumentError InverseRadialMap1D(0.0, 10.0)     # r_u must be positive
    @test_throws ArgumentError InverseRadialMap1D(-2.0, 10.0)
    @test_throws ArgumentError InverseRadialMap1D(10.0, 10.0)    # degenerate
    @test_throws ArgumentError InverseRadialMap1D(100.0, 10.0)   # orientation-reversing
end

@testset "jacobians: analytic, positive, FD-verified 1e-10" begin
    lin = LinearMap1D(1.7, 123.4)
    ru, R = 25.3, 1.0e4
    inv = InverseRadialMap1D(ru, R)
    for ξ in range(-1.0, 1.0, length=21)
        @test map_jacobian(lin, ξ) > 0
        @test map_jacobian(inv, ξ) > 0
        # closed-form identity from eq:compaction: dr/dξ = (R − r_u) r² / (2 r_u R)
        @test map_jacobian(inv, ξ) ≈ (R - ru) * map_forward(inv, ξ)^2 / (2 * ru * R) rtol = 1e-13
    end
    @test map_jacobian(lin, 0.3) == (123.4 - 1.7) / 2            # constant, exact
    @test map_jacobian(inv, -1.0) ≈ (R - ru) * ru / (2R) rtol = 1e-13
    @test map_jacobian(inv, 1.0) ≈ (R - ru) * R / (2ru) rtol = 1e-13
    # Float64 central FD on the linear map (truncation-free; calibrated rel err ~2e-13)
    for ξ in (-0.7, 0.0, 0.4)
        fd = (map_forward(lin, ξ + 1e-4) - map_forward(lin, ξ - 1e-4)) / 2e-4
        @test abs(fd - map_jacobian(lin, ξ)) <= 1e-10 * abs(map_jacobian(lin, ξ))
    end
    # compaction map: the 1e-10 analytic-vs-FD gate runs in BigFloat(256) (h = 1e-30) — the
    # Float64 central-FD floor exceeds 1e-10 at the steep outer end (see file header). The
    # formula under test is the SAME generic code path the Float64 map runs.
    setprecision(BigFloat, 256) do
        h = big"1e-30"
        for m in (LinearMap1D(big"1.7", big"123.4"), InverseRadialMap1D(big"25.3", big"1e4"))
            for ξ64 in range(-1.0, 1.0, length=17)[2:end-1]
                ξ = BigFloat(ξ64)
                fd = (map_forward(m, ξ + h) - map_forward(m, ξ - h)) / (2h)
                jac = map_jacobian(m, ξ)
                @test abs(fd - jac) <= big"1e-10" * abs(jac)
            end
        end
        # bridge: Float64 analytic values agree with BigFloat analytic values at 1e-13,
        # closing the FD chain for the Float64 path (big"25.3" ≠ 25.3 only at ~1e-16 rel)
        mb = InverseRadialMap1D(big"25.3", big"1e4")
        for ξ in (-0.875, -0.25, 0.0, 0.5, 0.9375)
            @test map_jacobian(inv, ξ) ≈ Float64(map_jacobian(mb, BigFloat(ξ))) rtol = 1e-13
            @test map_forward(inv, ξ) ≈ Float64(map_forward(mb, BigFloat(ξ))) rtol = 1e-13
        end
    end
end

@testset "BigFloat(256) pass-through on maps" begin
    setprecision(BigFloat, 256) do
        ru, R = big"25.3", big"1e4"
        m = InverseRadialMap1D(ru, R)
        @test map_forward(m, big"-1.0") == ru    # exact endpoints survive BigFloat
        @test map_forward(m, big"1.0") == R
        r0 = map_forward(m, big"0.125")
        @test r0 isa BigFloat
        @test precision(r0) == 256
        for r in (ru, big"77.7", big"1234.5", R)
            @test abs(map_forward(m, map_inverse(m, r))
                      - r) <= big"1e-70" * r     # genuinely 256-bit accuracy, not Float64's
        end
        lin = LinearMap1D(big"0.0", big"1.0")    # the angular z-block map
        @test map_forward(lin, big"-1.0") == 0
        @test map_forward(lin, big"1.0") == 1
        @test map_jacobian(lin, big"0.0") == big"0.5"
        for x in (big"0.1", big"0.5", big"0.925")
            @test abs(map_forward(lin, map_inverse(lin, x)) - x) <= big"1e-70"
        end
        # Float64-geometry map + BigFloat logical coordinate promotes (precision inherited
        # from the argument, per the DESIGN §4 BigFloat-compatibility discipline)
        mf = InverseRadialMap1D(25.3, 1.0e4)
        rb = map_forward(mf, big"0.25")
        @test rb isa BigFloat
        @test map_jacobian(mf, big"0.25") isa BigFloat
    end
end

@testset "Block: id + per-dimension maps" begin
    blk = Block(3, InverseRadialMap1D(25.3, 1.0e4), LinearMap1D(0.0, 1.0))
    @test blk.id == 3
    @test blockmap(blk, 1) isa InverseRadialMap1D
    @test blockmap(blk, 2) isa LinearMap1D
    @test map_forward(blockmap(blk, 1), -1.0) == 25.3
    @test map_forward(blockmap(blk, 2), 1.0) == 1.0
    # geometry real type must agree across dimensions (RT-consistency is a TYPE gate)
    @test_throws MethodError Block(0, LinearMap1D(0.0, 1.0), LinearMap1D(big"0.0", big"1.0"))
end

@testset "ElementId: refinement-tree arithmetic (exhaustive levels ≤ (4,3))" begin
    segs1 = [SegmentId(l, i) for l in 0:4 for i in 0:(1 << l)-1]
    segs2 = [SegmentId(l, i) for l in 0:3 for i in 0:(1 << l)-1]   # asymmetric: catches dim mixups
    for s1 in segs1, s2 in segs2
        id = ElementId(7, s1, s2)
        @test id.block == 7
        for dim in 1:2
            s = (s1, s2)[dim]
            other = dim == 1 ? 2 : 1
            for which in 0:1
                c = child_id(id, dim, which)
                @test c.segments[dim].level == s.level + 1
                @test c.segments[dim].index == 2 * s.index + which
                @test c.segments[other] == id.segments[other]     # other dim untouched
                @test c.block == id.block
                @test parent_id(c, dim) == id                     # parent ∘ child == identity
            end
            if s.level > 0
                @test child_id(parent_id(id, dim), dim, s.index & 1) == id   # child ∘ parent
            else
                @test_throws ArgumentError parent_id(id, dim)     # root has no parent
            end
            for dir in (-1, 1)
                nb = neighbor_id(id, dim, dir)
                at_edge = dir == -1 ? s.index == 0 : s.index == (1 << s.level) - 1
                if at_edge
                    @test nb === nothing                          # off the block boundary
                else
                    @test nb isa ElementId
                    @test nb.block == id.block
                    @test nb.segments[dim].level == s.level
                    @test nb.segments[other] == id.segments[other]
                    @test neighbor_id(nb, dim, -dir) == id        # neighbor symmetry
                end
            end
        end
    end
    # invariant violations throw
    @test_throws ArgumentError SegmentId(-1, 0)
    @test_throws ArgumentError SegmentId(2, 4)
    @test_throws ArgumentError SegmentId(2, -1)
    @test_throws ArgumentError SegmentId(0, 1)
    @test_throws ArgumentError SegmentId(63, 0)    # level cap keeps Int64 dyadic shifts exact
    root = ElementId(0, SegmentId(0, 0), SegmentId(0, 0))
    @test_throws ArgumentError child_id(root, 1, 2)
    @test_throws ArgumentError child_id(root, 1, -1)
    @test_throws ArgumentError neighbor_id(root, 1, 0)
    @test_throws ArgumentError neighbor_id(root, 1, 2)
    # the block id is part of the address
    @test ElementId(0, SegmentId(1, 0), SegmentId(0, 0)) != ElementId(1, SegmentId(1, 0), SegmentId(0, 0))
end

@testset "logical_interval: dyadic-exact bisection" begin
    root = ElementId(0, SegmentId(0, 0), SegmentId(0, 0))
    @test logical_interval(Float64, root, 1) == (-1.0, 1.0)
    @test logical_interval(Float64, SegmentId(1, 0)) == (-1.0, 0.0)
    @test logical_interval(Float64, SegmentId(1, 1)) == (0.0, 1.0)
    @test logical_interval(Float64, SegmentId(3, 5)) == (0.25, 0.5)
    for l in 0:6, i in 0:(1 << l)-1
        id = ElementId(0, SegmentId(l, i), SegmentId(0, 0))
        lo, hi = logical_interval(Float64, id, 1)
        @test -1.0 <= lo < hi <= 1.0
        mid = (lo + hi) / 2
        # children bisect the parent EXACTLY (shared face, no gap/overlap — dyadic arithmetic)
        @test logical_interval(Float64, child_id(id, 1, 0), 1) == (lo, mid)
        @test logical_interval(Float64, child_id(id, 1, 1), 1) == (mid, hi)
        # BigFloat intervals coincide exactly with the Float64 dyadics
        blo, bhi = logical_interval(BigFloat, id, 1)
        @test blo == lo && bhi == hi
    end
    # adjacent same-level segments share the face point exactly
    for l in 1:6, i in 0:(1 << l)-2
        @test logical_interval(Float64, SegmentId(l, i))[2] ==
              logical_interval(Float64, SegmentId(l, i + 1))[1]
    end
end

@testset "2:1 predicate: grounded hand cases" begin
    S = SegmentId
    E(s1, s2) = ElementId(0, s1, s2)
    @test two_to_one_balanced([E(S(0, 0), S(0, 0))])                    # single root
    quad = vec([E(S(1, i), S(1, j)) for i in 0:1, j in 0:1])
    @test two_to_one_balanced(quad)                                     # uniform level 1
    # refine the [−1,0]² quadrant once more: diff exactly 1 across faces — still balanced
    graded = vcat([E(S(1, 1), S(1, 0)), E(S(1, 0), S(1, 1)), E(S(1, 1), S(1, 1))],
                  vec([E(S(2, i), S(2, j)) for i in 0:1, j in 0:1]))
    @test two_to_one_balanced(graded)
    # normal-dim violation: [−1,0]×[−1,1] (level 1) faces [0,0.25]×[−1,1] (level 3) — diff 2
    @test !two_to_one_balanced([E(S(1, 0), S(0, 0)), E(S(3, 4), S(0, 0))])
    # …but level 2 ([0,0.5]) is the legality boundary: diff exactly 1 passes
    @test two_to_one_balanced([E(S(1, 0), S(0, 0)), E(S(2, 2), S(0, 0))])
    # tangential-dim violation: face contact in dim 1, but dim-2 levels differ by 2
    @test !two_to_one_balanced([E(S(1, 0), S(1, 0)), E(S(1, 1), S(3, 0))])
    # corner-only contact is NOT an interface (main.tex:554): level diff 2 is tolerated
    @test two_to_one_balanced([E(S(1, 0), S(1, 0)), E(S(3, 4), S(3, 4))])
    # a duplicated address never counts as its own neighbor
    @test two_to_one_balanced([E(S(2, 1), S(2, 1)), E(S(2, 1), S(2, 1))])
end

@testset "2:1 predicate: randomized battery vs independent reference" begin
    # seed + consumption order fixed by the calibration run (scratchpad/calibrate_mesh2d.jl):
    # free trees k∈1:6 → 36 balanced / 24 broken; k∈1:25 → 9 / 51; ripple 40/40 + breaker 40/40
    rng = MeshRNG(0x2b1a5eed2b1a5eed)
    nbal, nunbal = 0, 0
    for _ in 1:60
        t = random_tree(rng, rand_1n!(rng, 6))
        got = two_to_one_balanced(t)
        @test got == ref_two_to_one(t)         # exact agreement with the Rational reference
        got ? (nbal += 1) : (nunbal += 1)
    end
    for _ in 1:60
        t = random_tree(rng, rand_1n!(rng, 25))
        got = two_to_one_balanced(t)
        @test got == ref_two_to_one(t)
        got ? (nbal += 1) : (nunbal += 1)
    end
    @test nbal >= 20                            # battery genuinely exercises BOTH classes
    @test nunbal >= 20
    # constructive classes: ripple-balanced ⇒ true, then breaker ⇒ false — every tree
    for _ in 1:40
        leaves = [ElementId(0, SegmentId(0, 0), SegmentId(0, 0))]
        for _ in 1:rand_1n!(rng, 12)
            split_balanced!(leaves, leaves[rand_1n!(rng, length(leaves))])
        end
        @test ref_two_to_one(leaves)            # generator sanity (independent path)
        @test two_to_one_balanced(leaves)
        broken = break_tree(leaves)
        @test broken !== nothing
        @test !two_to_one_balanced(broken)
        @test !ref_two_to_one(broken)
    end
end

@testset "2:1 predicate: multi-block scoping" begin
    rng = MeshRNG(0x0e11_1b7c_0e11_1b7c)
    trees = map((0, 1)) do blk
        leaves = [ElementId(blk, SegmentId(0, 0), SegmentId(0, 0))]
        for _ in 1:6
            split_balanced!(leaves, leaves[rand_1n!(rng, length(leaves))])
        end
        leaves
    end
    @test two_to_one_balanced(vcat(trees...))            # blocks are independent
    @test !two_to_one_balanced(vcat(trees[1], break_tree(trees[2])))
    # same-block violation is NOT masked by other blocks' elements being present
    @test !two_to_one_balanced(vcat(break_tree(trees[1]), trees[2]))
    # cross-block addresses never interact (predicate scope is same-block adjacency;
    # cross-block balance needs block connectivity — the domain/mortar layer's job)
    a = ElementId(0, SegmentId(4, 0), SegmentId(4, 0))
    b = ElementId(1, SegmentId(0, 0), SegmentId(0, 0))
    @test two_to_one_balanced([a, b])
end

end # testset mesh2d
