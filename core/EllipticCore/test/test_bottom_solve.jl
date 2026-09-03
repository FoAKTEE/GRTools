# test_bottom_solve.jl — grounded gates for src/bottom_solve.jl (DAG node:
# multigrid-bottom-solver, elliptic-ssf GLOBAL_DAG row 20 + TOPOLOGICAL_SCHEDULE W8-A,
# CRITICAL PATH). Written BEFORE the implementation (TDD red run captured in
# evidence/multigrid-bottom-solver/tests_red.log).
#
# Target: CRP main.tex:674 — "At the coarsest level, the multigrid solver explicitly
# constructs the matrix A^bottom column-by-column and inverts it directly using the
# Eigen sparse LU decomposition routines (bottom solver). ... The LU decomposition of
# A^bottom is cached and reused in all multigrid V-cycles across AMR iterations, since
# the bottom grid remains the same even when AMR adds finer levels." (fig:precond
# caption, main.tex:632: "The coarsest level is solved directly (bottom solver)".)
# The yaml twin is anc/ScalarSelfForce.yaml:186-194 `UseBottomSolver: True`.
#
# This node fills the SEALED seam of `MultigridV(hier; bottom = (z, b) -> …)` — the
# multigrid-hierarchy node's `dense_bottom_solver` default is unchanged and stays the
# reference implementation these gates measure against. No API change to multigrid.jl.
#
# Oracle gates (constants calibrated PRE-RED in evidence/multigrid-bottom-solver/
# calibrate_bottom{1,2,3}.jl + calibration_round{1,2,3}.log — THROWAWAY prototypes over
# the sealed W5/W6/W7 APIs, NOT the implementation):
#   (1) BOTTOM RESIDUAL ‖A_c z − b‖/‖b‖ < 1e-12 on the coarsest level of a REAL
#       hierarchy (two affine blocks joined by a BlockConnection, complex coefficients
#       β = 1 + i/2): calibrated 6.28e-16 at ComplexF64 and 5.58e-16 at Float64, and
#       6.9e-16 on the CRP-shaped production coarsest grid (10 blocks × one 4×4 element
#       = 160 dofs, the fig:matrix caption's "Each element has 4×4 grid points"). At
#       Complex{BigFloat}(256) the generic dense fallback is calibrated 4.94e-77
#       (gated < 1e-60). `b` is not modified, and the aliased call `bottom!(z, z)` is
#       BITWISE the un-aliased one (calibrated 0.0).
#   (2) CACHING — the LU is built EXACTLY ONCE across a mock 8-step AMR run that refines
#       the FINEST level only (a graded puncture cascade in the corner of block 1; the
#       coarsest grid is the irreducible block structure and never changes): calibrated
#       `builds == 1`, `hits == 7`, with the coarsest element set / orders / ndofs
#       invariant across all 8 rebuilds. Across V-cycles the count never moves.
#       FALSIFIERS — each of these MUST rebuild: an `epoch` bump; a coefficient change
#       (β, or 𝔄¹¹ perturbed by 1e-12) at the SAME epoch; a penalty-`C` change; one
#       boundary face switched Dirichlet → Neumann; an order change (`p_refine`) or any
#       other change of the coarsest mesh (`maxlevels` truncation); and a Float64
#       operator on the SAME mesh as a ComplexF64 one (no cross-type collision). The
#       sharp form of the coefficient falsifier is numerical: after the change the
#       solver must return the NEW matrix's solution, residual < 1e-12.
#   (3) EQUIVALENCE — the sparse route equals the sealed `dense_bottom_solver` to
#       ≤ 1e-13 (calibrated 1.36e-15) on the same coarsest system; the V-cycle
#       contraction is identical to 1e-12; and GMRES ∘ MultigridV with the sparse bottom
#       reproduces the SEALED iteration counts (7 at both 3 and 4 levels,
#       ESSF::multigrid-hierarchy) on the manufactured discrete solution b := A·u_ex at
#       reltol 1e-10, converging to it with relative error ≤ 1e-9 (calibrated 4.39e-10
#       at 3 levels, 5.82e-10 at 4) and agreeing with the dense-bottom route to 1e-12.
#   (4) PRODUCTION-SHAPED COST is MEASURED, not gated (evidence/measure_landing_w8a.jl):
#       at the CRP coarsest shape (160 dofs, 8.88 % dense) sparse and dense are within a
#       factor ≈ 2 of each other in BOTH directions and each is ≈ 2 % of the
#       column-probe assembly that precedes them, so the sparse route buys nothing
#       measurable AT THAT SIZE; it wins by 10.7× (factor) / 6.8× (solve) at 640 dofs
#       and by 166× / 41× at 2560 dofs, which is why the paper's general-operator solver
#       uses a sparse LU. Correctness, not speed, is what gate (3) pins.
#   (5) APPLY PATH — after warm-up an apply allocates a CONSTANT number of bytes
#       (calibrated 0 with the three-argument `ldiv!(z, F, b)` route; the two-argument
#       `ldiv!(F, z)` copies internally and costs 2632 B) and NEVER rebuilds: the cache
#       counters are unchanged after 50 applies, and exactly one bottom solve runs per
#       V-cycle (counter).
#   (6) SCALAR-TYPE ROUTING — UMFPACK takes exactly `Float64` and `ComplexF64`; anything
#       else falls back to the generic DENSE LU. Both reasons are pinned here as
#       SparseArrays facts: `lu` of a `BigFloat` sparse matrix THROWS (ArgumentError,
#       "matrix type … not supported"), and `lu` of a `Float32`/`ComplexF32` sparse
#       matrix SILENTLY PROMOTES to `UmfpackLU{Float64}` / `UmfpackLU{ComplexF64}` —
#       so the dispatch is on the exact pair, never on `AbstractFloat`.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting siblings).

using Test
using LinearAlgebra
using SparseArrays
using REPL          # Julia >=1.11: @doc is a no-op unless REPL is loaded (runtests.jl convention)
using EllipticCore

# ---- calibrated gate constants (evidence/multigrid-bottom-solver/calibration_round*.log) ----
const BS_RES_TOL = 1e-12        # calibrated 6.28e-16 (ComplexF64) / 5.58e-16 (Float64)
const BS_RES_TOL_BIG = 1e-60    # calibrated 4.94e-77 at Complex{BigFloat}(256)
const BS_AGREE_TOL = 1e-13      # sparse vs dense_bottom_solver: calibrated 1.36e-15
const BS_RHO_TOL = 1e-12        # V-cycle contraction, sparse vs dense bottom
const BS_GMRES_ITERS = 7        # SEALED count (ESSF::multigrid-hierarchy) at 3 and 4 levels
const BS_SOLERR_TOL = 1e-9      # calibrated 4.39e-10 / 5.82e-10 at reltol 1e-10
const BS_XAGREE_TOL = 1e-12     # sparse-bottom vs dense-bottom GMRES solution
const BS_ALLOC_MAX = 128        # calibrated 0 B (three-argument ldiv! route)
const BS_AMR_STEPS = 8          # mock AMR run length: builds == 1, hits == 7

# ---- model problems (coefficient callables only — zero physics) --------------------------
const BS_C = 1.5                                   # CRP main.tex:595 penalty parameter
bs_A11(x, y) = 1.0
bs_A22(x, y) = 1.0
bs_zero2(x, y) = 0.0
const BS_BETA = 1.0 + 0.5im                        # definite complex model

function bs_zerobcs(mesh, T)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(T, nf))
    end
    return bcs
end

# operator factory: the SAME coefficient closures rediscretized on whatever mesh it is
# handed — what MultigridHierarchy needs, and what an AMR step re-runs
bs_factory(T, β; C = BS_C, A11 = bs_A11) =
    mesh -> MeshOperator(T, mesh, (A11, bs_A22), (x, y) -> β, (bs_zero2, bs_zero2);
                         C = C, bcs = bs_zerobcs(mesh, T))

# two affine unit blocks joined along dimension 1 (the sealed multigrid model domain)
function bs_twoblocks(RT)
    b1 = Block(1, LinearMap1D(RT(0), RT(1)), LinearMap1D(RT(0), RT(1)))
    b2 = Block(2, LinearMap1D(RT(1), RT(2)), LinearMap1D(RT(0), RT(1)))
    return [b1, b2], [BlockConnection(1, 1, 2)]
end

function bs_uniform_mesh(RT, k, N)
    blks, conns = bs_twoblocks(RT)
    els = ElementId[]
    for b in (1, 2), iy in 0:(2^k-1), ix in 0:(2^k-1)
        push!(els, ElementId(b, SegmentId(k, ix), SegmentId(k, iy)))
    end
    return Mesh2D(blks, els, (N, N); connections = conns)
end

# the CRP PRODUCTION coarsest grid, modelled: the irreducible block structure of a
# 10-block domain, one 4×4 element per block (fig:matrix caption, main.tex:632)
function bs_crp_coarse_mesh(RT, nb = 10)
    blks = [Block(b, LinearMap1D(RT(b - 1), RT(b)), LinearMap1D(RT(0), RT(1)))
            for b in 1:nb]
    conns = [BlockConnection(b, 1, b + 1) for b in 1:(nb-1)]
    els = [ElementId(b, SegmentId(0, 0), SegmentId(0, 0)) for b in 1:nb]
    return Mesh2D(blks, els, (4, 4); connections = conns)
end

# the graded "refine the finest level only" flag set: the deepest elements in the
# lower-left corner of block 1 (the CRP puncture cascade). Uniform refinement of every
# finest element explodes the mesh (32 → 131072 dofs in six steps, measured).
function bs_corner_flags(mesh)
    maxlev = maximum(e -> e.segments[1].level + e.segments[2].level, mesh.elements)
    return Set(el for el in mesh.elements
               if el.block == 1 &&
                  el.segments[1].level + el.segments[2].level == maxlev &&
                  el.segments[1].index == 0 && el.segments[2].index == 0)
end

# deterministic RHS / start vectors (Knuth-MMIX LCG; distinct name and seed per test file)
mutable struct BSRNG
    s::UInt64
end
bs_unit!(r::BSRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
function bs_vec!(r::BSRNG, T, n)
    RT = real(T)
    return T[T <: Complex ? T(RT(bs_unit!(r)), RT(bs_unit!(r))) : T(RT(bs_unit!(r)))
             for _ in 1:n]
end

bs_relres(A, z, b) = norm(A * z - b) / norm(b)

# stationary V-cycle iteration on A·u = 0 from a fixed start: the error IS the iterate
function bs_contraction(mg, T, n; niter = 6, seed = 0xbeef)
    x = bs_vec!(BSRNG(seed), T, n)
    b = zeros(T, n)
    ns = Float64[norm(x)]
    for _ in 1:niter
        vcycle!(x, mg, b)
        push!(ns, norm(x))
    end
    return [ns[i+1] / ns[i] for i in 1:length(ns)-1]
end

# the coarsest level of a fresh hierarchy over `mesh`, plus its dense oracle matrix
function bs_coarse(T, mesh; β = BS_BETA, C = BS_C, A11 = bs_A11)
    fac = bs_factory(T, β; C = C, A11 = A11)
    hier = MultigridHierarchy(fac(mesh), fac)
    return hier, hier.levels[end].op
end

# a signature of the coarsest grid — what must stay invariant across an AMR run
bs_signature(op) = (copy(op.mesh.elements), copy(op.mesh.orders), op.mesh.ndofs)


@testset verbose = true "bottom solver (multigrid-bottom-solver)" begin

    @testset "exports, docstrings, and the untouched seam" begin
        for name in (:BottomSolverCache, :sparse_bottom_solver)
            @test name in names(EllipticCore)
            doc = string(Base.Docs.doc(Base.Docs.Binding(EllipticCore, name)))
            @test !occursin("No documentation found", doc)
            @test length(doc) > 200
        end
        # the sealed seam and its default are UNCHANGED by this node
        @test :dense_bottom_solver in names(EllipticCore)
        @test :MultigridV in names(EllipticCore)
        cache = BottomSolverCache()
        @test cache.builds == 0 && cache.hits == 0
        @test isempty(cache.entries)
    end

    @testset "gate 1: bottom residual on a real hierarchy (ComplexF64 / Float64)" begin
        for (T, β) in ((ComplexF64, BS_BETA), (Float64, 1.0))
            hier, cop = bs_coarse(T, bs_uniform_mesh(Float64, 3, 4); β = β)
            @test length(hier.levels) >= 3
            n = cop.mesh.ndofs
            A = build_matrix(cop)
            b = bs_vec!(BSRNG(0x11), T, n)
            b0 = copy(b)
            z = zeros(T, n)
            bottom! = sparse_bottom_solver(cop)
            @test bottom!(z, b) === z
            @test bs_relres(A, z, b) < BS_RES_TOL
            @test b == b0                       # the right-hand side is not consumed
            # aliased call: z and b the same array
            za = copy(b)
            bottom!(za, za)
            @test za == z                       # BITWISE identical to the un-aliased route
        end
    end

    @testset "gate 1: the CRP production coarsest grid (10 blocks × 16 dofs)" begin
        T = ComplexF64
        mesh = bs_crp_coarse_mesh(Float64)
        @test mesh.ndofs == 160
        op = bs_factory(T, BS_BETA)(mesh)
        A = build_matrix(op)
        As = sparse(A)
        @test nnz(As) == 2272                    # 8.88 % dense (calibrated)
        n = mesh.ndofs
        b = bs_vec!(BSRNG(0x22), T, n)
        z = zeros(T, n)
        sparse_bottom_solver(op)(z, b)
        @test bs_relres(A, z, b) < BS_RES_TOL
        # …and it IS the coarsest level of its own hierarchy (one element per block)
        hier, cop = bs_coarse(T, mesh)
        @test length(hier.levels) == 1
        @test cop.mesh.ndofs == 160
    end

    @testset "gate 1: Complex{BigFloat}(256) generic dense fallback" begin
        setprecision(BigFloat, 256) do
            T = Complex{BigFloat}
            β = Complex{BigFloat}(1.0, 0.5)
            hier, cop = bs_coarse(T, bs_uniform_mesh(BigFloat, 2, 4); β = β)
            n = cop.mesh.ndofs
            A = build_matrix(cop)
            b = bs_vec!(BSRNG(0x33), T, n)
            z = zeros(T, n)
            sparse_bottom_solver(cop)(z, b)
            r = bs_relres(A, z, b)
            @test r < BS_RES_TOL_BIG
            @test precision(BigFloat) == 256     # the fallback keeps the working precision
            @test eltype(z) === T                # …and the scalar type (no Float64 demotion)
        end
    end

    @testset "gate 6: why the fallback exists (SparseArrays routing facts)" begin
        # UMFPACK cannot take BigFloat at all …
        @test_throws ArgumentError lu(sparse(BigFloat[4 1 0; 1 4 1; 0 1 4]))
        @test_throws ArgumentError lu(sparse(Complex{BigFloat}[4 1; 1 4]))
        # … and SILENTLY PROMOTES the small floats, which would change the scalar type
        # under the caller's feet — hence dispatch on the exact {Float64, ComplexF64} pair.
        @test lu(sparse(Float32[4 1; 1 4])) isa SparseArrays.UMFPACK.UmfpackLU{Float64}
        @test lu(sparse(ComplexF32[4 1; 1 4])) isa
              SparseArrays.UMFPACK.UmfpackLU{ComplexF64}
        # the generic dense LU takes all four
        for T in (BigFloat, Complex{BigFloat}, Float32, ComplexF32)
            @test lu(T[4 1; 1 4]) isa LU{T}
        end
    end

    @testset "gate 2: ONE LU across a mock 8-step AMR run (finest level only)" begin
        T = ComplexF64
        fac = bs_factory(T, BS_BETA)
        mesh = bs_uniform_mesh(Float64, 2, 4)
        cache = BottomSolverCache()
        subcache = SubdomainCache()
        sigs = []
        nlev = Int[]
        for it in 1:BS_AMR_STEPS
            op = fac(mesh)
            hier = MultigridHierarchy(op, fac; caches = subcache)
            cop = hier.levels[end].op
            push!(sigs, bs_signature(cop))
            push!(nlev, length(hier.levels))
            mg = MultigridV(hier; bottom = sparse_bottom_solver(cop; cache = cache))
            x = zeros(T, op.mesh.ndofs)
            vcycle!(x, mg, bs_vec!(BSRNG(0x44 + it), T, op.mesh.ndofs))
            @test all(isfinite, x)
            mesh = h_refine(mesh, bs_corner_flags(mesh))
        end
        @test cache.builds == 1                  # main.tex:674 — cached across AMR iterations
        @test cache.hits == BS_AMR_STEPS - 1
        @test length(cache.entries) == 1
        @test all(==(sigs[1]), sigs)             # the bottom grid really did not change …
        @test nlev[end] > nlev[1]                # … while AMR added finer levels
    end

    @testset "gate 2: caching across V-cycles (no rebuild on the apply path)" begin
        T = ComplexF64
        hier, cop = bs_coarse(T, bs_uniform_mesh(Float64, 3, 4))
        cache = BottomSolverCache()
        mg = MultigridV(hier; bottom = sparse_bottom_solver(cop; cache = cache))
        @test cache.builds == 1
        n = hier.levels[1].mesh.ndofs
        x = zeros(T, n)
        b = bs_vec!(BSRNG(0x55), T, n)
        for _ in 1:10
            vcycle!(x, mg, b)
        end
        @test cache.builds == 1
        @test cache.hits == 0
    end

    @testset "gate 2 falsifiers: every change of the coarse matrix rebuilds" begin
        T = ComplexF64
        mesh = bs_uniform_mesh(Float64, 2, 4)
        fac = bs_factory(T, BS_BETA)
        hier, cop = bs_coarse(T, mesh)
        cmesh = cop.mesh
        cache = BottomSolverCache()
        sparse_bottom_solver(cop; cache = cache)
        @test (cache.builds, cache.hits) == (1, 0)
        # the SAME operator object, and an identically rebuilt one: hits
        sparse_bottom_solver(cop; cache = cache)
        @test (cache.builds, cache.hits) == (1, 1)
        sparse_bottom_solver(fac(cmesh); cache = cache)
        @test (cache.builds, cache.hits) == (1, 2)
        # (a) explicit epoch bump
        sparse_bottom_solver(cop; cache = cache, epoch = 1)
        @test (cache.builds, cache.hits) == (2, 2)
        # (b) coefficient change at the SAME epoch — and the SHARP numerical falsifier
        op2 = bs_factory(T, -3.0 + 0.1im)(cmesh)
        bot2 = sparse_bottom_solver(op2; cache = cache)
        @test (cache.builds, cache.hits) == (3, 2)
        A2 = build_matrix(op2)
        b = bs_vec!(BSRNG(0x66), T, cmesh.ndofs)
        z = zeros(T, cmesh.ndofs)
        bot2(z, b)
        @test bs_relres(A2, z, b) < BS_RES_TOL           # NOT the stale β = 1 + i/2 solve
        # (c) a 1e-12 perturbation of a coefficient is still a different matrix
        sparse_bottom_solver(bs_factory(T, BS_BETA; A11 = (x, y) -> 1.0 + 1e-12)(cmesh);
                             cache = cache)
        @test (cache.builds, cache.hits) == (4, 2)
        # (d) penalty change
        sparse_bottom_solver(bs_factory(T, BS_BETA; C = 2.0)(cmesh); cache = cache)
        @test (cache.builds, cache.hits) == (5, 2)
        # (e) one boundary face Dirichlet -> Neumann
        bcs2 = bs_zerobcs(cmesh, T)
        k1 = first(sort(collect(keys(bcs2)); by = string))
        nf = cmesh.orders[cmesh.index[k1[1]]][3-k1[2]]
        bcs2[k1] = NeumannGhost(zeros(T, nf))
        sparse_bottom_solver(MeshOperator(T, cmesh, (bs_A11, bs_A22), (x, y) -> BS_BETA,
                                          (bs_zero2, bs_zero2); C = BS_C, bcs = bcs2);
                             cache = cache)
        @test (cache.builds, cache.hits) == (6, 2)
        # (f) an order change (p_refine) moves the coarsest grid
        _, pcop = bs_coarse(T, p_refine(mesh))
        sparse_bottom_solver(pcop; cache = cache)
        @test (cache.builds, cache.hits) == (7, 2)
        # (g) a truncated hierarchy has a DIFFERENT (finer) coarsest mesh
        h2 = MultigridHierarchy(fac(mesh), fac; maxlevels = 2)
        @test h2.levels[end].mesh.ndofs != cmesh.ndofs
        sparse_bottom_solver(h2.levels[end].op; cache = cache)
        @test (cache.builds, cache.hits) == (8, 2)
        # (h) a Float64 operator on the SAME mesh must not collide with the complex one
        sparse_bottom_solver(bs_factory(Float64, 1.0)(cmesh); cache = cache)
        @test (cache.builds, cache.hits) == (9, 2)
        @test length(cache.entries) == 9
        # a cache is caller-owned: a fresh one shares nothing
        c2 = BottomSolverCache()
        sparse_bottom_solver(cop; cache = c2)
        @test (c2.builds, c2.hits) == (1, 0)
        @test (cache.builds, cache.hits) == (9, 2)
    end

    @testset "gate 3: equivalence with the sealed dense_bottom_solver" begin
        T = ComplexF64
        hier, cop = bs_coarse(T, bs_uniform_mesh(Float64, 3, 4))
        n = cop.mesh.ndofs
        A = build_matrix(cop)
        b = bs_vec!(BSRNG(0x77), T, n)
        zs = zeros(T, n)
        zd = zeros(T, n)
        sparse_bottom_solver(cop)(zs, b)
        dense_bottom_solver(cop)(zd, b)
        @test norm(zs - zd) / norm(zd) < BS_AGREE_TOL
        @test bs_relres(A, zs, b) < BS_RES_TOL
        # …and the V-cycle built on either bottom contracts identically
        nf = hier.levels[1].mesh.ndofs
        ρs = bs_contraction(MultigridV(hier; bottom = sparse_bottom_solver(cop)), T, nf)
        ρd = bs_contraction(MultigridV(hier; bottom = dense_bottom_solver(cop)), T, nf)
        @test maximum(abs, ρs - ρd) < BS_RHO_TOL
        @test ρs[end] < 0.25
    end

    @testset "gate 3: GMRES ∘ MultigridV reproduces the sealed iteration counts" begin
        T = ComplexF64
        for k in (2, 3)
            mesh = bs_uniform_mesh(Float64, k, 4)
            hier, cop = bs_coarse(T, mesh)
            n = mesh.ndofs
            uex = bs_vec!(BSRNG(0xabc), T, n)
            b = similar(uex)
            apply!(b, hier.levels[1].op, uex)
            rs = gmres(hier.levels[1].op, b;
                       Minv! = MultigridV(hier; bottom = sparse_bottom_solver(cop)),
                       reltol = 1e-10, maxiter = 100)
            rd = gmres(hier.levels[1].op, b;
                       Minv! = MultigridV(hier; bottom = dense_bottom_solver(cop)),
                       reltol = 1e-10, maxiter = 100)
            @test rs.converged
            @test rs.iterations == BS_GMRES_ITERS          # SEALED count, level-independent
            @test rs.iterations == rd.iterations
            @test norm(rs.x - uex) / norm(uex) < BS_SOLERR_TOL
            @test norm(rs.x - rd.x) / norm(rd.x) < BS_XAGREE_TOL
        end
    end

    @testset "gate 5: the apply path allocates a constant amount and never rebuilds" begin
        for (T, β) in ((ComplexF64, BS_BETA), (Float64, 1.0))
            _, cop = bs_coarse(T, bs_uniform_mesh(Float64, 3, 4); β = β)
            n = cop.mesh.ndofs
            cache = BottomSolverCache()
            bottom! = sparse_bottom_solver(cop; cache = cache)
            b = bs_vec!(BSRNG(0x88), T, n)
            z = zeros(T, n)
            for _ in 1:3
                bottom!(z, b)                              # warm-up
            end
            allocs = [(@allocated bottom!(z, b)) for _ in 1:3]
            @test allocs[1] == allocs[2] == allocs[3]      # CONSTANT after warm-up
            @test allocs[1] <= BS_ALLOC_MAX
            for _ in 1:50
                bottom!(z, b)
            end
            @test (cache.builds, cache.hits) == (1, 0)     # the apply path never touches it
        end
    end

    @testset "gate 5: exactly one bottom solve per V-cycle" begin
        T = ComplexF64
        hier, cop = bs_coarse(T, bs_uniform_mesh(Float64, 3, 4))
        calls = Ref(0)
        inner = sparse_bottom_solver(cop)
        counted = (z, b) -> (calls[] += 1; inner(z, b))
        mg = MultigridV(hier; bottom = counted)
        n = hier.levels[1].mesh.ndofs
        x = zeros(T, n)
        b = bs_vec!(BSRNG(0x99), T, n)
        vcycle!(x, mg, b)
        @test calls[] == 1
        vcycle!(x, mg, b)
        vcycle!(x, mg, b)
        @test calls[] == 3
    end

    @testset "errors: length mismatches are refused" begin
        T = ComplexF64
        _, cop = bs_coarse(T, bs_uniform_mesh(Float64, 2, 4))
        n = cop.mesh.ndofs
        bottom! = sparse_bottom_solver(cop)
        @test_throws ArgumentError bottom!(zeros(T, n + 1), zeros(T, n))
        @test_throws ArgumentError bottom!(zeros(T, n), zeros(T, n - 1))
        @test_throws ArgumentError bottom!(zeros(T, 1), zeros(T, 1))
    end
end
