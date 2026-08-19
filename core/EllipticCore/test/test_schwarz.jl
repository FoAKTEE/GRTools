# test_schwarz.jl — grounded gates for src/schwarz.jl (DAG node: schwarz-smoother,
# elliptic-ssf GLOBAL_DAG row 17 + TOPOLOGICAL_SCHEDULE W6-B, CRITICAL PATH). Written
# BEFORE the implementation (TDD red run captured in
# evidence/schwarz-smoother/red_targeted.log).
#
# Target: the element-centered additive Schwarz smoother of CRP main.tex:676 — one
# subdomain per DG element covering the element plus a configurable `overlap` of points
# into each face neighbor (spectre element-centered subdomain geometry, Vu_2022);
# subdomain matrices built PATCH-LOCALLY by RESTRICTED column probes of the global
# matrix-free LINEARIZED apply (correct off-diagonal couplings within the patch, zero
# outside — per-element matrices, NOT translation-class caching: curved maps, variable
# coefficients, hooks and mortars are all in the patch by construction); cached complex
# LU factorizations keyed (element id, subdomain shape, coefficient epoch, shift) with
# builds/hits counters; weighted overlap combination with EXACT partition-of-unity
# weights; one smoothing pass = solve all subdomains + weighted combine, threaded with
# a deterministic (bitwise == serial) gather-fold reduction; diagonal complex shift
# −i·s threaded through for the schwarz-complex-shift node.
#
# Oracle gates (worker node contract; constants calibrated PRE-RED in scratchpad
# calibrate_schwarz.jl + calibrate_schwarz_scan.jl from the sealed W6-A layer plus a
# THROWAWAY dense-restriction prototype — NOT the implementation):
#   1. subdomain matrix == restriction of build_matrix(A) to the patch, rel ≤ 1e-13,
#      on the mixed h/p mesh (one 2:1 mortar interface patch, p-jump patches, variable
#      𝔄, complex β, nonzero γ) AND on a hook-modified patch (vtu-style LinearizedHook
#      terms must be IN the subdomain matrix — spectre SubdomainOperator.hpp:628-680);
#      calibrated 0.0 (the patch probes reproduce the global probes bitwise).
#   2. PoU weights sum to 1 EXACTLY on every dof (incl. overlaps and corners): weights
#      are quantized to the 2^-30 grid with the largest-weight anchor taking the exact
#      residual, so every fold order is exact; calibrated 0 non-exact dofs.
#   3. 2 passes damp high-frequency error > 10x on the Poisson model (4x4 conforming,
#      N=5, ov=2), residual ratio ‖A·u0‖/‖A·u2‖: calibrated 77.2 (seeds 69.7-78.3) for
#      the random top-Legendre-mode vector and 90.0 for white noise. The white-noise
#      gate is the WEIGHTING falsifier: degenerate RAS weights (the SpecCore
#      schwarz.jl:62 _wintr sign bug — band weights identically zero after
#      normalization) measure 7.6 < 10 on the same model.
#   4. cache-reuse counter: construction builds one LU per subdomain, smoothing passes
#      build ZERO new LUs, an identically-shaped second smoother on the same cache is
#      all hits, an epoch bump or shift change rebuilds.
#   5. THREADED sweep bitwise == serial (JULIA_NUM_THREADS=4).
#   6. complex-coefficient patch (Helmholtz β = −25 − 0.4i) solves match the dense
#      patch backslash ≤ 1e-12 (calibrated 0.0) and one smoothing pass from u=0 equals
#      the test-side weighted-combine reference ≤ 1e-12 (calibrated 0.0).
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting siblings).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test
# ---- files, distinct name/seed) --------------------------------------------------------
mutable struct SWRNG
    s::UInt64
end
sw_unit!(r::SWRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
sw_cvec!(r::SWRNG, n) = [complex(sw_unit!(r), sw_unit!(r)) for _ in 1:n]

const SW_C = 1.5                                     # CRP main.tex:595 penalty parameter

sw_one2(x, y) = 1.0
sw_zero2(x, y) = 0.0
sw_A11(x, y) = 1 + 0.3x + 0.1y^2                     # variable diagonal flux matrix
sw_A22(x, y) = 1 + 0.2y
const SW_βc = 0.3 + 0.2im
const SW_βH = -25.0 - 0.4im                          # complex Helmholtz β = −k² − iε
sw_γ1(x, y) = 0.05
sw_γ2(x, y) = -0.03

# ---- the mixed h/p mesh (one 2:1 interface W | E11,E12 and two p-jumps at orders
# ---- 4 vs 5) on the affine unit block — the wave-6-A headline mesh shape ---------------
sw_blk() = Block(71, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
const SWW   = ElementId(71, SegmentId(1, 0), SegmentId(0, 0))
const SWE11 = ElementId(71, SegmentId(2, 2), SegmentId(1, 0))
const SWE12 = ElementId(71, SegmentId(2, 2), SegmentId(1, 1))
const SWE21 = ElementId(71, SegmentId(2, 3), SegmentId(1, 0))
const SWE22 = ElementId(71, SegmentId(2, 3), SegmentId(1, 1))
sw_els() = [SWW, SWE11, SWE12, SWE21, SWE22]
sw_ords() = Dict(SWW => (4, 4), SWE11 => (4, 4), SWE12 => (4, 4),
                 SWE21 => (5, 5), SWE22 => (5, 5))
sw_mixed_mesh() = Mesh2D(sw_blk(), sw_els(), sw_ords())

# conforming 4x4 Poisson model mesh (the damping-gate model)
sw_pblk() = Block(72, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
sw_pmesh(N) = Mesh2D(sw_pblk(),
                     vec([ElementId(72, SegmentId(2, ix), SegmentId(2, iy))
                          for ix in 0:3, iy in 0:3]), (N, N))

# zero-Dirichlet ghost bcs on every boundary face (through the public topology)
function sw_zerobcs(mesh)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(ComplexF64, nf))
    end
    return bcs
end

sw_mixed_op(β) = MeshOperator(ComplexF64, sw_mixed_mesh(), (sw_A11, sw_A22),
                              (x, y) -> β, (sw_γ1, sw_γ2); C = SW_C,
                              bcs = sw_zerobcs(sw_mixed_mesh()))

# ---- random high-frequency vector: random complex coefficients on the top per-element
# ---- Legendre modes (max(ix,iy) > N - keep), synthesized at the Gauss nodes ------------
function sw_legvals(n, ξ)
    V = zeros(length(ξ), n)
    V[:, 1] .= 1.0
    n > 1 && (V[:, 2] .= ξ)
    for k in 2:n-1
        V[:, k+1] .= ((2k - 1) .* ξ .* V[:, k] .- (k - 1) .* V[:, k-1]) ./ k
    end
    return V
end
function sw_hfvec(mesh, rng; keep = 2)
    u = zeros(ComplexF64, mesh.ndofs)
    for e in eachindex(mesh.elements)
        Nx, Ny = mesh.orders[e]
        Vx = sw_legvals(Nx, EllipticCore.refelement(Float64, Nx).nodes)
        Vy = sw_legvals(Ny, EllipticCore.refelement(Float64, Ny).nodes)
        Ue = zeros(ComplexF64, Nx, Ny)
        for iy in 1:Ny, ix in 1:Nx
            c = complex(sw_unit!(rng), sw_unit!(rng))          # draw ALWAYS (determinism)
            if max(ix, iy) > max(Nx, Ny) - keep
                Ue .+= c .* (Vx[:, ix] * transpose(Vy[:, iy]))
            end
        end
        u[mesh.offsets[e]+1:mesh.offsets[e]+Nx*Ny] .= vec(Ue)
    end
    return u
end

sw_resnorm(op, u) = norm(apply!(zeros(ComplexF64, op.mesh.ndofs), op, u))

@testset "schwarz (schwarz-smoother)" begin

@testset "element-centered subdomain geometry + EXACT partition of unity" begin
    op = sw_mixed_op(SW_βc)
    mesh = op.mesh
    sm = SchwarzSmoother(op; overlap = 2)
    @test length(sm.subdomains) == 5                 # one subdomain per element
    # calibrated patch sizes: central prod(orders) + ov × transverse order per face
    # neighbor (full transverse extent, incl. BOTH halves across the 2:1 face of W)
    @test [length(s.gdofs) for s in sm.subdomains] == [32, 42, 42, 43, 43]
    sm1 = SchwarzSmoother(op; overlap = 1)
    @test [length(s.gdofs) for s in sm1.subdomains] == [24, 29, 29, 34, 34]
    # central element block comes first, in the global layout order
    @test sm.subdomains[1].gdofs[1:16] == collect(1:16)          # W owns dofs 1:16
    @test sm.subdomains[4].gdofs[1:25] == collect(49:73)         # E21 owns dofs 49:73
    # the W patch reaches exactly ov=2 low-x planes into E11 (offset 16, i ∈ 1:2)
    band = [16 + tensor_lin(i, j, 4) for j in 1:4 for i in 1:2]
    @test issubset(band, sm.subdomains[1].gdofs)
    @test !any(g -> g in sm.subdomains[1].gdofs,
               [16 + tensor_lin(3, j, 4) for j in 1:4])          # and no deeper
    # no duplicate dofs in any patch
    @test all(s -> length(unique(s.gdofs)) == length(s.gdofs), sm.subdomains)

    # NODE GATE: PoU weights sum to 1 EXACTLY on every dof — including overlap bands
    # and corners (quantized 2^-30 weights + exact-residual anchor; calibrated 0
    # non-exact dofs). The fold below is the same subdomain-order accumulation the
    # sweep uses; quantization makes every other fold order exact too.
    wsum = zeros(mesh.ndofs)
    for s in eachindex(sm.subdomains), k in eachindex(sm.subdomains[s].gdofs)
        wsum[sm.subdomains[s].gdofs[k]] += sm.subdomains[s].weights[k]
    end
    @test all(==(1.0), wsum)
    @test all(s -> all(w -> 0 ≤ w ≤ 1, s.weights), sm.subdomains)
    # overlap dofs genuinely shared: every band dof carries weight < 1 somewhere
    @test any(s -> any(<(1), s.weights), sm.subdomains)

    @test_throws ArgumentError SchwarzSmoother(op; overlap = -1)
end

@testset "patch matrix == restriction of build_matrix (incl. 2:1 mortar patch)" begin
    op = sw_mixed_op(SW_βc)
    A = build_matrix(op)
    sc = maximum(abs, A)
    sm = SchwarzSmoother(op; overlap = 2)
    # NODE GATE rel ≤ 1e-13 (calibrated 0.0 — restricted probes == global probes
    # bitwise). Subdomain 1 is the 2:1 mortar interface patch (W + both half-face
    # neighbors); subdomains 2-5 cover parent-neighbor and p-jump (4|5) patches.
    for s in 1:5
        g = sm.subdomains[s].gdofs
        As = subdomain_matrix(sm, s)
        @test size(As) == (length(g), length(g))
        @test maximum(abs, As .- A[g, g]) ≤ 1e-13 * sc
    end
end

@testset "hook-modified patch: LinearizedHook terms are IN the subdomain matrix" begin
    # vtu-style flux jump −iω⟨⟨ψ⟩⟩ on both sides of the 2:1 interface (the spectre
    # SubdomainOperator includes modify_boundary_data in the subdomain matrices)
    ω = 0.9
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    hooks = Dict((SWW, 1, 2) => LinearizedHook(fω),
                 (SWE11, 1, 1) => LinearizedHook(fω),
                 (SWE12, 1, 1) => LinearizedHook(fω))
    mesh = sw_mixed_mesh()
    oph = MeshOperator(ComplexF64, mesh, (sw_A11, sw_A22), (x, y) -> SW_βc,
                       (sw_γ1, sw_γ2); C = SW_C, hooks = hooks, bcs = sw_zerobcs(mesh))
    op0 = sw_mixed_op(SW_βc)
    Ah = build_matrix(oph)
    A0 = build_matrix(op0)
    smh = SchwarzSmoother(oph; overlap = 2)
    for s in 1:5
        g = smh.subdomains[s].gdofs
        @test maximum(abs, subdomain_matrix(smh, s) .- Ah[g, g]) ≤ 1e-13 * maximum(abs, Ah)
    end
    # falsifier: the hooked W patch must DIFFER from the unhooked restriction (the
    # same ω=0.9 hook is visible in the W6-A global matrix at rel 1.2e-3 > 1e-4)
    g1 = smh.subdomains[1].gdofs
    @test maximum(abs, subdomain_matrix(smh, 1) .- A0[g1, g1]) > 1e-4 * maximum(abs, A0)
end

@testset "diagonal complex shift hook (schwarz-complex-shift seam)" begin
    op = sw_mixed_op(SW_βc)
    A = build_matrix(op)
    cache = SubdomainCache()
    sm0 = SchwarzSmoother(op; overlap = 2, cache = cache)
    sms = SchwarzSmoother(op; overlap = 2, shift = 0.7, cache = cache)
    # the shift is a NEW cache key (element, shape, epoch, s): all rebuilds, no hits
    @test cache.builds == 10
    @test cache.hits == 0
    rng = SWRNG(0x0000000073770001)
    for s in 1:5
        A0s = subdomain_matrix(sm0, s)
        Ass = subdomain_matrix(sms, s)
        # A_ss(shift) == A_ss(0) − i·s·Id BITWISE (spectre ExplicitInverse.hpp:249-254)
        @test Ass == A0s .- (im * 0.7) .* Matrix(I, size(A0s)...)
        # and the cached factorization solves the shifted dense restriction ≤ 1e-12
        g = sms.subdomains[s].gdofs
        rs = sw_cvec!(rng, length(g))
        ds = sms.fact[s] \ rs
        dref = (A[g, g] .- (im * 0.7) .* Matrix(I, length(g), length(g))) \ rs
        @test maximum(abs, ds .- dref) ≤ 1e-12 * maximum(abs, dref)
    end
end

@testset "cached LU: builds/hits counters, pass reuse, epoch invalidation" begin
    op = sw_mixed_op(SW_βc)
    n = op.mesh.ndofs
    cache = SubdomainCache()
    sm1 = SchwarzSmoother(op; overlap = 2, cache = cache)
    @test cache.builds == 5                          # one LU per subdomain
    @test cache.hits == 0
    # NODE GATE: smoothing passes build ZERO new LUs
    rng = SWRNG(0x0000000073770002)
    u = sw_cvec!(rng, n)
    b = sw_cvec!(rng, n)
    smooth!(u, sm1, b; passes = 2, threaded = false)
    @test cache.builds == 5
    # identically-shaped smoother on the same cache: all hits (AMR-reuse pattern)
    sm2 = SchwarzSmoother(op; overlap = 2, cache = cache)
    @test cache.builds == 5
    @test cache.hits == 5
    # coefficient epoch bump invalidates: full rebuild
    sm3 = SchwarzSmoother(op; overlap = 2, epoch = 1, cache = cache)
    @test cache.builds == 10
    @test cache.hits == 5
    # different overlap = different subdomain shape = different keys
    sm4 = SchwarzSmoother(op; overlap = 1, cache = cache)
    @test cache.builds == 15
    @test cache.hits == 5
end

@testset "2 passes damp high-frequency error > 10x (Poisson model)" begin
    mesh = sw_pmesh(5)
    op = MeshOperator(ComplexF64, mesh, (sw_one2, sw_one2), sw_zero2,
                      (sw_zero2, sw_zero2); C = SW_C, bcs = sw_zerobcs(mesh))
    sm = SchwarzSmoother(op; overlap = 2)
    b = zeros(ComplexF64, mesh.ndofs)
    # random top-Legendre-mode error: NODE GATE ‖A·u0‖/‖A·u2‖ > 10 (calibrated 77.2,
    # seed spread 69.7-78.3)
    u = sw_hfvec(mesh, SWRNG(0x0000000073770003); keep = 2)
    r0 = sw_resnorm(op, u)
    smooth!(u, sm, b; passes = 2, threaded = false)
    @test r0 / sw_resnorm(op, u) > 10
    # white-noise error, same gate (calibrated 90.0) — the WEIGHTING falsifier:
    # degenerate RAS weights (zero band weights, the SpecCore _wintr sign bug)
    # measure 7.6 on this model and FAIL this gate
    uw = sw_cvec!(SWRNG(0x0000000073770004), mesh.ndofs)
    r0w = sw_resnorm(op, uw)
    smooth!(uw, sm, b; passes = 2, threaded = false)
    @test r0w / sw_resnorm(op, uw) > 10
end

@testset "THREADED sweep bitwise == serial (JULIA_NUM_THREADS=4)" begin
    if Threads.nthreads() < 2
        @test_skip "threaded sweep == serial needs -t2"    # single-threaded run: no threaded/serial pair to compare
    else
        op = sw_mixed_op(SW_βH)                          # complex Helmholtz, overlaps live
        n = op.mesh.ndofs
        sm = SchwarzSmoother(op; overlap = 2)
        rng = SWRNG(0x0000000073770005)
        for _ in 1:3
            u0 = sw_cvec!(rng, n)
            b = sw_cvec!(rng, n)
            us = smooth!(copy(u0), sm, b; passes = 2, threaded = false)
            ut = smooth!(copy(u0), sm, b; passes = 2, threaded = true)
            @test ut == us                               # bitwise
        end
    end
end

@testset "complex Helmholtz patch solves == dense patch backslash + combine" begin
    op = sw_mixed_op(SW_βH)
    n = op.mesh.ndofs
    A = build_matrix(op)
    sm = SchwarzSmoother(op; overlap = 2)
    rng = SWRNG(0x0000000073770006)
    r = sw_cvec!(rng, n)
    # NODE GATE per-subdomain cached-LU solve == dense patch backslash ≤ 1e-12
    # (calibrated 0.0)
    for s in 1:5
        g = sm.subdomains[s].gdofs
        rs = r[g]
        dref = A[g, g] \ rs
        @test maximum(abs, (sm.fact[s] \ rs) .- dref) ≤ 1e-12 * maximum(abs, dref)
    end
    # one smoothing pass from u = 0 (so the pass residual is exactly b = r) equals the
    # test-side weighted combine of the dense patch solves ≤ 1e-12 (calibrated 0.0)
    u = smooth!(zeros(ComplexF64, n), sm, r; passes = 1, threaded = false)
    ref = zeros(ComplexF64, n)
    for s in 1:5
        g = sm.subdomains[s].gdofs
        w = sm.subdomains[s].weights
        ds = A[g, g] \ r[g]
        for k in eachindex(g)
            ref[g[k]] += w[k] * ds[k]
        end
    end
    @test maximum(abs, u .- ref) ≤ 1e-12 * maximum(abs, ref)
end

end # testset schwarz
