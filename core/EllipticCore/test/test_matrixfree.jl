# test_matrixfree.jl — grounded gates for src/matrixfree.jl (DAG node: matrix-free-operator,
# elliptic-ssf GLOBAL_DAG row 15 + TOPOLOGICAL_SCHEDULE W6-A, CRITICAL PATH). Written BEFORE
# the implementation (TDD red run captured in evidence/matrix-free-operator/pkg_test_red.log).
#
# Target: the OWNED global composition that waves 4-5 kept test-side — Mesh2D-level
# connectivity containers (elements, per-element orders, global nodal layout, neighbor/
# mortar topology incl. 2:1 and cross-block connections) and the two-pass sweep
# (pass 1: every element builds its ElementState and face traces are sender-projected to
# the interface mortars; pass 2: receivers apply the ModifyBoundaryData hooks on the RAW
# received mortar data, evaluate the eq:numflux mortar algebra with the pairwise (σ, Am)
# constructed IDENTICALLY on both sides — the W5-A contract: Am evaluated ONCE per
# interface mortar and shared, σ from mirrored Mortar constructor args — and lift),
# composing dg_operator + dg_penalty + mortars + boundary_hooks into one matrix-free
# A(u) over the global nodal vector layout. EVERY matrix application routes through the
# LINEARIZED hook/ghost seam (W5-B contract); fixed data reach the RHS only via
# operator-on-zero b -= A_full(0). Conventions N1-N7 (dg_operator.jl) inherited verbatim.
#
# Oracle gates (worker node contract; constants calibrated PRE-RED in scratchpad
# calibrate_matrixfree.jl from the sealed wave-2..5 layers plus a THROWAWAY prototype of
# the composition — NOT the implementation):
#   1. HEADLINE — build_matrix(A) column probes == an INDEPENDENT quadrature-level dense
#      assembly (eq:numflux mortar algebra from the wave-3 matrices + refelement
#      primitives; topology hand-coded) on a small mixed h/p mesh with one 2:1 interface
#      and two p-jumps, variable 𝔄, complex β, nonzero γ: rel ≤ 1e-13 (calibrated 4.2e-16).
#   2. complex linearity A(αu + v) = αA(u) + A(v) rel ≤ 1e-13 (calibrated 2.4e-16).
#   3. direct-solve bypass == GMRES(A): manufactured Poisson |x_direct − x_gmres|∞ ≤ 1e-10
#      (calibrated 1.5e-14, 66 iterations) and complex Helmholtz β = −25 − 0.4i ≤ 1e-10
#      (calibrated 6.8e-14, 66 iterations).
#   4. hook routing falsifiers BOTH ways: a LinearizedHook is VISIBLE in build_matrix
#      (rel diff 1.2e-3 calibrated > 1e-4 gate) while operator_on_zero stays exactly 0;
#      a FixedDataHook leaves build_matrix BITWISE unchanged while operator_on_zero equals
#      the W5-B per-face closed form rel ≤ 1e-13 (calibrated 0.0). End-to-end twins of the
#      W5-B node gates THROUGH the operator: prescribed-jump solve L2 ≤ 1e-10 (calibrated
#      3.4e-15), composed junction [u]=J, [u′]=K−iω⟨⟨u⟩⟩ solve L2 ≤ 1e-10 (calibrated
#      4.9e-15), and dropping the linearized stage leaves L2 > 1e-3 (calibrated 4.9e-2).
#   5. no-hook/no-mortar meshes reproduce the wave-4 (test_dg_penalty) and W5-A
#      (test_mortars 2:1) test-side sweeps BITWISE (calibrated exactly true, 3 random
#      complex vectors each).
#   6. THREADED apply (JULIA_NUM_THREADS=4) bitwise == serial on 3 random complex vectors
#      (element-parallel two-pass; no cross-element reductions).
#   7. cross-block closure (mesh2d-blocks caveat): a two-block mesh joined by a
#      BlockConnection (identical transverse maps, exact endpoint continuity) with a
#      cross-block 2:1 interface builds the SAME operator as the single-block twin mesh,
#      rel ≤ 1e-13 (calibrated 0.0 — bitwise); SIPG symmetry ≤ 1e-12 (calibrated 2.5e-16),
#      min sym eig > 0.05 (calibrated 0.263); transverse 2:1 violations and map/endpoint
#      mismatches across connections are construction errors.
#   8. BigFloat(256) pass-through: exactly-represented-polynomial residual == M·S rel
#      ≤ 1e-60 (calibrated 1.0e-75); direct_solve error ≤ 1e-60 (calibrated 2.1e-77).
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test files,
# ---- distinct name/seed) ---------------------------------------------------------------
mutable struct MFRNG
    s::UInt64
end
mf_unit!(r::MFRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
mf_vec!(r::MFRNG, n) = [mf_unit!(r) for _ in 1:n]
mf_cvec!(r::MFRNG, n) = [complex(mf_unit!(r), mf_unit!(r)) for _ in 1:n]

const MF_SIDES = (-1, 1)
const MF_C = 1.5                                     # CRP main.tex:595 penalty parameter

mf_one2(x, y) = 1.0
mf_zero2(x, y) = 0.0
mf_uex(x, y) = sinpi(x) * sinpi(y)                   # zero Dirichlet on [0,1]²
mf_S(x, y) = 2 * pi^2 * mf_uex(x, y)                 # −∇²u = S
mf_A11(x, y) = 1 + 0.3x + 0.1y^2                     # variable diagonal flux matrix
mf_A22(x, y) = 1 + 0.2y
const MF_βc = 0.3 + 0.2im
mf_γ1(x, y) = 0.05
mf_γ2(x, y) = -0.03

# ---- the mixed h/p mesh of the headline gate: one 2:1 interface (W | E11,E12) and two
# ---- p-jumps (E11|E21, E12|E22 at orders 4 vs 5) on the affine unit block ---------------
mf_blk() = Block(91, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
const MFW   = ElementId(91, SegmentId(1, 0), SegmentId(0, 0))
const MFE11 = ElementId(91, SegmentId(2, 2), SegmentId(1, 0))
const MFE12 = ElementId(91, SegmentId(2, 2), SegmentId(1, 1))
const MFE21 = ElementId(91, SegmentId(2, 3), SegmentId(1, 0))
const MFE22 = ElementId(91, SegmentId(2, 3), SegmentId(1, 1))
mf_els() = [MFW, MFE11, MFE12, MFE21, MFE22]
mf_ords() = Dict(MFW => (4, 4), MFE11 => (4, 4), MFE12 => (4, 4),
                 MFE21 => (5, 5), MFE22 => (5, 5))

# hand-coded topology of that mesh, INDEPENDENT of the Mesh2D containers under test:
# interior interfaces listed once as (hi-side element, dim, [(lo-side element, mortar
# transverse segment)…]) and the boundary faces as (element, dim, sidx)
mf_ifaces() = [
    (MFW, 1, [(MFE11, SegmentId(1, 0)), (MFE12, SegmentId(1, 1))]),   # the 2:1 interface
    (MFE11, 1, [(MFE21, SegmentId(1, 0))]),                            # p-jump 4|5
    (MFE12, 1, [(MFE22, SegmentId(1, 1))]),                            # p-jump 4|5
    (MFE11, 2, [(MFE12, SegmentId(2, 2))]),                            # conforming
    (MFE21, 2, [(MFE22, SegmentId(2, 3))]),                            # conforming
]
mf_bdry() = [(MFW, 1, 1), (MFW, 2, 1), (MFW, 2, 2), (MFE11, 2, 1), (MFE12, 2, 2),
             (MFE21, 1, 2), (MFE21, 2, 1), (MFE22, 1, 2), (MFE22, 2, 2)]

# test-side global layout for a (ids, orders) list — must agree with the Mesh2D contract
function mf_off(ids, ords)
    off = Dict{Any,Int}()
    o = 0
    for id in ids
        off[id] = o
        o += prod(ords[id])
    end
    return off, o
end

# zero-Dirichlet ghost bcs on every boundary face of a mesh (through the public topology)
function mf_zerobcs(mesh)
    bcs = Dict{Any,Any}()
    for (el, d, sidx) in boundary_faces(mesh)
        nf = mesh.orders[mesh.index[el]][3-d]
        bcs[(el, d, sidx)] = DirichletGhost(zeros(ComplexF64, nf))
    end
    return bcs
end

# massive right-hand side b = M·S and mass-weighted L2 error, test-side per element
function mf_rhs(blk, ids, ords, off, n, S)
    b = zeros(ComplexF64, n)
    for id in ids
        em = element_matrices(Float64, blk, id, ords[id])
        Nx, Ny = em.N
        b[off[id]+1:off[id]+Nx*Ny] .=
            em.M * vec([S(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
    end
    return b
end
function mf_l2err(blk, ids, ords, off, u, ue)
    s = 0.0
    for id in ids
        em = element_matrices(Float64, blk, id, ords[id])
        Nx, Ny = em.N
        e = u[off[id]+1:off[id]+Nx*Ny] .-
            vec([ue(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
        s += real(dot(e, em.M * e))
    end
    return sqrt(s)
end

# ---- INDEPENDENT quadrature-level dense assembly (the headline oracle): the eq:numflux
# ---- mortar algebra written directly in the wave-3 matrices + refelement primitives.
# ---- P from barycentric interpolation at geometry-derived points, R = jac·W_f⁻¹PᵀW_m,
# ---- σ from the Vu:2024 law with PHYSICAL element extents, Am evaluated from the
# ---- coefficient callable at geometry-derived mortar points. No Mesh2D / MeshOperator /
# ---- Mortar / numerical_flux / strong_residual code paths. -----------------------------
function mfo_dphys(em, d)
    Nx, Ny = em.N
    if d == 1
        D = EllipticCore.refelement(Float64, Nx).D
        M = kron(Matrix(1.0I, Ny, Ny), D)
        for j in 1:Ny, i in 1:Nx
            M[tensor_lin(i, j, Nx), :] ./= em.jac[1][i]
        end
        return M
    else
        D = EllipticCore.refelement(Float64, Ny).D
        M = kron(D, Matrix(1.0I, Nx, Nx))
        for j in 1:Ny, i in 1:Nx
            M[tensor_lin(i, j, Nx), :] ./= em.jac[2][j]
        end
        return M
    end
end
mfo_avol(f, em) = vec([f(em.x[1][i], em.x[2][j]) for i in 1:em.N[1], j in 1:em.N[2]])
function mfo_facecoeff(blk, id, em, d, sidx, Afun)
    lo, hi = logical_interval(Float64, id, d)
    xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
    xt = em.x[3-d]
    return d == 1 ? [Afun[1](xf, t) for t in xt] : [Afun[2](t, xf) for t in xt]
end
function mfo_matrix(blk, ids, ords, off, n, ifaces, bdry, Afun, βfun, γfun)
    A = zeros(ComplexF64, n, n)
    ems = Dict(id => element_matrices(Float64, blk, id, ords[id]) for id in ids)
    Dp = Dict(id => (mfo_dphys(ems[id], 1), mfo_dphys(ems[id], 2)) for id in ids)
    rng(id) = off[id]+1:off[id]+prod(ords[id])
    # volume terms of eq:dgres_strong_matrix
    for id in ids
        em = ems[id]
        A[rng(id), rng(id)] .+= em.M * Diagonal(mfo_avol(βfun, em)) .+
            em.M * Diagonal(mfo_avol(γfun[1], em)) * Dp[id][1] .+
            em.M * Diagonal(mfo_avol(γfun[2], em)) * Dp[id][2] .-
            em.MD[1] * Diagonal(mfo_avol(Afun[1], em)) * Dp[id][1] .-
            em.MD[2] * Diagonal(mfo_avol(Afun[2], em)) * Dp[id][2]
    end
    # zero-Dirichlet ghost faces: −MLDT·diag(s𝔄f)·Iᵀ + 2σ·ML·diag(𝔄f)·Iᵀ, receiver σ
    for (id, d, sidx) in bdry
        em = ems[id]
        s = MF_SIDES[sidx]
        fm = face_matrices(em, d, s)
        nvol, nf = size(fm.I)
        p = nvol ÷ nf - 1
        lo, hi = logical_interval(Float64, id, d)
        h = map_forward(blockmap(blk, d), hi) - map_forward(blockmap(blk, d), lo)
        σ = MF_C * (p + 1)^2 / h
        Af = mfo_facecoeff(blk, id, em, d, sidx, Afun)
        Ti = transpose(fm.I)
        A[rng(id), rng(id)] .+= .-fm.MLDT * Diagonal(s .* Af) * Ti .+
                                2σ .* (fm.ML * Diagonal(Af) * Ti)
    end
    # interior receiver faces: the once-per-face self terms −MLDT·diag(s𝔄f)·Iᵀ + ML·Bi
    for (elhi, d, lows) in ifaces
        recv = [(elhi, +1, 2)]
        append!(recv, [(nb, -1, 1) for (nb, _) in lows])
        for (id, s, sidx) in recv
            em = ems[id]
            fm = face_matrices(em, d, s)
            Af = mfo_facecoeff(blk, id, em, d, sidx, Afun)
            Ti = transpose(fm.I)
            Bi = s .* (Ti * Diagonal(mfo_avol(Afun[d == 1 ? 1 : 2], em)) * Dp[id][d])
            A[rng(id), rng(id)] .+= .-fm.MLDT * Diagonal(s .* Af) * Ti .+ fm.ML * Bi
        end
    end
    # interface mortars: the mortar-sum parts of Ψ*, (nF)* on BOTH sides
    for (elhi, d, lows) in ifaces, (nbel, span) in lows
        e, f = elhi, nbel
        em_e, em_f = ems[e], ems[f]
        fm_e = face_matrices(em_e, d, +1)
        fm_f = face_matrices(em_f, d, -1)
        nf_e = size(fm_e.I, 2); nf_f = size(fm_f.I, 2)
        Nm = max(nf_e, nf_f)
        t = 3 - d
        lom, him = logical_interval(Float64, span)
        elm = EllipticCore.refelement(Float64, Nm)
        ζm = [(lom * (1 - x) + him * (1 + x)) / 2 for x in elm.nodes]
        tmap = blockmap(blk, t)
        xtm = [map_forward(tmap, z) for z in ζm]
        function PR(id, nf)
            elf = EllipticCore.refelement(Float64, nf)
            lof, hif = logical_interval(Float64, id, t)
            τ = [(2z - (lof + hif)) / (hif - lof) for z in ζm]
            P = interpolation_matrix(elf, τ)
            jac = (him - lom) / (hif - lof)
            return P, jac .* (inv(Diagonal(elf.weights)) * transpose(P) * Diagonal(elm.weights))
        end
        P_e, R_e = PR(e, nf_e)
        P_f, R_f = PR(f, nf_f)
        dmap = blockmap(blk, d)
        lo_e, hi_e = logical_interval(Float64, e, d)
        lo_f, hi_f = logical_interval(Float64, f, d)
        h_e = map_forward(dmap, hi_e) - map_forward(dmap, lo_e)
        h_f = map_forward(dmap, hi_f) - map_forward(dmap, lo_f)
        σ = MF_C * (max(em_e.N[d] - 1, em_f.N[d] - 1) + 1)^2 / min(h_e, h_f)
        xfc = map_forward(dmap, hi_e)
        Am = d == 1 ? [Afun[1](xfc, tt) for tt in xtm] : [Afun[2](tt, xfc) for tt in xtm]
        DAm = Diagonal(Am)
        Ti_e = transpose(fm_e.I); Ti_f = transpose(fm_f.I)
        Bi_e = (+1) .* (Ti_e * Diagonal(mfo_avol(Afun[d == 1 ? 1 : 2], em_e)) * Dp[e][d])
        Bi_f = (-1) .* (Ti_f * Diagonal(mfo_avol(Afun[d == 1 ? 1 : 2], em_f)) * Dp[f][d])
        Af_e = mfo_facecoeff(blk, e, em_e, d, 2, Afun)
        Af_f = mfo_facecoeff(blk, f, em_f, d, 1, Afun)
        re, rf = rng(e), rng(f)
        A[re, re] .+= fm_e.MLDT * Diagonal(+1 .* Af_e) * (R_e * P_e * Ti_e) ./ 2 .-
                      fm_e.ML * (R_e * P_e * Bi_e ./ 2) .+ σ .* (fm_e.ML * (R_e * DAm * P_e) * Ti_e)
        A[re, rf] .+= fm_e.MLDT * Diagonal(+1 .* Af_e) * (R_e * P_f * Ti_f) ./ 2 .+
                      fm_e.ML * (R_e * P_f * Bi_f) ./ 2 .- σ .* (fm_e.ML * (R_e * DAm * P_f) * Ti_f)
        A[rf, rf] .+= fm_f.MLDT * Diagonal(-1 .* Af_f) * (R_f * P_f * Ti_f) ./ 2 .-
                      fm_f.ML * (R_f * P_f * Bi_f ./ 2) .+ σ .* (fm_f.ML * (R_f * DAm * P_f) * Ti_f)
        A[rf, re] .+= fm_f.MLDT * Diagonal(-1 .* Af_f) * (R_f * P_e * Ti_e) ./ 2 .+
                      fm_f.ML * (R_f * P_e * Bi_e) ./ 2 .- σ .* (fm_f.ML * (R_f * DAm * P_e) * Ti_e)
    end
    return A
end

@testset "matrixfree (matrix-free-operator)" begin

@testset "Mesh2D containers + topology validation" begin
    blk = mf_blk()
    mesh = Mesh2D(blk, mf_els(), mf_ords())
    @test mesh.ndofs == 3 * 16 + 2 * 25
    @test mesh.elements == mf_els()                  # layout order preserved
    @test mesh.offsets == [0, 16, 32, 48, 73]        # documented global nodal layout
    @test mesh.orders == [(4, 4), (4, 4), (4, 4), (5, 5), (5, 5)]
    @test Set(boundary_faces(mesh)) == Set(mf_bdry())
    # uniform-order convenience constructor
    ids4 = vec([ElementId(91, SegmentId(2, ix), SegmentId(2, iy)) for ix in 0:3, iy in 0:3])
    m4 = Mesh2D(blk, ids4, (5, 5))
    @test m4.ndofs == 16 * 25
    @test length(boundary_faces(m4)) == 16           # 4 faces × 4 boundary edges

    # within-block 2:1 violation (transverse level jump of 2 across a face)
    bad = [ElementId(91, SegmentId(1, 0), SegmentId(0, 0)),
           ElementId(91, SegmentId(2, 2), SegmentId(2, 0))]
    @test_throws ArgumentError Mesh2D(blk, bad, (4, 4))
    # tiling hole: E12 missing leaves E11's north face uncovered (interior, no neighbor)
    @test_throws ArgumentError Mesh2D(blk, [MFW, MFE11, MFE21, MFE22], (4, 4))
    # duplicate element / unknown block id
    @test_throws ArgumentError Mesh2D(blk, [MFW, MFW], (4, 4))
    @test_throws ArgumentError Mesh2D(blk, [ElementId(77, SegmentId(0, 0), SegmentId(0, 0))],
                                      (4, 4))

    # MeshOperator seam validation: bcs must cover EXACTLY the boundary faces; hooks
    # bind only to interior faces
    A = (mf_one2, mf_one2)
    γ = (mf_zero2, mf_zero2)
    bcs = mf_zerobcs(mesh)
    badbcs = copy(bcs)
    delete!(badbcs, (MFW, 1, 1))
    @test_throws ArgumentError MeshOperator(ComplexF64, mesh, A, mf_zero2, γ; bcs = badbcs)
    extrabcs = copy(bcs)
    extrabcs[(MFW, 1, 2)] = DirichletGhost(zeros(ComplexF64, 4))     # interior face
    @test_throws ArgumentError MeshOperator(ComplexF64, mesh, A, mf_zero2, γ; bcs = extrabcs)
    @test_throws ArgumentError MeshOperator(ComplexF64, mesh, A, mf_zero2, γ; bcs = bcs,
        hooks = Dict((MFW, 1, 1) => NoHook()))                       # boundary face
end

@testset "HEADLINE: build_matrix == independent dense assembly (mixed h/p mesh)" begin
    blk = mf_blk()
    ids = mf_els()
    ords = mf_ords()
    off, n = mf_off(ids, ords)
    mesh = Mesh2D(blk, ids, ords)
    @test [off[id] for id in ids] == mesh.offsets    # test-side layout == container layout
    op = MeshOperator(ComplexF64, mesh, (mf_A11, mf_A22), (x, y) -> MF_βc, (mf_γ1, mf_γ2);
                      C = MF_C, bcs = mf_zerobcs(mesh))
    Aprobe = build_matrix(op)
    @test size(Aprobe) == (n, n)
    Aora = mfo_matrix(blk, ids, ords, off, n, mf_ifaces(), mf_bdry(),
                      (mf_A11, mf_A22), (x, y) -> MF_βc, (mf_γ1, mf_γ2))
    sc = maximum(abs, Aora)
    # NODE GATE 1e-13 (calibrated 4.2e-16)
    @test maximum(abs, Aprobe .- Aora) ≤ 1e-13 * sc
end

@testset "complex linearity + affine hook decomposition" begin
    blk = mf_blk()
    mesh = Mesh2D(blk, mf_els(), mf_ords())
    op = MeshOperator(ComplexF64, mesh, (mf_A11, mf_A22), (x, y) -> MF_βc, (mf_γ1, mf_γ2);
                      C = MF_C, bcs = mf_zerobcs(mesh))
    n = mesh.ndofs
    rng = MFRNG(0x00000000feed0071)
    u = mf_cvec!(rng, n)
    v = mf_cvec!(rng, n)
    α = 0.37 - 1.21im
    yu = apply!(zeros(ComplexF64, n), op, u)
    yv = apply!(zeros(ComplexF64, n), op, v)
    yav = apply!(zeros(ComplexF64, n), op, α .* u .+ v)
    # NODE GATE A(αu + v) == αA(u) + A(v) at 1e-13 (calibrated 2.4e-16)
    @test maximum(abs, yav .- α .* yu .- yv) ≤ 1e-13 * maximum(abs, yav)

    # with fixed hooks on the 2:1 interface faces: A_full(u) = A_lin(u) + A_full(0)
    # (rides the per-mortar hook seam; deltas live at the mortar points)
    rngh = MFRNG(0x00000000feed0072)
    hooks = Dict((MFW, 1, 2) => FixedDataHook(mf_cvec!(rngh, 4), mf_cvec!(rngh, 4)),
                 (MFE11, 1, 1) => FixedDataHook(mf_cvec!(rngh, 4), mf_cvec!(rngh, 4)),
                 (MFE12, 1, 1) => FixedDataHook(mf_cvec!(rngh, 4), mf_cvec!(rngh, 4)))
    opf = MeshOperator(ComplexF64, mesh, (mf_A11, mf_A22), (x, y) -> MF_βc, (mf_γ1, mf_γ2);
                       C = MF_C, hooks = hooks, bcs = mf_zerobcs(mesh))
    R0 = operator_on_zero(opf)
    @test maximum(abs, R0) > 0.1                     # genuinely nonzero fixed sources
    Ru = apply_full!(zeros(ComplexF64, n), opf, u)
    Lu = apply!(zeros(ComplexF64, n), opf, u)
    @test maximum(abs, Ru .- Lu .- R0) ≤ 1e-13 * maximum(abs, Ru)
end

@testset "solve gates: direct-solve bypass == GMRES (Poisson + complex Helmholtz)" begin
    blk = mf_blk()
    ids = mf_els()
    ords = mf_ords()
    off, n = mf_off(ids, ords)
    mesh = Mesh2D(blk, ids, ords)
    γ0 = (mf_zero2, mf_zero2)
    # manufactured Poisson
    opp = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), mf_zero2, γ0;
                       C = MF_C, bcs = mf_zerobcs(mesh))
    b = mf_rhs(blk, ids, ords, off, n, mf_S)
    xd = direct_solve(opp, b)
    res = gmres(opp, b; abstol = 1e-12, reltol = 0.0, maxiter = 2000, restart = 2000)
    @test res.converged
    @test res.iterations ≤ 300                       # calibrated 66
    # NODE GATE bypass == GMRES at 1e-10 (calibrated 1.5e-14)
    @test maximum(abs, xd .- res.x) ≤ 1e-10
    @test mf_l2err(blk, ids, ords, off, xd, mf_uex) < 5e-2     # calibrated 8.2e-3
    # complex Helmholtz β = −k² − iε
    βH = -25.0 - 0.4im
    SH(x, y) = (2 * pi^2 + βH) * mf_uex(x, y)
    oph = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), (x, y) -> βH, γ0;
                       C = MF_C, bcs = mf_zerobcs(mesh))
    bh = mf_rhs(blk, ids, ords, off, n, SH)
    @test eltype(bh) == ComplexF64
    xdh = direct_solve(oph, bh)
    resh = gmres(oph, bh; abstol = 1e-12, reltol = 0.0, maxiter = 4000, restart = 4000)
    @test resh.converged
    @test resh.iterations ≤ 500                      # calibrated 66
    @test maximum(abs, xdh .- resh.x) ≤ 1e-10        # calibrated 6.8e-14
    @test mf_l2err(blk, ids, ords, off, xdh, mf_uex) < 1e-1    # calibrated 1.1e-2
end

@testset "hook routing: linearized in the matrix, fixed only in b (falsifiers)" begin
    blk = Block(92, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    N = 6
    L = ElementId(92, SegmentId(1, 0), SegmentId(0, 0))
    R = ElementId(92, SegmentId(1, 1), SegmentId(0, 0))
    mesh = Mesh2D(blk, [L, R], (N, N))
    γ0 = (mf_zero2, mf_zero2)
    mkop(; hooks = Dict(), bcs = mf_zerobcs(mesh)) =
        MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), (x, y) -> MF_βc, γ0;
                     C = MF_C, hooks = hooks, bcs = bcs)
    op0 = mkop()
    M0 = build_matrix(op0)
    # (a) LinearizedHook (vtu-style −iω⟨⟨ψ⟩⟩ flux jump) is VISIBLE in the matrix and
    # leaves operator_on_zero exactly zero
    ω = 0.9
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    opl = mkop(hooks = Dict((L, 1, 2) => LinearizedHook(fω), (R, 1, 1) => LinearizedHook(fω)))
    Ml = build_matrix(opl)
    @test maximum(abs, Ml .- M0) > 1e-4 * maximum(abs, M0)     # calibrated 1.2e-3
    @test all(iszero, operator_on_zero(opl))                   # calibrated exactly 0
    # (b) FixedDataHook leaves the matrix BITWISE unchanged (linearized seam) and lands
    # exactly in operator_on_zero == the W5-B per-face closed form
    J, K = 0.7 + 0.0im, -1.3 + 0.0im
    onesf = ones(N)
    opf = mkop(hooks = Dict((L, 1, 2) => FixedDataHook(-J .* onesf, K .* onesf),
                            (R, 1, 1) => FixedDataHook(+J .* onesf, K .* onesf)))
    @test build_matrix(opf) == M0                              # calibrated bitwise-equal
    R0 = operator_on_zero(opf)
    direct = zeros(ComplexF64, mesh.ndofs)
    for (id, s, ψsign) in ((L, +1, -1), (R, -1, +1))
        em = element_matrices(Float64, blk, id, (N, N))
        fm = face_matrices(em, 1, s)
        nvol, nf = size(fm.I)
        σ = penalty_sigma(MF_C, nvol ÷ nf - 1, nvol ÷ nf - 1,
                          perpendicular_size(fm), perpendicular_size(fm))
        Δψ = ψsign * J .* onesf
        Δφ = K .* onesf
        o = id == L ? 0 : N * N
        direct[o+1:o+N*N] .+= fm.MLDT * (s .* Δψ ./ 2) .- fm.ML * (σ .* Δψ .- Δφ ./ 2)
    end
    @test maximum(abs, R0) > 1                                 # nonzero fixed sources
    @test maximum(abs, R0 .- direct) ≤ 1e-13 * maximum(abs, R0)  # calibrated 0.0
end

@testset "W5-B node-gate twins THROUGH the operator (jump + composed junction)" begin
    blk = Block(93, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    N = 6
    xm = 0.5
    L = ElementId(93, SegmentId(1, 0), SegmentId(0, 0))
    R = ElementId(93, SegmentId(1, 1), SegmentId(0, 0))
    mesh = Mesh2D(blk, [L, R], (N, N))
    ids = [L, R]
    ords = Dict(L => (N, N), R => (N, N))
    off, n = mf_off(ids, ords)
    onesf = ones(N)
    γ0 = (mf_zero2, mf_zero2)
    # Dirichlet data of uex at the x-boundaries, Neumann-zero at the y-boundaries
    function dirneu_bcs(uexf)
        bcs = Dict{Any,Any}()
        for (el, d, sidx) in boundary_faces(mesh)
            em = element_matrices(Float64, blk, el, (N, N))
            lo, hi = logical_interval(Float64, el, d)
            xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
            xt = em.x[3-d]
            data = d == 1 ? [uexf(xf, t) for t in xt] : [uexf(t, xf) for t in xt]
            bcs[(el, d, sidx)] = d == 1 ? DirichletGhost(ComplexF64.(data)) :
                                          NeumannGhost(zeros(ComplexF64, N))
        end
        return bcs
    end

    # (i) prescribed-jump solve [Ψ]=J, [∂Ψ]=K (piecewise cubic, exact): sign table from
    # the W5-B header — Δψ OPPOSITE (−J left / +J right), Δφ SAME (+K)
    J, K = 0.7, -1.3
    uL(x) = x^3 - x + 0.25
    uR(x) = uL(x) + J + K * (x - xm)
    uexj(x, y) = x < xm ? uL(x) : uR(x)
    Sj(x, y) = -6x
    opj = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), mf_zero2, γ0; C = MF_C,
                       hooks = Dict((L, 1, 2) => FixedDataHook(-J .* onesf, K .* onesf),
                                    (R, 1, 1) => FixedDataHook(+J .* onesf, K .* onesf)),
                       bcs = dirneu_bcs(uexj))
    bj = mf_rhs(blk, ids, ords, off, n, Sj)
    R0j = operator_on_zero(opj)
    xj = direct_solve(opj, bj .- R0j)
    # NODE GATE 1e-10 (calibrated 3.4e-15)
    @test mf_l2err(blk, ids, ords, off, xj, uexj) ≤ 1e-10

    # (ii) composed junction [u] = J, [u′] = K − iω⟨⟨u⟩⟩ (linearized stage on the RAW
    # received data BEFORE the fixed shifts) — manufactured complex piecewise cubic
    ω = 0.9
    Jc = 0.35 - 0.15im
    Kc = 0.6 + 0.25im
    uLc(x) = (0.8 + 0.5im) * x^3 - (0.4 - 0.2im) * x + 0.3
    K̃ = Kc - im * ω * (uLc(xm) + Jc / 2)
    uRc(x) = uLc(x) + Jc + K̃ * (x - xm)
    uexc(x, y) = x < xm ? uLc(x) : uRc(x)
    Sc(x, y) = -6 * (0.8 + 0.5im) * x + MF_βc * uexc(x, y)
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    hooksc = Dict(
        (L, 1, 2) => ComposedHook(LinearizedHook(fω), FixedDataHook(-Jc .* onesf, Kc .* onesf)),
        (R, 1, 1) => ComposedHook(LinearizedHook(fω), FixedDataHook(+Jc .* onesf, Kc .* onesf)))
    opc = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), (x, y) -> MF_βc, γ0; C = MF_C,
                       hooks = hooksc, bcs = dirneu_bcs(uexc))
    bc = mf_rhs(blk, ids, ords, off, n, Sc)
    R0c = operator_on_zero(opc)
    xc = direct_solve(opc, bc .- R0c)
    @test mf_l2err(blk, ids, ords, off, xc, uexc) ≤ 1e-10      # calibrated 4.9e-15
    # bypass == GMRES through the SAME operator callable (calibrated 2.5e-14, 53 its)
    resc = gmres(opc, bc .- R0c; abstol = 1e-13, reltol = 0.0, maxiter = 500, restart = 500)
    @test resc.converged
    @test maximum(abs, xc .- resc.x) ≤ 1e-10
    # FALSIFIER: dropping the linearized stage leaves the matrix without −iω⟨⟨ψ⟩⟩ and the
    # solve misses at O(1) (calibrated 4.9e-2 ≫ 1e-3); R0 is unchanged by the drop
    opx = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), (x, y) -> MF_βc, γ0; C = MF_C,
                       hooks = Dict((L, 1, 2) => FixedDataHook(-Jc .* onesf, Kc .* onesf),
                                    (R, 1, 1) => FixedDataHook(+Jc .* onesf, Kc .* onesf)),
                       bcs = dirneu_bcs(uexc))
    R0x = operator_on_zero(opx)
    @test maximum(abs, R0x .- R0c) ≤ 1e-13 * maximum(abs, R0c)
    xx = direct_solve(opx, bc .- R0x)
    @test mf_l2err(blk, ids, ords, off, xx, uexc) > 1e-3
end

@testset "bitwise reproduction of the wave-4 and W5-A test-side sweeps" begin
    # (i) conforming 4×4 mesh — the dgp_make_apply twin of test_dg_penalty.jl
    blk = Block(94, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    N = 5
    ids = vec([ElementId(94, SegmentId(2, ix), SegmentId(2, iy)) for ix in 0:3, iy in 0:3])
    ems = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in ids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, (mf_one2, mf_one2), mf_zero2,
                                         (mf_zero2, mf_zero2)) for id in ids)
    off = Dict(id => (k - 1) * N * N for (k, id) in enumerate(ids))
    n = length(ids) * N * N
    function w4apply!(y, u)
        sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+N*N))
                   for id in ids)
        for id in ids
            ext = [[begin
                        nb = neighbor_id(id, d, MF_SIDES[sidx])
                        nb === nothing ? zero_dirichlet_ghost(sts[id].faces[d][sidx]) :
                            exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1])
                    end for sidx in 1:2] for d in 1:2]
            y[off[id]+1:off[id]+N*N] .= strong_residual(ems[id], cs[id], sts[id], ext,
                                                        PenaltyFlux(MF_C))
        end
        return y
    end
    mesh = Mesh2D(blk, ids, (N, N))
    op = MeshOperator(ComplexF64, mesh, (mf_one2, mf_one2), mf_zero2, (mf_zero2, mf_zero2);
                      C = MF_C, bcs = mf_zerobcs(mesh))
    rng = MFRNG(0x00000000feed0074)
    for _ in 1:3
        u = mf_cvec!(rng, n)
        yref = w4apply!(zeros(ComplexF64, n), u)
        y = apply!(zeros(ComplexF64, n), op, u; threaded = false)
        @test y == yref                              # bitwise (calibrated true)
        yfull = apply_full!(zeros(ComplexF64, n), op, u; threaded = false)
        @test yfull == yref                          # zero-data ghosts: full == linearized
    end

    # (ii) the W5-A 2:1 mesh — the mor_apply21 twin of test_mortars.jl (uniform N = 5)
    blk2 = mf_blk()
    ids2 = mf_els()
    N2 = 5
    ems2 = Dict(id => element_matrices(Float64, blk2, id, (N2, N2)) for id in ids2)
    cs2 = Dict(id => element_coefficients(ems2[id], blk2, id, (mf_one2, mf_one2), mf_zero2,
                                          (mf_zero2, mf_zero2)) for id in ids2)
    off2 = Dict(id => (k - 1) * N2 * N2 for (k, id) in enumerate(ids2))
    n2 = length(ids2) * N2 * N2
    fmW = face_matrices(ems2[MFW], 1, +1)
    fmE = face_matrices(ems2[MFE11], 1, -1)
    Am = ones(N2)
    mWlo = Mortar(fmW, fmE, MORTAR_LOWER, Am; C = MF_C)
    mWup = Mortar(fmW, fmE, MORTAR_UPPER, Am; C = MF_C)
    mE = Mortar(fmE, fmW, MORTAR_FULL, Am; C = MF_C)
    pen = PenaltyFlux(MF_C)
    flux = Dict(MFW   => PerFaceFlux(((pen, MortarFlux([mWlo, mWup])), (pen, pen))),
                MFE11 => PerFaceFlux(((MortarFlux([mE]), pen), (pen, pen))),
                MFE12 => PerFaceFlux(((MortarFlux([mE]), pen), (pen, pen))),
                MFE21 => PerFaceFlux(((pen, pen), (pen, pen))),
                MFE22 => PerFaceFlux(((pen, pen), (pen, pen))))
    conf = [(1, MFE11, MFE21), (1, MFE12, MFE22), (2, MFE11, MFE12), (2, MFE21, MFE22)]
    function ext21(sts, id)
        function confnb(d, sidx)
            for (dd, lo, hi) in conf
                dd == d || continue
                id == lo && sidx == 2 && return sts[hi].faces[d][1]
                id == hi && sidx == 1 && return sts[lo].faces[d][2]
            end
            return nothing
        end
        return [[begin
                     if id == MFW && d == 1 && sidx == 2
                         lo = project_to_mortar(mE.P, sts[MFE11].faces[1][1])
                         up = project_to_mortar(mE.P, sts[MFE12].faces[1][1])
                         FaceState(vcat(lo.psi, up.psi), vcat(lo.nflux, up.nflux))
                     elseif (id == MFE11 || id == MFE12) && d == 1 && sidx == 1
                         P = id == MFE11 ? mWlo.P : mWup.P
                         project_to_mortar(P, sts[MFW].faces[1][2])
                     else
                         nb = confnb(d, sidx)
                         nb === nothing ? zero_dirichlet_ghost(sts[id].faces[d][sidx]) :
                             exterior_state(nb)
                     end
                 end for sidx in 1:2] for d in 1:2]
    end
    function w5apply!(y, u)
        sts = Dict(id => element_state(ems2[id], cs2[id],
                                       view(u, off2[id]+1:off2[id]+N2*N2)) for id in ids2)
        for id in ids2
            y[off2[id]+1:off2[id]+N2*N2] .= strong_residual(ems2[id], cs2[id], sts[id],
                                                            ext21(sts, id), flux[id])
        end
        return y
    end
    mesh2 = Mesh2D(blk2, ids2, (N2, N2))
    op2 = MeshOperator(ComplexF64, mesh2, (mf_one2, mf_one2), mf_zero2,
                       (mf_zero2, mf_zero2); C = MF_C, bcs = mf_zerobcs(mesh2))
    for _ in 1:3
        u = mf_cvec!(rng, n2)
        yref = w5apply!(zeros(ComplexF64, n2), u)
        y = apply!(zeros(ComplexF64, n2), op2, u; threaded = false)
        @test y == yref                              # bitwise (calibrated true)
    end
end

@testset "threaded apply (JULIA_NUM_THREADS=4) bitwise == serial" begin
    # the two-pass sweep is element-parallel with no cross-element reductions, so the
    # threaded apply must be BITWISE identical to serial; the verification command runs
    # the suite under JULIA_NUM_THREADS=4
    @test Threads.nthreads() ≥ 2
    blk = mf_blk()
    mesh = Mesh2D(blk, mf_els(), mf_ords())
    ω = 0.9
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    op = MeshOperator(ComplexF64, mesh, (mf_A11, mf_A22), (x, y) -> MF_βc, (mf_γ1, mf_γ2);
                      C = MF_C,
                      hooks = Dict((MFE11, 2, 2) => LinearizedHook(fω),
                                   (MFE12, 2, 1) => LinearizedHook(fω)),
                      bcs = mf_zerobcs(mesh))
    n = mesh.ndofs
    rng = MFRNG(0x00000000feed0075)
    for _ in 1:3
        u = mf_cvec!(rng, n)
        ys = apply!(zeros(ComplexF64, n), op, u; threaded = false)
        yt = apply!(zeros(ComplexF64, n), op, u; threaded = true)
        @test yt == ys                               # bitwise
        fs = apply_full!(zeros(ComplexF64, n), op, u; threaded = false)
        ft = apply_full!(zeros(ComplexF64, n), op, u; threaded = true)
        @test ft == fs
    end
end

@testset "cross-block: connection == single-block twin + 2:1 validation" begin
    # two blocks joined at x = 0.5 (identical transverse maps, exact endpoint continuity)
    # with a cross-block 2:1 interface; the single-block twin tiles the same physical
    # mesh, so the operators must agree (calibrated bitwise-equal; gate 1e-13)
    Lb = Block(95, LinearMap1D(0.0, 0.5), LinearMap1D(0.0, 1.0))
    Rb = Block(96, LinearMap1D(0.5, 1.0), LinearMap1D(0.0, 1.0))
    elL = ElementId(95, SegmentId(0, 0), SegmentId(0, 0))
    elR1 = ElementId(96, SegmentId(0, 0), SegmentId(1, 0))
    elR2 = ElementId(96, SegmentId(0, 0), SegmentId(1, 1))
    mesh_x = Mesh2D([Lb, Rb], [elL, elR1, elR2], (5, 5);
                    connections = [BlockConnection(95, 1, 96)])
    @test mesh_x.ndofs == 75
    @test length(boundary_faces(mesh_x)) == 7        # the shared face is interior:
    # left 1 (elL spans full y) + right 2 + bottom 2 + top 2 perimeter faces
    op_x = MeshOperator(ComplexF64, mesh_x, (mf_A11, mf_A22), (x, y) -> MF_βc,
                        (mf_zero2, mf_zero2); C = MF_C, bcs = mf_zerobcs(mesh_x))
    Ax = build_matrix(op_x)
    Sb = Block(97, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    mesh_s = Mesh2D(Sb, [ElementId(97, SegmentId(1, 0), SegmentId(0, 0)),
                         ElementId(97, SegmentId(1, 1), SegmentId(1, 0)),
                         ElementId(97, SegmentId(1, 1), SegmentId(1, 1))], (5, 5))
    op_s = MeshOperator(ComplexF64, mesh_s, (mf_A11, mf_A22), (x, y) -> MF_βc,
                        (mf_zero2, mf_zero2); C = MF_C, bcs = mf_zerobcs(mesh_s))
    As = build_matrix(op_s)
    @test maximum(abs, Ax .- As) ≤ 1e-13 * maximum(abs, As)    # calibrated 0.0
    # cross-block SIPG symmetry (real Poisson): the pairwise (σ, Am) contract holds
    # across the connection (calibrated asym 2.5e-16, min sym eig 0.263)
    op_xp = MeshOperator(ComplexF64, mesh_x, (mf_one2, mf_one2), mf_zero2,
                         (mf_zero2, mf_zero2); C = MF_C, bcs = mf_zerobcs(mesh_x))
    Axp = build_matrix(op_xp)
    @test maximum(abs, imag.(Axp)) == 0
    Ar = real.(Axp)
    @test opnorm(Ar - transpose(Ar), Inf) / opnorm(Ar, Inf) < 1e-12
    @test minimum(eigvals(Symmetric((Ar + transpose(Ar)) / 2))) > 0.05

    # cross-block transverse 2:1 violation is a construction error
    bad = [elL,
           ElementId(96, SegmentId(0, 0), SegmentId(2, 0)),
           ElementId(96, SegmentId(0, 0), SegmentId(2, 1)),
           ElementId(96, SegmentId(0, 0), SegmentId(2, 2)),
           ElementId(96, SegmentId(0, 0), SegmentId(2, 3))]
    @test_throws ArgumentError Mesh2D([Lb, Rb], bad, (5, 5);
                                      connections = [BlockConnection(95, 1, 96)])
    # mismatched transverse maps across a connection are a construction error
    Rb2 = Block(88, LinearMap1D(0.5, 1.0), LinearMap1D(0.0, 2.0))
    @test_throws ArgumentError Mesh2D([Lb, Rb2],
        [elL, ElementId(88, SegmentId(0, 0), SegmentId(0, 0))], (5, 5);
        connections = [BlockConnection(95, 1, 88)])
    # endpoint discontinuity across a connection is a construction error
    Rb3 = Block(87, LinearMap1D(0.6, 1.0), LinearMap1D(0.0, 1.0))
    @test_throws ArgumentError Mesh2D([Lb, Rb3],
        [elL, ElementId(87, SegmentId(0, 0), SegmentId(0, 0))], (5, 5);
        connections = [BlockConnection(95, 1, 87)])
end

@testset "BigFloat(256) pass-through" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        blk = Block(98, LinearMap1D(big"0.0", big"1.0"), LinearMap1D(big"0.0", big"1.0"))
        e1 = ElementId(98, SegmentId(1, 0), SegmentId(0, 0))
        e2 = ElementId(98, SegmentId(1, 1), SegmentId(0, 0))
        mesh = Mesh2D(blk, [e1, e2], (4, 4))
        oneB(x, y) = big"1.0"
        zeroB(x, y) = big"0.0"
        bcs = Dict{Any,Any}()
        for (el, d, sidx) in boundary_faces(mesh)
            bcs[(el, d, sidx)] = DirichletGhost(zeros(BigFloat, 4))
        end
        op = MeshOperator(BigFloat, mesh, (oneB, oneB), zeroB, (zeroB, zeroB);
                          C = MF_C, bcs = bcs)
        # exactly-represented zero-boundary polynomial: R·u == M·S at BigFloat precision
        up(x, y) = x * (1 - x) * y * (1 - y)
        Sp(x, y) = 2 * x * (1 - x) + 2 * y * (1 - y)
        u = zeros(BigFloat, mesh.ndofs)
        b = zeros(BigFloat, mesh.ndofs)
        for (k, el) in enumerate([e1, e2])
            em = element_matrices(BigFloat, blk, el, (4, 4))
            o = mesh.offsets[k]
            u[o+1:o+16] .= vec([up(em.x[1][i], em.x[2][j]) for i in 1:4, j in 1:4])
            b[o+1:o+16] .= em.M * vec([Sp(em.x[1][i], em.x[2][j]) for i in 1:4, j in 1:4])
        end
        y = apply!(zeros(BigFloat, mesh.ndofs), op, u)
        @test eltype(y) == BigFloat && precision(real(eltype(y))) == 256
        @test maximum(abs, y .- b) ≤ big"1e-60" * maximum(abs, b)   # calibrated 1.0e-75
        # generic direct-solve bypass (dense generic LU) preserves precision
        x = direct_solve(op, b)
        @test eltype(x) == BigFloat
        @test maximum(abs, x .- u) ≤ big"1e-60"                     # calibrated 2.1e-77
    finally
        setprecision(BigFloat, old)
    end
end

end # testset matrixfree
