# test_boundary_hooks.jl — grounded gates for src/boundary_hooks.jl (DAG node:
# boundary-data-hooks, elliptic-ssf GLOBAL_DAG row 14 + TOPOLOGICAL_SCHEDULE W5-B). Written
# BEFORE the implementation (TDD red run captured in evidence/boundary-data-hooks/
# pkg_test_red.log).
#
# Target: the ModifyBoundaryData seam (spectre GetModifyBoundaryData pattern,
# audit-spectre-ssf §1.2/§2) — typed hooks that post-process RECEIVED exterior FaceState
# data before the numerical flux, over the N6 plumbing of dg_operator.jl (conventions
# N1-N7 inherited verbatim), in TWO modes:
#
#   (a) LINEARIZED hooks (`LinearizedHook`) — modifications LINEAR-HOMOGENEOUS in the
#       trial-field trace data (the vtu −imΩ⟨⟨Ψ⟩⟩ flux jump of spectre Equations.cpp:92-116
#       is the downstream instance, dg-jump-conditions node) that act inside EVERY operator
#       application — they are part of the GMRES/Schwarz/multigrid matrix
#       (spectre DgOperator.hpp:708-728).
#   (b) FIXED-DATA hooks (`FixedDataHook`) — additive shifts of the received data
#       INDEPENDENT of the trial field (the worldtube ±Ψ^P jumps of spectre
#       Equations.cpp:70-91 are the downstream instance), dropped from the linearized
#       application and moved to the RHS via OPERATOR-ON-ZERO: R_full(0) is exactly the
#       fixed-source contribution, b -= R_full(0) keeps the matrix linear-homogeneous
#       (spectre DgOperator.hpp:1081-1130, impose_inhomogeneous_boundary_conditions_on_source).
#
# plus the external ghost-BC seam generalized from the W4-B zero-Dirichlet ghost: the ghost
# exterior as a function of interior trace + boundary data (`ghost_state`) with its own
# linearized twin (`ghost_state_linearized`, data zeroed) so inhomogeneous BC data ride the
# same operator-on-zero path:
#   DirichletGhost(g): u_ext = 2g − u_int, F_ext = F_int (Vu:2024 eq:dirichlet) — in the
#     SENDER convention N3: FaceState(2g − ψ_int, −φ_int); g = 0 reproduces
#     zero_dirichlet_ghost exactly.
#   NeumannGhost(h): u_ext = u_int, (n·F)_ext = 2h − (n·F)_int w.r.t. the RECEIVER normal
#     (Vu:2024 eq:neumann) — sender convention: FaceState(ψ_int, φ_int − 2h); through the
#     flux, Ψ* = ψ_int exactly and (nF)* = h.
#
# Jump algebra gated here (derived from N1-N7; the physical-jump signs mirror spectre's
# worldtube hook — ψ shifts with OPPOSITE sign on the two sides, φ shifts with the SAME
# sign because n·F is w.r.t. the SENDER's normal): a prescribed interface jump
# [Ψ] = J, [∂Ψ] = K is imposed by FixedDataHooks (Δψ, Δφ) = (−J, +K) on the left receiver
# and (+J, +K) on the right receiver; the hooked exterior then equals the receiver's own
# continuous trace for the exact piecewise solution, so all face corrections vanish and
# the exact solution is an exact discrete solution.
#
# Oracle gates (worker node contract, constants calibrated PRE-RED in scratchpad
# calibrate_boundary_hooks.jl from the sealed wave-2/3/4 layers plus a throwaway
# closure-pair prototype — NOT the final API):
#   1. two-element prescribed-jump solve reproduces the exact piecewise solution:
#      L2 ≤ 1e-10 (NODE GATE; calibrated 9.5e-15), and the SOLVED field's trace jumps
#      match: |[Ψ]−J| ≤ 1e-12 (calibrated 4.4e-16), |[∂Ψ]−K| ≤ 1e-10 (calibrated 3.6e-14).
#   2. operator-on-zero fixed-source path ≡ direct inhomogeneous assembly (independent
#      per-face closed forms from the wave-3 matrices) rel ≤ 1e-13 (NODE GATE; calibrated
#      1.5e-16); affine decomposition R_full(u) = R_lin(u) + R_full(0) rel ≤ 1e-13
#      (calibrated 2.4e-16, composed variant 2.0e-16).
#   3. linearized-hook linearity: A(u+v) = A(u) + A(v) rel ≤ 1e-13 (NODE GATE; calibrated
#      2.4e-16), A(αu) = αA(u) rel ≤ 1e-13 (calibrated 3.0e-16); for a pure-linearized hook
#      set, full ≡ linearized apply and R_full(0) ≡ 0 (calibrated exactly 0.0).
#   4. hook composition (linearized + fixed simultaneously) on a manufactured complex
#      problem with junction conditions [u] = J, [u'] = K − iω⟨⟨u⟩⟩: solve L2 ≤ 1e-10
#      (calibrated 1.1e-14); FALSIFIER — dropping the linearized hook from the matrix
#      leaves an O(1) error > 1e-3 (calibrated 4.9e-2): the hook term is load-bearing in
#      the matrix.
#   5. no-hook identity: NoHook + zero-data ghosts == the wave-4 sweep (exterior_state +
#      zero_dirichlet_ghost) exactly (==; calibrated maxdiff 0.0), full and linearized.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test files,
# ---- distinct name/seed) ---------------------------------------------------------------
mutable struct BDHRNG
    s::UInt64
end
bdh_unit!(r::BDHRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                        Int64(r.s >> 11) / 2.0^53 * 2 - 1)
bdh_vec!(r::BDHRNG, n) = [bdh_unit!(r) for _ in 1:n]
bdh_cvec!(r::BDHRNG, n) = [complex(bdh_unit!(r), bdh_unit!(r)) for _ in 1:n]

const BDH_SIDES = (-1, 1)
const BDH_N = 6                                     # piecewise-cubic exact at N ≥ 4
const BDH_XM = 0.5                                  # interface abscissa of the 2×1 mesh

# ---- two-element (2 in x, 1 in y) harness: mesh-level orchestration stays deliberately
# ---- test-side here, as in test_dg_penalty.jl (connectivity containers belong to the
# ---- matrix-free-operator node) ---------------------------------------------------------
function bdh_setup(blk, N, A, β, γ)
    ids = [ElementId(blk.id, SegmentId(1, ix), SegmentId(0, 0)) for ix in 0:1]
    ems = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in ids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, A, β, γ) for id in ids)
    off = Dict(id => (k - 1) * N * N for (k, id) in enumerate(ids))
    return ids, ems, cs, off
end

# global matrix-free apply through the HOOK seams: interior interfaces route the received
# neighbor state through modify_boundary_data[(_linearized)], block-boundary faces through
# ghost_state[(_linearized)]. `hooks`/`bcs`: Dict{(id, d, sidx) => hook/ghost object}.
function bdh_make_apply(ids, ems, cs, off, N, flux, hooks, bcs, linearized)
    npe = N * N
    return function (y, u)
        mod = linearized ? modify_boundary_data_linearized : modify_boundary_data
        gho = linearized ? ghost_state_linearized : ghost_state
        sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+npe))
                   for id in ids)
        for id in ids
            ext = [[begin
                        nb = neighbor_id(id, d, BDH_SIDES[sidx])
                        int = sts[id].faces[d][sidx]
                        nb === nothing ? gho(bcs[(id, d, sidx)], int) :
                            mod(hooks[(id, d, sidx)], int,
                                exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1]))
                    end for sidx in 1:2] for d in 1:2]
            y[off[id]+1:off[id]+npe] .= strong_residual(ems[id], cs[id], sts[id], ext, flux)
        end
        return y
    end
end

# massive right-hand side b = M·S per element
function bdh_rhs(ids, ems, off, N, S)
    b = zeros(ComplexF64, length(ids) * N * N)
    for id in ids
        em = ems[id]
        b[off[id]+1:off[id]+N*N] .=
            em.M * vec([S(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
    end
    return b
end

# discrete (mass-weighted) L2 error of a global solution vector vs the exact solution
function bdh_l2err(ids, ems, off, N, u, uex)
    s = 0.0
    for id in ids
        em = ems[id]
        e = u[off[id]+1:off[id]+N*N] .-
            vec([uex(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
        s += real(dot(e, em.M * e))
    end
    return sqrt(s)
end

# boundary data of a function of (x, y) at the face points of face (d, sidx)
function bdh_facedata(em, blk, id, d, sidx, f)
    lo, hi = logical_interval(Float64, id, d)
    xf = map_forward(blockmap(blk, d), sidx == 1 ? lo : hi)
    xt = em.x[3 - d]
    return d == 1 ? [f(xf, t) for t in xt] : [f(t, xf) for t in xt]
end

@testset "boundary_hooks (boundary-data-hooks)" begin

bdh_one2(x, y) = 1.0
bdh_zero2(x, y) = 0.0
bdh_blk() = Block(41, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))

@testset "hook seam: typed hooks over received FaceState data" begin
    rng = BDHRNG(0x00000000beef0051)
    nf = 5
    int = FaceState(bdh_cvec!(rng, nf), bdh_cvec!(rng, nf))
    ext = FaceState(bdh_cvec!(rng, nf), bdh_cvec!(rng, nf))

    # NoHook: the default seam is the identity on the received data — the SAME object
    # (bitwise wave-4 behavior), in both applications
    @test NoHook() isa BoundaryDataHook
    @test modify_boundary_data(NoHook(), int, ext) === ext
    @test modify_boundary_data_linearized(NoHook(), int, ext) === ext

    # LinearizedHook: wraps a callable (int, ext) -> FaceState that must be
    # linear-homogeneous in the trace data; acts identically in BOTH applications
    # (it belongs to the matrix). vtu-style instance: φ_ext -= iω({{ψ}}).
    ω = 0.9
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    lh = LinearizedHook(fω)
    @test lh isa BoundaryDataHook
    h1 = modify_boundary_data(lh, int, ext)
    h2 = modify_boundary_data_linearized(lh, int, ext)
    @test h1.psi == ext.psi && h2.psi == ext.psi
    @test h1.nflux == ext.nflux .- im * ω .* (int.psi .+ ext.psi) ./ 2
    @test h2.nflux == h1.nflux

    # FixedDataHook: additive trial-independent shift of the received data (sender
    # convention); IDENTITY in the linearized application (it is RHS data, not matrix)
    Δψ = bdh_cvec!(rng, nf)
    Δφ = bdh_cvec!(rng, nf)
    fh = FixedDataHook(Δψ, Δφ)
    @test fh isa BoundaryDataHook
    hf = modify_boundary_data(fh, int, ext)
    @test hf.psi == ext.psi .+ Δψ
    @test hf.nflux == ext.nflux .+ Δφ
    @test modify_boundary_data_linearized(fh, int, ext) === ext
    # the delta is FaceState-shaped data: mismatched lengths must throw
    @test_throws ArgumentError FixedDataHook(bdh_cvec!(rng, nf), bdh_cvec!(rng, nf + 1))
    # type promotion: complex shift on real received data promotes cleanly
    rint = FaceState(bdh_vec!(rng, nf), bdh_vec!(rng, nf))
    rext = FaceState(bdh_vec!(rng, nf), bdh_vec!(rng, nf))
    @test eltype(modify_boundary_data(fh, rint, rext).psi) == ComplexF64

    # ComposedHook: stages apply IN ORDER, each receiving the ORIGINAL interior trace
    # (hooks modify received exterior data only, seam N6) and the previous stage's output
    ch = ComposedHook(lh, fh)
    @test ch isa BoundaryDataHook
    cf = modify_boundary_data(ch, int, ext)
    step = modify_boundary_data(lh, int, ext)
    ref = modify_boundary_data(fh, int, step)
    @test cf.psi == ref.psi && cf.nflux == ref.nflux
    # order matters when stages interact through ψ: fixed-then-linear shifts the {{ψ}}
    ch2 = ComposedHook(fh, lh)
    cf2 = modify_boundary_data(ch2, int, ext)
    @test maximum(abs, cf2.nflux .- cf.nflux) > 1e-3        # genuinely different order
    # linearized composition drops exactly the fixed stages
    cl = modify_boundary_data_linearized(ch, int, ext)
    ll = modify_boundary_data_linearized(lh, int, ext)
    @test cl.psi == ll.psi && cl.nflux == ll.nflux

    # BigFloat pass-through (precision inherited, DESIGN §4)
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        bint = FaceState(BigFloat.(bdh_vec!(rng, 3)), BigFloat.(bdh_vec!(rng, 3)))
        bext = FaceState(BigFloat.(bdh_vec!(rng, 3)), BigFloat.(bdh_vec!(rng, 3)))
        bh = FixedDataHook(BigFloat.(bdh_vec!(rng, 3)), BigFloat.(bdh_vec!(rng, 3)))
        bout = modify_boundary_data(bh, bint, bext)
        @test eltype(bout.psi) == BigFloat && precision(real(eltype(bout.psi))) == 256
    finally
        setprecision(BigFloat, old)
    end
end

@testset "ghost-BC seam: DirichletGhost/NeumannGhost generalize the W4-B ghost" begin
    rng = BDHRNG(0x00000000beef0052)
    nf = BDH_N
    ψ = bdh_cvec!(rng, nf)
    φ = bdh_cvec!(rng, nf)
    int = FaceState(ψ, φ)
    g = bdh_cvec!(rng, nf)
    h = bdh_cvec!(rng, nf)

    # DirichletGhost: u_ext = 2g − u_int, F_ext = F_int (Vu:2024 eq:dirichlet) in the N3
    # sender convention; zero data REPRODUCES the W4-B zero_dirichlet_ghost exactly
    d0 = DirichletGhost(zeros(ComplexF64, nf))
    @test d0 isa GhostBoundary
    z = zero_dirichlet_ghost(int)
    gz = ghost_state(d0, int)
    @test gz.psi == z.psi && gz.nflux == z.nflux
    dg = DirichletGhost(g)
    gd = ghost_state(dg, int)
    @test gd.psi == 2 .* g .- ψ
    @test gd.nflux == -φ
    # the LINEARIZED ghost zeroes the data: it IS the zero-Dirichlet ghost, exactly
    gl = ghost_state_linearized(dg, int)
    @test gl.psi == z.psi && gl.nflux == z.nflux

    # through the penalty flux: Ψ* == g (Dirichlet imposed on the average), {{nF}} = φ_int,
    # penalty riding on the jump 2(ψ_int − g)
    blk = bdh_blk()
    em = element_matrices(Float64, blk, ElementId(41, SegmentId(1, 0), SegmentId(0, 0)),
                          (nf, nf))
    fm = em.faces[1][1]
    Af = bdh_cvec!(rng, nf)
    σ = penalty_sigma(1.5, nf - 1, nf - 1, perpendicular_size(fm), perpendicular_size(fm))
    ψs, φs = numerical_flux(PenaltyFlux(), fm, Af, int, gd)
    sc = maximum(abs, ψ)
    @test maximum(abs, ψs .- g) ≤ 1e-13 * sc
    φref = φ .- 2σ .* (Af .* (ψ .- g))
    @test maximum(abs, φs .- φref) ≤ 1e-13 * maximum(abs, φref)

    # NeumannGhost: u_ext = u_int, (n·F)_ext = 2h − (n·F)_int w.r.t. the RECEIVER normal
    # (Vu:2024 eq:neumann) — sender convention FaceState(ψ_int, φ_int − 2h); through the
    # flux Ψ* == ψ_int EXACTLY (zero jump, penalty off) and (nF)* == h
    ng = NeumannGhost(h)
    @test ng isa GhostBoundary
    gn = ghost_state(ng, int)
    @test gn.psi == ψ
    @test gn.nflux == φ .- 2 .* h
    ψn, φn = numerical_flux(PenaltyFlux(), fm, Af, int, gn)
    @test ψn == ψ
    @test maximum(abs, φn .- h) ≤ 1e-13 * max(1.0, maximum(abs, h))
    # zero-data Neumann is trial-only data: linearized == full
    n0 = NeumannGhost(zeros(ComplexF64, nf))
    g0f = ghost_state(n0, int)
    g0l = ghost_state_linearized(n0, int)
    @test g0f.psi == g0l.psi && g0f.nflux == g0l.nflux

    # real trial data + real BC data stay real (no spurious complexification)
    rr = ghost_state(DirichletGhost(bdh_vec!(rng, 4)),
                     FaceState(bdh_vec!(rng, 4), bdh_vec!(rng, 4)))
    @test eltype(rr.psi) == Float64
end

@testset "no-hook identity: default hooks == wave-4 sweep, exact (==)" begin
    N = BDH_N
    rng = BDHRNG(0x00000000beef0053)
    A11(x, y) = 1 + 0.3x + 0.1y^2
    A22(x, y) = 1 + 0.2y
    βc = 0.3 + 0.2im
    ids, ems, cs, off = bdh_setup(bdh_blk(), N, (A11, A22), (x, y) -> βc,
                                  (bdh_zero2, bdh_zero2))
    L, R = ids
    n = 2 * N * N
    u = bdh_cvec!(rng, n)

    hooks = Dict((L, 1, 2) => NoHook(), (R, 1, 1) => NoHook())
    bcs = Dict{Any,Any}()
    for id in ids, d in 1:2, sidx in 1:2
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        nb === nothing || continue
        bcs[(id, d, sidx)] = DirichletGhost(zeros(ComplexF64, N))
    end

    # the wave-4 twin: exterior_state plumbing + zero_dirichlet_ghost, as in
    # test_dg_penalty.jl (the sweep this node's seams must reproduce bit-for-bit)
    yref = zeros(ComplexF64, n)
    sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+N*N))
               for id in ids)
    for id in ids
        ext = [[begin
                    nb = neighbor_id(id, d, BDH_SIDES[sidx])
                    nb === nothing ? zero_dirichlet_ghost(sts[id].faces[d][sidx]) :
                        exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1])
                end for sidx in 1:2] for d in 1:2]
        yref[off[id]+1:off[id]+N*N] .=
            strong_residual(ems[id], cs[id], sts[id], ext, PenaltyFlux())
    end

    yfull = zeros(ComplexF64, n)
    bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, false)(yfull, u)
    ylin = zeros(ComplexF64, n)
    bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, true)(ylin, u)
    @test yfull == yref                    # calibrated maxdiff 0.0
    @test ylin == yref                     # zero-data ghosts: linearized == full == wave-4
end

@testset "NODE GATE: two-element prescribed-jump solve, [Ψ]=J, [∂Ψ]=K" begin
    # exact piecewise-cubic model on [0,1]²: u_L = x³ − x + 1/4, u_R = u_L + J + K(x − x_m);
    # −∇²u = −6x on both sides, u independent of y. Jump hooks per the header sign table;
    # inhomogeneous DirichletGhost at x ∈ {0,1} (data u(0) = 1/4 ≠ 0, u(1) = 1/4 + J + K/2),
    # NeumannGhost(0) at y ∈ {0,1}.
    N = BDH_N
    J, K = 0.7, -1.3
    uL(x) = x^3 - x + 0.25
    uR(x) = uL(x) + J + K * (x - BDH_XM)
    uex(x, y) = x < BDH_XM ? uL(x) : uR(x)
    S(x, y) = -6x
    blk = bdh_blk()
    ids, ems, cs, off = bdh_setup(blk, N, (bdh_one2, bdh_one2), bdh_zero2,
                                  (bdh_zero2, bdh_zero2))
    L, R = ids
    n = 2 * N * N
    onesf = ones(N)

    # receiver-side hook data: ψ shifts OPPOSITE (−J left, +J right), φ shifts SAME (+K)
    # — the sender-normal flip cancels the side flip, exactly as in the spectre worldtube
    # hook (audit-spectre-ssf §1.2, Equations.cpp:70-91)
    hooks = Dict((L, 1, 2) => FixedDataHook(-J .* onesf, K .* onesf),
                 (R, 1, 1) => FixedDataHook(+J .* onesf, K .* onesf))
    bcs = Dict{Any,Any}()
    for id in ids, d in 1:2, sidx in 1:2
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        nb === nothing || continue
        bcs[(id, d, sidx)] = d == 1 ?
            DirichletGhost(ComplexF64.(bdh_facedata(ems[id], blk, id, d, sidx, uex))) :
            NeumannGhost(zeros(ComplexF64, N))
    end

    apl_full = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, false)
    apl_lin = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, true)
    b = bdh_rhs(ids, ems, off, N, S)
    R0 = zeros(ComplexF64, n)
    apl_full(R0, zeros(ComplexF64, n))                 # operator-on-zero → fixed sources
    res = gmres(apl_lin, b .- R0; abstol = 1e-13, reltol = 0.0, maxiter = 500, restart = 500)
    @test res.converged

    # NODE GATE 1e-10 (calibrated L2 9.5e-15, max 1.7e-14)
    @test bdh_l2err(ids, ems, off, N, res.x, uex) ≤ 1e-10
    uexact = vcat([vec([uex(ems[id].x[1][i], ems[id].x[2][j]) for i in 1:N, j in 1:N])
                   for id in ids]...)
    @test maximum(abs, res.x .- uexact) ≤ 1e-10

    # the SOLVED field's interface trace jumps reproduce the prescribed data:
    # [Ψ] = ψ_R − ψ_L == J (calibrated 4.4e-16); [∂Ψ] = F_R − F_L == K (calibrated 3.6e-14;
    # both fluxes re-expressed w.r.t. the LEFT receiver's +x̂ normal: F_R = −φ_R sent)
    stL = element_state(ems[L], cs[L], res.x[off[L]+1:off[L]+N*N])
    stR = element_state(ems[R], cs[R], res.x[off[R]+1:off[R]+N*N])
    @test maximum(abs, stR.faces[1][1].psi .- stL.faces[1][2].psi .- J) ≤ 1e-12
    @test maximum(abs, (-stR.faces[1][1].nflux) .- stL.faces[1][2].nflux .- K) ≤ 1e-10
end

@testset "NODE GATE: operator-on-zero ≡ direct inhomogeneous assembly" begin
    # random complex fixed data on EVERY seam (jump hooks on the interface, inhomogeneous
    # Dirichlet at x-boundaries, inhomogeneous Neumann at y-boundaries), variable real 𝔄,
    # complex β — R_full(0) must equal the per-face closed forms assembled DIRECTLY from
    # the wave-3 matrices (independent of numerical_flux/strong_residual):
    #   fixed jump Δ:   MLDT·(s𝔄(Δψ/2)) − ML·(σ𝔄Δψ − Δφ/2)
    #   Dirichlet g:    MLDT·(s𝔄g) − ML·(2σ𝔄g)
    #   Neumann h:      −ML·h
    N = BDH_N
    rng = BDHRNG(0x00000000beef0054)
    A11(x, y) = 1 + 0.3x + 0.1y^2
    A22(x, y) = 1 + 0.2y
    βc = 0.3 + 0.2im
    ids, ems, cs, off = bdh_setup(bdh_blk(), N, (A11, A22), (x, y) -> βc,
                                  (bdh_zero2, bdh_zero2))
    L, R = ids
    n = 2 * N * N

    hdata = Dict((L, 1, 2) => (bdh_cvec!(rng, N), bdh_cvec!(rng, N)),
                 (R, 1, 1) => (bdh_cvec!(rng, N), bdh_cvec!(rng, N)))
    hooks = Dict{Any,Any}(k => FixedDataHook(v...) for (k, v) in hdata)
    bcs = Dict{Any,Any}()
    bdata = Dict{Any,Any}()
    for id in ids, d in 1:2, sidx in 1:2
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        nb === nothing || continue
        v = bdh_cvec!(rng, N)
        bdata[(id, d, sidx)] = v
        bcs[(id, d, sidx)] = d == 1 ? DirichletGhost(v) : NeumannGhost(v)
    end
    apf = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, false)
    apl = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, true)
    R0 = zeros(ComplexF64, n)
    apf(R0, zeros(ComplexF64, n))

    direct = zeros(ComplexF64, n)
    for id in ids, d in 1:2, sidx in 1:2
        em = ems[id]
        fm = em.faces[d][sidx]
        s = BDH_SIDES[sidx]
        nvol, nf = size(fm.I)
        σ = penalty_sigma(1.5, nvol ÷ nf - 1, nvol ÷ nf - 1,
                          perpendicular_size(fm), perpendicular_size(fm))
        Af = cs[id].Aface[d][sidx]
        rows = off[id]+1:off[id]+N*N
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        if nb === nothing
            v = bdata[(id, d, sidx)]
            if d == 1
                direct[rows] .+= fm.MLDT * (s .* Af .* v) .- fm.ML * (2σ .* Af .* v)
            else
                direct[rows] .-= fm.ML * v
            end
        else
            Δψ, Δφ = hdata[(id, d, sidx)]
            direct[rows] .+= fm.MLDT * (s .* Af .* Δψ ./ 2) .-
                             fm.ML * (σ .* Af .* Δψ .- Δφ ./ 2)
        end
    end
    sc = maximum(abs, R0)
    @test sc > 1                                       # genuinely nonzero fixed sources
    @test maximum(abs, R0 .- direct) ≤ 1e-13 * sc      # NODE GATE (calibrated 1.5e-16)

    # affine decomposition: R_full(u) − R_lin(u) == R_full(0) — b −= A_hooks(0) keeps the
    # matrix linear-homogeneous (calibrated 2.4e-16)
    u = bdh_cvec!(rng, n)
    Ru = zeros(ComplexF64, n)
    apf(Ru, u)
    Lu = zeros(ComplexF64, n)
    apl(Lu, u)
    @test maximum(abs, Ru .- Lu .- R0) ≤ 1e-13 * maximum(abs, Ru)
end

@testset "NODE GATE: linearized-hook linearity in every operator application" begin
    N = BDH_N
    rng = BDHRNG(0x00000000beef0055)
    A11(x, y) = 1 + 0.3x + 0.1y^2
    A22(x, y) = 1 + 0.2y
    βc = 0.3 + 0.2im
    ids, ems, cs, off = bdh_setup(bdh_blk(), N, (A11, A22), (x, y) -> βc,
                                  (bdh_zero2, bdh_zero2))
    L, R = ids
    n = 2 * N * N
    ω = 0.9
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)
    hooks = Dict((L, 1, 2) => LinearizedHook(fω), (R, 1, 1) => LinearizedHook(fω))
    bcs = Dict{Any,Any}()
    for id in ids, d in 1:2, sidx in 1:2
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        nb === nothing || continue
        bcs[(id, d, sidx)] = d == 1 ? DirichletGhost(zeros(ComplexF64, N)) :
                                      NeumannGhost(zeros(ComplexF64, N))
    end
    apl = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, true)
    apf = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, false)

    u = bdh_cvec!(rng, n)
    v = bdh_cvec!(rng, n)
    Au = zeros(ComplexF64, n); apl(Au, u)
    Av = zeros(ComplexF64, n); apl(Av, v)
    Auv = zeros(ComplexF64, n); apl(Auv, u .+ v)
    # NODE GATE A(u+v) = A(u) + A(v) at 1e-13 (calibrated 2.4e-16)
    @test maximum(abs, Auv .- Au .- Av) ≤ 1e-13 * maximum(abs, Auv)
    # complex-scalar homogeneity (calibrated 3.0e-16)
    α = 0.37 - 1.21im
    Aαu = zeros(ComplexF64, n); apl(Aαu, α .* u)
    @test maximum(abs, Aαu .- α .* Au) ≤ 1e-13 * maximum(abs, Aαu)
    # a pure-linearized hook set acts in EVERY application: full == linearized exactly,
    # and the operator stays linear-HOMOGENEOUS: R_full(0) == 0 (calibrated exactly 0.0)
    Fu = zeros(ComplexF64, n); apf(Fu, u)
    @test Fu == Au
    F0 = zeros(ComplexF64, n); apf(F0, zeros(ComplexF64, n))
    @test all(iszero, F0)
end

@testset "hook composition (linearized + fixed) on a manufactured complex problem" begin
    # composed junction conditions (derivation in the file header / N1-N7): with the
    # vtu-style hook acting on the RAW received data followed by the fixed jump shifts,
    # the modified-flux system imposes [u] = J, [u'] = K − iω⟨⟨u⟩⟩ at the interface.
    # Manufactured piecewise cubic satisfying exactly those conditions; (−∇² + β)u = S.
    N = BDH_N
    rng = BDHRNG(0x00000000beef0056)
    ω = 0.9
    βc = 0.3 + 0.2im
    J = 0.35 - 0.15im
    K = 0.6 + 0.25im
    uL(x) = (0.8 + 0.5im) * x^3 - (0.4 - 0.2im) * x + 0.3
    K̃ = K - im * ω * (uL(BDH_XM) + J / 2)             # ⟨⟨u⟩⟩ = u_L(x_m) + J/2 exactly
    uR(x) = uL(x) + J + K̃ * (x - BDH_XM)
    uex(x, y) = x < BDH_XM ? uL(x) : uR(x)
    S(x, y) = -6 * (0.8 + 0.5im) * x + βc * uex(x, y)
    blk = bdh_blk()
    ids, ems, cs, off = bdh_setup(blk, N, (bdh_one2, bdh_one2), (x, y) -> βc,
                                  (bdh_zero2, bdh_zero2))
    L, R = ids
    n = 2 * N * N
    onesf = ones(N)
    fω(i, e) = FaceState(e.psi, e.nflux .- im * ω .* (i.psi .+ e.psi) ./ 2)

    hooks = Dict((L, 1, 2) => ComposedHook(LinearizedHook(fω),
                                           FixedDataHook(-J .* onesf, K .* onesf)),
                 (R, 1, 1) => ComposedHook(LinearizedHook(fω),
                                           FixedDataHook(+J .* onesf, K .* onesf)))
    bcs = Dict{Any,Any}()
    for id in ids, d in 1:2, sidx in 1:2
        nb = neighbor_id(id, d, BDH_SIDES[sidx])
        nb === nothing || continue
        bcs[(id, d, sidx)] = d == 1 ?
            DirichletGhost(ComplexF64.(bdh_facedata(ems[id], blk, id, d, sidx, uex))) :
            NeumannGhost(zeros(ComplexF64, N))
    end
    apf = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, false)
    apl = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooks, bcs, true)
    b = bdh_rhs(ids, ems, off, N, S)
    R0 = zeros(ComplexF64, n)
    apf(R0, zeros(ComplexF64, n))
    res = gmres(apl, b .- R0; abstol = 1e-13, reltol = 0.0, maxiter = 500, restart = 500)
    @test res.converged
    # composed solve reproduces the manufactured solution (calibrated 1.1e-14)
    @test bdh_l2err(ids, ems, off, N, res.x, uex) ≤ 1e-10

    # affine decomposition holds for COMPOSED hooks too (calibrated 2.0e-16)
    u = bdh_cvec!(rng, n)
    Ru = zeros(ComplexF64, n); apf(Ru, u)
    Lu = zeros(ComplexF64, n); apl(Lu, u)
    @test maximum(abs, Ru .- Lu .- R0) ≤ 1e-13 * maximum(abs, Ru)

    # FALSIFIER: drop the linearized stage from the hooks — R_full(0) is UNCHANGED (the
    # linear stage vanishes on zero data), but the matrix loses the −iω⟨⟨ψ⟩⟩ term and the
    # solve misses the manufactured solution at O(1) (calibrated 4.9e-2 ≫ gate 1e-3):
    # the linearized hook is load-bearing INSIDE the matrix.
    hooksx = Dict((L, 1, 2) => FixedDataHook(-J .* onesf, K .* onesf),
                  (R, 1, 1) => FixedDataHook(+J .* onesf, K .* onesf))
    apxf = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooksx, bcs, false)
    apxl = bdh_make_apply(ids, ems, cs, off, N, PenaltyFlux(), hooksx, bcs, true)
    R0x = zeros(ComplexF64, n)
    apxf(R0x, zeros(ComplexF64, n))
    @test maximum(abs, R0x .- R0) ≤ 1e-13 * maximum(abs, R0)
    resx = gmres(apxl, b .- R0x; abstol = 1e-13, reltol = 0.0, maxiter = 500, restart = 500)
    @test resx.converged
    @test bdh_l2err(ids, ems, off, N, resx.x, uex) > 1e-3
end

end # testset boundary_hooks
