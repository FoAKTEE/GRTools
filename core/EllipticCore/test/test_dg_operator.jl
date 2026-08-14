# test_dg_operator.jl — grounded gates for src/dg_operator.jl (DAG node: dg-strong-residuals,
# elliptic-ssf GLOBAL_DAG row 11 + TOPOLOGICAL_SCHEDULE W4-A, CRITICAL PATH). Written BEFORE
# the implementation (TDD red run captured in evidence/dg-strong-residuals/pkg_test_red.log).
#
# Target: the elementwise strong-form DG residual of CRP eq:dgres_strong (main.tex:565-574)
# in its matrix form eq:dgres_strong_matrix (main.tex:636-654), over the sealed wave-3 matrix
# layer (dg_matrices.jl):
#
#   R = − MD_i·F^i + M·(β Ψ + γ^i (DΨ)_i)
#       + MLD^T_i · (𝔄^{ij} n_j (Ψ* − Ψ))  −  ML · ((n_i F^i)* − n_i F^i),
#
# with F^i = 𝔄^{ii} ∂_i Ψ (𝔄 diagonal, main.tex:575), coefficients 𝔄/β/γ CALLABLES of the
# physical coordinates evaluated per element (volume nodes + face points), and the numerical
# flux (Ψ*, (nF)*) pluggable — only the central/average instance (Ψ* = {{Ψ}}, (nF)* = {{nF}},
# the σ→0 limit of eq:numflux main.tex:584-593) ships with this node.
#
# BINDING sign/normal conventions (owned by this node, asserted below):
#   * tensor-product maps ⇒ faces are axis-aligned; outward normal of face (d, s) is s·ê_d;
#     n_i F^i = s F^d and 𝔄^{ij} n_j = s 𝔄^{dd} on that face.
#   * FaceState carries (ψ, φ = n·F) w.r.t. the OWNER's outward normal, at the face's
#     transverse collocation points (wave-3 face-matrix column ordering).
#   * exchange is SENDER-convention: ext.nflux is w.r.t. the NEIGHBOR's outward normal
#     (= minus the receiver's), so CentralFlux gives Ψ* = (ψint+ψext)/2 and, w.r.t. the
#     receiver's normal, (nF)* = (φint − φext)/2.
#   * interface plumbing is a pure function FaceState → FaceState (`exterior_state` = the
#     identity default); jump/BC hooks later substitute their own (spectre
#     GetModifyBoundaryData pattern, audit-spectre-ssf §2).
#   * transmissive/no-op exterior = FaceState(ψint, −φint): BOTH face corrections vanish
#     identically (gated bitwise-level below).
#
# Oracle gates (worker node contract, constants calibrated PRE-RED in scratchpad
# calibrate_dg_operator.jl from the sealed wave-2/3 layers only):
#   1. MMS Poisson (𝔄=I, β=γ=0) on a conforming 2×2-element affine mesh: residual of the
#      exact interpolant vs M·S converges spectrally over N=4..10 — every ratio ≤ 0.5
#      (measured ≤ 0.19), LSQ log10 slope ≤ −0.7/order (measured −1.01), err(10) ≤ 1e-6
#      (measured 4.4e-8; roundoff floor not reached).
#   2. discrete divergence theorem through the residual layer's OWN face states, random
#      nodal fields and smooth non-polynomial 𝔄 callables, affine AND compaction elements:
#      Σ 1ᵀMD_d F^d == Σ_faces 1ᵀ ML φ to 1e-12 relative (algebraic identity: MD annihilates
#      constants + wave-3 SBP; measured ~8e-16); sign falsifier: unsigned φ breaks it > 1e-3.
#   3. complex smoke: complex β/γ/𝔄 ⇒ complex residual with nonzero imaginary part;
#      conjugation consistency R(conj Ψ; conj coeffs, conj ext) == conj R(Ψ) to 1e-15·scale
#      (real matrices commute with conj); global-mesh complex linearity to 1e-13.
#   4. compaction-mapped element consistency via the wave-3 SBP identity:
#      −MD_d F == MD_dᵀ F − Σ_s ML_s φ_s per dimension to 1e-12 (measured ~8e-16), and
#      mirror-exterior residual == volume-only reference to 1e-13.
#   5. face-normal antisymmetry between conforming neighbors: ψ*_L == ψ*_R and
#      φ*_L == −φ*_R; ML-weighted interface contributions telescope in the global sum to
#      1e-12 (measured ~1e-15); for globally-continuous polynomial flux both face
#      corrections vanish elementwise to 1e-12.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test files,
# ---- distinct name/seed; keeps stdlib Random OUT of the test env) --------------------------
mutable struct DGORNG
    s::UInt64
end
dgo_unit!(r::DGORNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                        Int64(r.s >> 11) / 2.0^53 * 2 - 1)
dgo_vec!(r::DGORNG, n) = [dgo_unit!(r) for _ in 1:n]
dgo_cvec!(r::DGORNG, n) = [complex(dgo_unit!(r), dgo_unit!(r)) for _ in 1:n]

const DGO_SIDES = (-1, 1)

dgo_lin(i, j, Nx) = i + (j - 1) * Nx
dgo_maxabs(v) = maximum(abs, v)

# fixtures (same geometry family as test_dg_matrices.jl, distinct block ids)
dgo_affine_blk() = Block(11, LinearMap1D(0.0, 2.0), LinearMap1D(0.0, 1.0))
dgo_affine_ids() = [ElementId(11, SegmentId(1, ix), SegmentId(1, iy)) for ix in 0:1, iy in 0:1]
dgo_compaction_blk() = Block(13, InverseRadialMap1D(25.3, 1.0e4), LinearMap1D(0.0, 0.25))
dgo_compaction_id() = ElementId(13, SegmentId(2, 3), SegmentId(1, 1))

# independent face physical coordinates of face (d, s): normal endpoint + transverse nodes
function dgo_face_coords(blk, id, em, d, s)
    RT = eltype(em.jac[1])
    lo, hi = logical_interval(RT, id, d)
    xf = map_forward(blockmap(blk, d), s == -1 ? lo : hi)
    xt = em.x[3 - d]
    return d == 1 ? [(xf, t) for t in xt] : [(t, xf) for t in xt]
end

# independent physical-derivative reference through the sealed wave-3 identity Dphys = M⁻¹MD
dgo_dphys_ref(em, Ψ, d) = em.M \ (em.MD[d] * Ψ)

# mirror (transmissive) exterior of an interior face state
dgo_mirror(fs) = FaceState(fs.psi, -fs.nflux)

# assemble per-element states over a conforming one-block mesh; returns Dicts keyed by id
function dgo_mesh_states(blk, ids, N, A, β, γ, Ψs)
    ems = Dict(id => element_matrices(Float64, blk, id, N) for id in ids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, A, β, γ) for id in ids)
    sts = Dict(id => element_state(ems[id], cs[id], Ψs[id]) for id in ids)
    return ems, cs, sts
end

# residual of every element: interior ext states via exterior_state ∘ neighbor face state
# (identity plumbing), block-boundary ext states from the `bext(id, d, sidx, em)` callback
function dgo_mesh_residuals(blk, ids, ems, cs, sts, bext)
    Rs = Dict{ElementId,Vector}()
    for id in ids
        ext = [[begin
                    s = DGO_SIDES[sidx]
                    nb = neighbor_id(id, d, s)
                    nb === nothing ? bext(id, d, sidx, ems[id]) :
                        exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1])
                end for sidx in 1:2] for d in 1:2]
        Rs[id] = strong_residual(ems[id], cs[id], sts[id], ext)
    end
    return Rs
end

# volume-only reference residual (no face terms), from wave-3 objects + state fields
dgo_volume_ref(em, c, st) =
    em.M * (c.beta .* st.psi .+ c.gamma[1] .* st.dpsi[1] .+ c.gamma[2] .* st.dpsi[2]) .-
    em.MD[1] * st.F[1] .- em.MD[2] * st.F[2]

@testset "dg_operator (dg-strong-residuals)" begin

# constant / zero / one closures used throughout
one2(x, y) = 1.0
zero2(x, y) = 0.0

@testset "API basics: coefficients, state, face conventions, validation" begin
    rng = DGORNG(0xd60_2026_0728_0004)
    blk = dgo_affine_blk()
    id = dgo_affine_ids()[2, 1]                      # (ix, iy) = (1, 0)
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))

    A1(x, y) = 1.3 + 0.21x - 0.11y
    A2(x, y) = 0.8 + 0.07x + 0.19y
    B(x, y) = 0.4 - 0.3x + 0.2y
    G1(x, y) = 0.9 + 0.1x
    G2(x, y) = -0.6 + 0.05y
    c = element_coefficients(em, blk, id, (A1, A2), B, (G1, G2))

    # volume values at the tensor nodes, i-fastest ordering
    for j in 1:Ny, i in 1:Nx
        p = dgo_lin(i, j, Nx)
        x, y = em.x[1][i], em.x[2][j]
        @test c.A[1][p] ≈ A1(x, y) rtol = 1e-14
        @test c.A[2][p] ≈ A2(x, y) rtol = 1e-14
        @test c.beta[p] ≈ B(x, y) rtol = 1e-14
        @test c.gamma[1][p] ≈ G1(x, y) rtol = 1e-14
        @test c.gamma[2][p] ≈ G2(x, y) rtol = 1e-14
    end
    # face coefficient values AT the face points (normal endpoint × transverse nodes)
    for d in 1:2, sidx in 1:2
        pts = dgo_face_coords(blk, id, em, d, DGO_SIDES[sidx])
        Af = d == 1 ? A1 : A2
        @test c.Aface[d][sidx] ≈ [Af(p...) for p in pts] rtol = 1e-14
    end

    # element state: derivatives, fluxes, face states with the OWN-normal sign convention
    Ψ = dgo_vec!(rng, Nx * Ny)
    st = element_state(em, c, Ψ)
    @test st.psi == Ψ
    for d in 1:2
        @test dgo_maxabs(st.dpsi[d] - dgo_dphys_ref(em, Ψ, d)) <=
              1e-12 * max(1.0, dgo_maxabs(st.dpsi[d]))
        @test st.F[d] ≈ c.A[d] .* st.dpsi[d] rtol = 1e-14
        for sidx in 1:2
            s = DGO_SIDES[sidx]
            fm = face_matrices(em, d, s)
            fs = st.faces[d][sidx]
            @test fs.psi ≈ transpose(fm.I) * Ψ rtol = 1e-13
            @test fs.nflux ≈ s .* (transpose(fm.I) * st.F[d]) rtol = 1e-13
        end
    end

    # central flux formulas + exterior_state identity plumbing
    nf = length(st.faces[1][2].psi)
    e = FaceState(dgo_vec!(rng, nf), dgo_vec!(rng, nf))
    fs = st.faces[1][2]
    ψstar, φstar = numerical_flux(CentralFlux(), face_matrices(em, 1, +1),
                                  c.Aface[1][2], fs, e)
    @test ψstar ≈ (fs.psi .+ e.psi) ./ 2 rtol = 1e-15
    @test φstar ≈ (fs.nflux .- e.nflux) ./ 2 rtol = 1e-15
    @test exterior_state(e) === e
    @test CentralFlux() isa NumericalFlux

    # transmissive (mirror) exterior kills BOTH face corrections: residual == volume-only
    ext = [[dgo_mirror(st.faces[d][sidx]) for sidx in 1:2] for d in 1:2]
    R = strong_residual(em, c, st, ext)
    Rvol = dgo_volume_ref(em, c, st)
    @test dgo_maxabs(R - Rvol) <= 1e-13 * max(1.0, dgo_maxabs(Rvol))
    # in-place variant agrees
    R2 = zeros(Float64, Nx * Ny)
    strong_residual!(R2, em, c, st, ext)
    @test R2 == R
    # β=1, γ=0, constant Ψ ⇒ F ~ 0 and R == M·Ψ
    c1 = element_coefficients(em, blk, id, (A1, A2), one2, (zero2, zero2))
    st1 = element_state(em, c1, fill(0.7, Nx * Ny))
    ext1 = [[dgo_mirror(st1.faces[d][sidx]) for sidx in 1:2] for d in 1:2]
    R1 = strong_residual(em, c1, st1, ext1)
    @test R1 ≈ em.M * fill(0.7, Nx * Ny) rtol = 1e-12

    # BigFloat(256) smoke: precision inherited end-to-end through state + residual
    setprecision(BigFloat, 256) do
        blkb = Block(11, LinearMap1D(big"0.0", big"2.0"), LinearMap1D(big"0.0", big"1.0"))
        emb = element_matrices(BigFloat, blkb, id, (3, 3))
        cb = element_coefficients(emb, blkb, id, ((x, y) -> big"1.0", (x, y) -> big"1.0"),
                                  (x, y) -> big"0.0", ((x, y) -> big"0.0", (x, y) -> big"0.0"))
        stb = element_state(emb, cb, [big(i)^2 / 10 for i in 1:9])
        @test eltype(stb.psi) === BigFloat && precision(stb.dpsi[1][1]) == 256
        extb = [[dgo_mirror(stb.faces[d][sidx]) for sidx in 1:2] for d in 1:2]
        Rb = strong_residual(emb, cb, stb, extb)
        @test eltype(Rb) === BigFloat && precision(Rb[1]) == 256
    end

    # validation
    @test_throws ArgumentError element_coefficients(em, blk,
        ElementId(12, SegmentId(0, 0), SegmentId(0, 0)), (A1, A2), B, (G1, G2))
    @test_throws ArgumentError element_state(em, c, ones(Nx * Ny + 1))
    @test_throws ArgumentError strong_residual!(zeros(3), em, c, st, ext)
    @test_throws ArgumentError FaceState([1.0, 2.0], [1.0])
end

@testset "MMS: Poisson spectral p-convergence, conforming 2×2 affine mesh, N=4..10" begin
    a1, b1, c1, d1 = 1.7, 0.3, 1.3, -0.4
    u(x, y) = sin(a1 * x + b1) * cos(c1 * y + d1)
    ux(x, y) = a1 * cos(a1 * x + b1) * cos(c1 * y + d1)
    uy(x, y) = -c1 * sin(a1 * x + b1) * sin(c1 * y + d1)
    S(x, y) = (a1^2 + c1^2) * u(x, y)                 # S = −Δu (eq:first-order with 𝔄=I)

    blk = dgo_affine_blk()
    ids = dgo_affine_ids()
    errs = Float64[]
    for N in 4:10
        Ψs = Dict{ElementId,Vector{Float64}}()
        emtmp = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in ids)
        for id in ids
            em = emtmp[id]
            Ψs[id] = vec([u(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
        end
        ems, cs, sts = dgo_mesh_states(blk, ids, (N, N), (one2, one2), zero2,
                                       (zero2, zero2), Ψs)
        # exterior ghost of the EXACT solution at block boundaries (sender convention:
        # ghost outward normal = −receiver normal ⇒ φ_ghost = −s·∂_d u)
        bext = (id, d, sidx, em) -> begin
            s = DGO_SIDES[sidx]
            pts = dgo_face_coords(blk, id, em, d, s)
            g = d == 1 ? ux : uy
            FaceState([u(p...) for p in pts], [-s * g(p...) for p in pts])
        end
        Rs = dgo_mesh_residuals(blk, ids, ems, cs, sts, bext)
        emax, smax = 0.0, 0.0
        for id in ids
            em = ems[id]
            MS = em.M * vec([S(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
            emax = max(emax, dgo_maxabs(Rs[id] - MS))
            smax = max(smax, dgo_maxabs(MS))
        end
        push!(errs, emax / smax)
    end
    # spectral decay gates (constants calibrated pre-red; measured ratios ≤ 0.19,
    # slope −1.01/order, err(10) = 4.4e-8 — margins ×2.6, ×1.4, ×20)
    @test all(errs[2:end] ./ errs[1:end-1] .<= 0.5)
    Nv = collect(4:10)
    le = log10.(errs)
    slope = (length(Nv) * sum(Nv .* le) - sum(Nv) * sum(le)) /
            (length(Nv) * sum(Nv .^ 2) - sum(Nv)^2)
    @test slope <= -0.7
    @test errs[end] <= 1e-6
end

divg_fixtures = (("affine", dgo_affine_blk(), dgo_affine_ids()[1, 2], (5, 4),
                  0xd60_2026_0728_0011),
                 ("compaction", dgo_compaction_blk(), dgo_compaction_id(), (6, 5),
                  0xd60_2026_0728_0012))
@testset "divergence theorem via residual face states [$name]" for (name, blk, id, N, seed) in divg_fixtures
    rng = DGORNG(seed)
    em = element_matrices(Float64, blk, id, N)
    # smooth NON-polynomial coefficients: the identity is algebraic (holds for any nodal F)
    A1(x, y) = 0.8 + 0.3 * sin(x / 500) * cos(y)
    A2(x, y) = 1.1 - 0.2 * cos(x / 900) * sin(y)
    c = element_coefficients(em, blk, id, (A1, A2), zero2, (zero2, zero2))
    Ψ = dgo_vec!(rng, N[1] * N[2])
    st = element_state(em, c, Ψ)
    lhs = sum(em.MD[1] * st.F[1]) + sum(em.MD[2] * st.F[2])
    rhs, wrong, scale = 0.0, 0.0, abs(lhs)
    for d in 1:2, sidx in 1:2
        s = DGO_SIDES[sidx]
        fm = face_matrices(em, d, s)
        t = sum(fm.ML * st.faces[d][sidx].nflux)
        rhs += t
        wrong += s * t                                # unsigned-φ falsifier (drops n_d)
        scale = max(scale, abs(t))
    end
    @test abs(lhs - rhs) <= 1e-12 * scale
    @test abs(lhs - wrong) > 1e-3 * scale
end

@testset "complex coefficients: smoke, conjugation consistency, mesh linearity" begin
    rng = DGORNG(0xd60_2026_0728_0021)
    blk = dgo_affine_blk()
    id = dgo_affine_ids()[1, 1]
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))

    A1(x, y) = 1.0 + 0.1x
    A2(x, y) = 0.9 + 0.05im * y
    B(x, y) = 0.3 - 0.2x + (1.1 + 0.4y) * im
    G1(x, y) = 0.1 + 0.25im * x
    G2(x, y) = -0.2 + 0.15im * y
    c = element_coefficients(em, blk, id, (A1, A2), B, (G1, G2))
    @test eltype(c.beta) <: Complex

    # complex β/γ/𝔄 ⇒ complex residual with nonzero imaginary part (REAL Ψ, mirror ext:
    # the imaginary part comes from the coefficients alone)
    Ψ = dgo_vec!(rng, Nx * Ny)
    st = element_state(em, c, Ψ)
    extm = [[dgo_mirror(st.faces[d][sidx]) for sidx in 1:2] for d in 1:2]
    R0 = strong_residual(em, c, st, extm)
    @test eltype(R0) <: Complex
    @test dgo_maxabs(imag.(R0)) > 1e-3 * dgo_maxabs(R0)

    # conjugation consistency on fully complex data: R(conj Ψ; conj coeffs, conj ext)
    # == conj R(Ψ)  (all matrices are REAL — conjugation commutes with the assembly)
    nf1, nf2 = Ny, Nx                        # face point counts: dim-1 face ↔ Ny, dim-2 ↔ Nx
    Ψc = dgo_cvec!(rng, Nx * Ny)
    stc0 = element_state(em, c, Ψc)
    ext = [[FaceState(dgo_cvec!(rng, d == 1 ? nf1 : nf2), dgo_cvec!(rng, d == 1 ? nf1 : nf2))
            for _ in 1:2] for d in 1:2]
    R = strong_residual(em, c, stc0, ext)
    cc = element_coefficients(em, blk, id,
                              ((x, y) -> conj(A1(x, y)), (x, y) -> conj(A2(x, y))),
                              (x, y) -> conj(B(x, y)),
                              ((x, y) -> conj(G1(x, y)), (x, y) -> conj(G2(x, y))))
    stc = element_state(em, cc, conj.(Ψc))
    extc = [[FaceState(conj.(ext[d][sidx].psi), conj.(ext[d][sidx].nflux))
             for sidx in 1:2] for d in 1:2]
    Rc = strong_residual(em, cc, stc, extc)
    @test dgo_maxabs(Rc - conj.(R)) <= 1e-15 * max(1.0, dgo_maxabs(R))

    # global 2×2-mesh complex linearity (zero boundary exterior, identity plumbing)
    ids = dgo_affine_ids()
    α, μ = 0.3 - 1.2im, -0.7 + 0.4im
    us = Dict(idk => dgo_cvec!(rng, Nx * Ny) for idk in ids)
    vs = Dict(idk => dgo_cvec!(rng, Nx * Ny) for idk in ids)
    ws = Dict(idk => α .* us[idk] .+ μ .* vs[idk] for idk in ids)
    bext = (idk, d, sidx, emk) -> FaceState(zeros(ComplexF64, d == 1 ? Ny : Nx),
                                            zeros(ComplexF64, d == 1 ? Ny : Nx))
    applyR(Ψs) = begin
        ems, cs, sts = dgo_mesh_states(blk, ids, (Nx, Ny), (A1, A2), B, (G1, G2), Ψs)
        dgo_mesh_residuals(blk, ids, ems, cs, sts, bext)
    end
    Ru, Rv, Rw = applyR(us), applyR(vs), applyR(ws)
    for idk in ids
        comb = α .* Ru[idk] .+ μ .* Rv[idk]
        @test dgo_maxabs(Rw[idk] - comb) <= 1e-13 * max(1.0, dgo_maxabs(comb))
    end
end

@testset "compaction element: wave-3 SBP identity through the residual layer" begin
    rng = DGORNG(0xd60_2026_0728_0031)
    blk, id = dgo_compaction_blk(), dgo_compaction_id()
    Nx, Ny = 6, 5
    em = element_matrices(Float64, blk, id, (Nx, Ny))
    A1(x, y) = 0.8 + 0.3 * sin(x / 500) * cos(y)
    A2(x, y) = 1.1 - 0.2 * cos(x / 900) * sin(y)
    B(x, y) = 0.5 + 0.1 * y
    G1(x, y) = 0.02 - 0.01 * y
    G2(x, y) = 0.3 + 0.2 * y
    c = element_coefficients(em, blk, id, (A1, A2), B, (G1, G2))
    Ψ = dgo_vec!(rng, Nx * Ny)
    st = element_state(em, c, Ψ)

    # strong ↔ weak: −MD_d F == MD_dᵀ F − Σ_s ML_s φ_s (algebraic via the wave-3 2D SBP
    # identity MD_d + MD_dᵀ == ML⁽ᵈ'⁺⁾(I⁽ᵈ'⁺⁾)ᵀ − ML⁽ᵈ'⁻⁾(I⁽ᵈ'⁻⁾)ᵀ and φ = s·IᵀF)
    for d in 1:2
        v1 = -(em.MD[d] * st.F[d])
        v2 = transpose(em.MD[d]) * st.F[d]
        for sidx in 1:2
            v2 -= face_matrices(em, d, DGO_SIDES[sidx]).ML * st.faces[d][sidx].nflux
        end
        scale = max(dgo_maxabs(v1), dgo_maxabs(transpose(em.MD[d]) * st.F[d]))
        @test dgo_maxabs(v1 - v2) <= 1e-12 * scale
    end

    # mirror exterior ⇒ residual == volume-only reference (face corrections vanish)
    ext = [[dgo_mirror(st.faces[d][sidx]) for sidx in 1:2] for d in 1:2]
    R = strong_residual(em, c, st, ext)
    Rvol = dgo_volume_ref(em, c, st)
    @test dgo_maxabs(R - Rvol) <= 1e-13 * max(1.0, dgo_maxabs(Rvol))
end

@testset "face-normal antisymmetry + interface telescoping (central flux)" begin
    rng = DGORNG(0xd60_2026_0728_0041)
    blk = dgo_affine_blk()
    ids = dgo_affine_ids()
    Nx, Ny = 5, 5
    A1(x, y) = 0.9 + 0.13x - 0.07y
    A2(x, y) = 1.2 - 0.09x + 0.11y

    # (a) two-neighbor exchange: idL upper-x face ↔ idR lower-x face
    idL, idR = ids[1, 1], ids[2, 1]
    emL = element_matrices(Float64, blk, idL, (Nx, Ny))
    emR = element_matrices(Float64, blk, idR, (Nx, Ny))
    cL = element_coefficients(emL, blk, idL, (A1, A2), zero2, (zero2, zero2))
    cR = element_coefficients(emR, blk, idR, (A1, A2), zero2, (zero2, zero2))
    stL = element_state(emL, cL, dgo_vec!(rng, Nx * Ny))
    stR = element_state(emR, cR, dgo_vec!(rng, Nx * Ny))
    fsL, fsR = stL.faces[1][2], stR.faces[1][1]
    ψL, φL = numerical_flux(CentralFlux(), face_matrices(emL, 1, +1), cL.Aface[1][2],
                            fsL, exterior_state(fsR))
    ψR, φR = numerical_flux(CentralFlux(), face_matrices(emR, 1, -1), cR.Aface[1][1],
                            fsR, exterior_state(fsL))
    @test dgo_maxabs(ψL - ψR) <= 1e-15 * max(1.0, dgo_maxabs(ψL))
    @test dgo_maxabs(φL + φR) <= 1e-15 * max(1.0, dgo_maxabs(φL))   # (nF)*_L == −(nF)*_R
    # ML-weighted contributions cancel in the global sum (conforming face: same ML column sums)
    tL = sum(face_matrices(emL, 1, +1).ML * φL)
    tR = sum(face_matrices(emR, 1, -1).ML * φR)
    @test abs(tL + tR) <= 1e-12 * max(1.0, abs(tL))

    # (b) global conservation on the 2×2 mesh: Σ_k 1ᵀR_k == −Σ_{boundary faces} 1ᵀ ML φ
    # (interior interface lifts telescope; MLDT columns sum to zero since Σ_k ℓ′_k ≡ 0)
    Ψs = Dict(idk => dgo_vec!(rng, Nx * Ny) for idk in ids)
    ems, cs, sts = dgo_mesh_states(blk, ids, (Nx, Ny), (A1, A2), zero2, (zero2, zero2), Ψs)
    bext = (idk, d, sidx, emk) -> dgo_mirror(sts[idk].faces[d][sidx])
    Rs = dgo_mesh_residuals(blk, ids, ems, cs, sts, bext)
    tot = sum(sum(Rs[idk]) for idk in ids)
    bnd, scale = 0.0, 1.0
    for idk in ids, d in 1:2, sidx in 1:2
        if neighbor_id(idk, d, DGO_SIDES[sidx]) === nothing
            t = sum(face_matrices(ems[idk], d, DGO_SIDES[sidx]).ML *
                    sts[idk].faces[d][sidx].nflux)
            bnd += t
            scale = max(scale, abs(t))
        end
    end
    @test abs(tot + bnd) <= 1e-12 * scale

    # (c) globally-continuous polynomial flux (u global poly of per-dim degree ≤ N−1,
    # affine 𝔄 ⇒ F degree ≤ N−1 ⇒ exact matching face data): the two interface corrections
    # vanish elementwise — every element's residual equals its volume-only reference
    up(x, y) = (0.3 + 0.7x - 0.2x^2 + 0.05x^3) * (1.1 - 0.4y + 0.3y^2)
    Ψp = Dict{ElementId,Vector{Float64}}()
    for idk in ids
        em = ems[idk]
        Ψp[idk] = vec([up(em.x[1][i], em.x[2][j]) for i in 1:Nx, j in 1:Ny])
    end
    emsp, csp, stsp = dgo_mesh_states(blk, ids, (Nx, Ny), (A1, A2), zero2,
                                      (zero2, zero2), Ψp)
    bextp = (idk, d, sidx, emk) -> dgo_mirror(stsp[idk].faces[d][sidx])
    Rsp = dgo_mesh_residuals(blk, ids, emsp, csp, stsp, bextp)
    for idk in ids
        Rvol = dgo_volume_ref(emsp[idk], csp[idk], stsp[idk])
        @test dgo_maxabs(Rsp[idk] - Rvol) <= 1e-12 * max(1.0, dgo_maxabs(Rvol))
    end
end

end # testset dg_operator
