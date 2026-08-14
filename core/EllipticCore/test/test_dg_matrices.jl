# test_dg_matrices.jl — grounded gates for src/dg_matrices.jl (DAG node: dg-operator-matrices,
# elliptic-ssf GLOBAL_DAG row 10 + TOPOLOGICAL_SCHEDULE W3-A, CRITICAL PATH). Written BEFORE
# the implementation (TDD red run captured in evidence/dg-operator-matrices/pkg_test_red.log).
#
# Target: the per-element DG matrix layer of CRP main.tex:647-664 (eq:dgres_strong_matrix
# constituents), specialized to the CRP's tensor-product block maps (one CoordinateMap1D per
# dimension, mesh2d.jl), whose Jacobian matrix is DIAGONAL: invjac^î_i = δ^î_i / j_i and
# J = j_1 j_2 with j_d the composed (block map ∘ dyadic segment) 1D Jacobian.
#
#   M_pq      = δ_pq J_p w_{p1} w_{p2}                               [main.tex:649]
#   MD_{i,pq} = D_{î,pq} J_q (invjac^î_i)_q w_{p1} w_{p2}            [main.tex:650]
#   ML_{p q̄}  = I_{p q̄} Ĵ_{q̄} w_{q̄}                                  [main.tex:651]
#   MLDT_{i,p q̄} = Dᵀ_{î,pq} I_{q q̄} Ĵ_{q̄} (invjac^î_i)_{q̄} w_{q̄}    [main.tex:652]
#   D_{î,pq}  = ℓ′_{q_î}(ξ_{p_î}) ∏_{ĵ≠î} δ_{p_ĵ q_ĵ}                [eq:diffmat, main.tex:657-659]
#   I_{p q̄}   = ℓ_{p_î}(±1) ∏_{ĵ≠î} δ_{p_ĵ q̄_ĵ}                      [main.tex:661-664]
#
# with Ĵ (surface Jacobian) = the transverse-dimension Jacobian at the face collocation
# point, and the face-normal invjac factor of MLDT evaluated AT the face (ξ_î = ±1), not at
# volume nodes. Only the face-NORMAL derivative direction of MLDT is delivered: the flux
# matrix 𝔄^{ij} of eq:dgres_strong is diagonal (main.tex:575) and the term carries n_j, so
# for tensor-product maps every other component is multiplied by zero.
#
# Oracle gates (GLOBAL_DAG row 10 + worker node contract):
#   * dense fine-quadrature references on a random AFFINE element: every matrix entry of
#     M, MD_1, MD_2 and all four faces' I, ML, MLDT vs braindead-loop quadrature assembly of
#     the defining inner products (eq:dgproj) — all integrands polynomial, hence the fine rule
#     is exact — to 1e-13 relative. Includes the mass-lumping THEOREM: the exact integral
#     ∫ψ_p ψ_q J dξ² is itself diagonal for affine maps (deg 2N−2 ≤ 2N−1 Gauss exactness).
#   * compaction-mapped element (paper contrast r_u=25.3, R=1e4, eq:compaction): M vs the
#     pointwise definition with independently recomputed geometry; MD_1 vs an EXACT analytic
#     fine-quadrature reference (for i = radial, ψ_p ∂_x ψ_q J = ℓ_p1 ℓ′_q1 · ℓ_p2 ℓ_q2 j_2 is
#     polynomial — the 1/j_1 chain rule cancels J's non-polynomial factor); MD_2 and dim-2
#     faces vs definitional references with independent geometry (the collocation rule IS the
#     definition there, main.tex:636: "evaluate the integrals ... with Gauss quadrature on the
#     Legendre-Gauss collocation points").
#   * 2D Jacobian-weighted generalized SBP (the precise statement, extending the wave-2 1D
#     identity W·D + Dᵀ·W == face_hi·face_hiᵀ − face_lo·face_loᵀ): for each direction i,
#         MD_i + MD_iᵀ == ML⁽ⁱ,⁺⁾·(I⁽ⁱ,⁺⁾)ᵀ − ML⁽ⁱ,⁻⁾·(I⁽ⁱ,⁻⁾)ᵀ         to 1e-13 relative,
#     on affine AND compaction elements. Proof: with tensor-product maps both sides factor per
#     transverse index t into [W D + Dᵀ W − B] · j_t(t) w_t(t); the bracket is the wave-2
#     identity (nodal vectors ↔ degree ≤ N−1 polynomials, u′v+uv′ has degree ≤ 2N−3 which the
#     rule integrates exactly, boundary values interpolated). The volume J and the 1/j_i
#     inverse-Jacobian combine into exactly the surface Jacobian — that is the 2D content.
#     Falsifier: swapping the two face lifts (wrong outward orientation) must break it badly.
#   * consistency MD_i == M · Dphys_i with Dphys_i = diag(1/j_i) (D_î ⊗ δ) (the GLOBAL_DAG
#     row-10 reading "M·D + (M·D)ᵀ = boundary lift"), and MD_i · 1 == 0 (D annihilates
#     constants ⇒ discrete divergence of a constant flux vanishes identically).
#   * polynomial exactness: for polynomial fluxes F^i of per-dim degree ≤ N−1 (physical coords
#     on the affine element; logical coords on the compaction element, radial direction), the
#     projected divergence (MD·F)_p equals the analytic ∫ψ_p ∂_iF^i J dξ² to 1e-13, and its
#     grid sum telescopes to the exact surface integral ∮ F·n Ĵ (discrete divergence theorem).
#   * BigFloat(256) pass-through on one small case: all of the above scale down to ≤ 1e-70
#     relative (Float64 arithmetic cannot pass), precision(entries) == 256.
#   * cache: element assembly draws its 1D constituents (D, weights, face vectors) from the
#     (RT, N, precision)-keyed refelement cache — exact hit/miss counts over repeated builds;
#     a ComplexF64 request normalizes to the Float64 entry and yields identical real matrices.
#
# Test helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore
using NumericsCore: gl_nodes_weights

# ---- tiny deterministic generator (keeps stdlib Random OUT of the test env; same Knuth-MMIX
# ---- LCG convention as test_refelement.jl / test_mesh2d.jl, distinct name/seed) -------------
mutable struct DGMRNG
    s::UInt64
end
dgm_unit!(r::DGMRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                        Int64(r.s >> 11) / 2.0^53 * 2 - 1)

# local linear index (independent of the implementation's tensor_lin): i fastest, column-major
lin(i, j, Nx) = i + (j - 1) * Nx

# relative max-abs gate against a reference matrix/vector
relerr(A, R) = maximum(abs, A - R) / max(one(eltype(R)), maximum(abs, R))

# ---- independently composed per-dimension element geometry ---------------------------------
# block map ∘ dyadic segment map, from the SEALED wave-2 layers only (mesh2d + refelement):
# returns closures (phys, jac) of the logical coordinate ξ ∈ [-1,1] of the ELEMENT.
function el_geometry(RT, m, id, dim)
    lo, hi = logical_interval(RT, id, dim)
    half = (hi - lo) / 2
    ζ = ξ -> (lo * (1 - ξ) + hi * (1 + ξ)) / 2
    phys = ξ -> map_forward(m, ζ(ξ))
    jac = ξ -> map_jacobian(m, ζ(ξ)) * half
    return phys, jac
end

# basis-evaluation tables on a fine Gauss rule: L[k, q] = ℓ_q(ξf_k), dL[k, q] = ℓ′_q(ξf_k)
# (columns of the differentiation matrix are the exact nodal values of ℓ′_q, degree ≤ N−2,
# so interpolating them is exact — both tables are exact-polynomial objects).
function basis_tables(el, ξf)
    L = interpolation_matrix(el, ξf)
    return L, L * el.D
end

# dense fine-quadrature volume references: entries of ∫ ψ_p ψ_q J dξ² ("M") and of
# ∫ ψ_p ∂_iψ_q J dξ² ("MD_i") assembled by braindead loops; exact whenever the integrand
# is polynomial on the fine rule.
function ref_volume(RT, el1, el2, jac1, jac2, Nx, Ny; NQ = 2 * max(Nx, Ny) + 4)
    ξf, Wf = gl_nodes_weights(RT, NQ)
    L1, dL1 = basis_tables(el1, ξf)
    L2, dL2 = basis_tables(el2, ξf)
    j1f, j2f = jac1.(ξf), jac2.(ξf)
    n = Nx * Ny
    Mref = zeros(RT, n, n)
    MD1ref = zeros(RT, n, n)
    MD2ref = zeros(RT, n, n)
    for q2 in 1:Ny, q1 in 1:Nx, p2 in 1:Ny, p1 in 1:Nx
        p, q = lin(p1, p2, Nx), lin(q1, q2, Nx)
        sM = sD1 = sD2 = zero(RT)
        for k2 in 1:NQ, k1 in 1:NQ
            wJ = Wf[k1] * Wf[k2] * j1f[k1] * j2f[k2]
            ψp = L1[k1, p1] * L2[k2, p2]
            sM += wJ * ψp * L1[k1, q1] * L2[k2, q2]
            # physical derivatives by the diagonal chain rule ∂_x = (1/j1) ∂_ξ, ∂_y = (1/j2) ∂_η
            sD1 += wJ * ψp * dL1[k1, q1] / j1f[k1] * L2[k2, q2]
            sD2 += wJ * ψp * L1[k1, q1] * dL2[k2, q2] / j2f[k2]
        end
        Mref[p, q] = sM
        MD1ref[p, q] = sD1
        MD2ref[p, q] = sD2
    end
    return Mref, MD1ref, MD2ref
end

# dense fine-quadrature face references for a face of normal dimension `dim`, side s = ±1:
# Iref[p,q̄] = ψ_p(face node q̄), MLref[p,q̄] = ∫ ψ_p ℓ̄_q̄ Ĵ dt, MLDTref[p,q̄] = ∫ ∂_nψ_p ℓ̄_q̄ Ĵ dt
# with Ĵ(t) = transverse Jacobian, ∂_n = (1/j_n(face)) ∂_ξn. Exact when Ĵ is polynomial
# (affine transverse map, or any map on the fine rule when the transverse Jacobian is).
function ref_face(RT, dim, s, eln, elt, jacn, jact, Nx, Ny; NQ = 2 * max(Nx, Ny) + 4)
    ξf, Wf = gl_nodes_weights(RT, NQ)
    Lt, _ = basis_tables(elt, ξf)
    rowf = vec(interpolation_matrix(eln, [RT(s)]))          # ℓ_k(±1), normal dim
    drowf = vec(interpolation_matrix(eln, [RT(s)]) * eln.D) # ℓ′_k(±1), normal dim
    jn_face = jacn(RT(s))
    jtf = jact.(ξf)
    Nn, Nt = length(rowf), length(elt.nodes)
    n = Nx * Ny
    Iref = zeros(RT, n, Nt)
    MLref = zeros(RT, n, Nt)
    MLDTref = zeros(RT, n, Nt)
    Ttab = interpolation_matrix(elt, elt.nodes)             # exact unit rows (collision path)
    for q̄ in 1:Nt, pt in 1:Nt, pn in 1:Nn
        p = dim == 1 ? lin(pn, pt, Nx) : lin(pt, pn, Nx)
        sML = sMLDT = zero(RT)
        for k in 1:NQ
            w = Wf[k] * jtf[k] * Lt[k, pt] * Lt[k, q̄]
            sML += w * rowf[pn]
            sMLDT += w * drowf[pn] / jn_face
        end
        Iref[p, q̄] = rowf[pn] * Ttab[q̄, pt]
        MLref[p, q̄] = sML
        MLDTref[p, q̄] = sMLDT
    end
    return Iref, MLref, MLDTref
end

# definitional (collocation-rule) face reference, for faces whose transverse Jacobian is not
# polynomial (main.tex:636 defines the discrete integrals BY the collocation Gauss rule):
# entries directly from the eq. formulas with independently recomputed geometry.
function ref_face_defn(RT, dim, s, eln, elt, jacn, jact, Nx, Ny)
    rowf = vec(interpolation_matrix(eln, [RT(s)]))
    drowf = vec(interpolation_matrix(eln, [RT(s)]) * eln.D)
    jn_face = jacn(RT(s))
    Nn, Nt = length(rowf), length(elt.nodes)
    n = Nx * Ny
    Iref = zeros(RT, n, Nt)
    MLref = zeros(RT, n, Nt)
    MLDTref = zeros(RT, n, Nt)
    for q̄ in 1:Nt, pn in 1:Nn
        p = dim == 1 ? lin(pn, q̄, Nx) : lin(q̄, pn, Nx)
        sJw = jact(elt.nodes[q̄]) * elt.weights[q̄]
        Iref[p, q̄] = rowf[pn]
        MLref[p, q̄] = rowf[pn] * sJw
        MLDTref[p, q̄] = drowf[pn] * sJw / jn_face
    end
    return Iref, MLref, MLDTref
end

# 2D polynomial in ascending powers: val(c, x, y) = Σ_ab c[a+1,b+1] x^a y^b
polyval2(c, x, y) = sum(c[a+1, b+1] * x^a * y^b for a in 0:size(c, 1)-1, b in 0:size(c, 2)-1)
polydx2(c, x, y) = sum(a * c[a+1, b+1] * x^(a - 1) * y^b for a in 1:size(c, 1)-1, b in 0:size(c, 2)-1;
                       init = zero(x))
polydy2(c, x, y) = sum(b * c[a+1, b+1] * x^a * y^(b - 1) for a in 0:size(c, 1)-1, b in 1:size(c, 2)-1;
                       init = zero(x))
randmat!(rng, n1, n2) = [dgm_unit!(rng) for _ in 1:n1, _ in 1:n2]

@testset "dg_matrices (dg-operator-matrices)" begin

# fixed fixtures: one affine and one compaction-mapped element, asymmetric sizes AND
# asymmetric refinement levels per dimension to catch any index/dimension transposition.
affine_blk() = Block(3, LinearMap1D(2.8, 5.2), LinearMap1D(0.0, 0.25))
affine_id() = ElementId(3, SegmentId(2, 1), SegmentId(1, 0))
compaction_blk() = Block(5, InverseRadialMap1D(25.3, 1.0e4), LinearMap1D(0.0, 0.25))
compaction_id() = ElementId(5, SegmentId(2, 3), SegmentId(1, 1))   # steepest radial quarter

@testset "affine element: every matrix vs dense fine-quadrature reference" begin
    blk, id = affine_blk(), affine_id()
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))
    el1, el2 = refelement(Float64, Nx), refelement(Float64, Ny)
    phys1, jac1 = el_geometry(Float64, blockmap(blk, 1), id, 1)
    phys2, jac2 = el_geometry(Float64, blockmap(blk, 2), id, 2)

    @test em.N == (Nx, Ny)
    @test size(em.M) == (Nx * Ny, Nx * Ny)
    @test em.M isa Diagonal
    @test all(>(0), em.M.diag)
    @test em.x[1] ≈ phys1.(el1.nodes) rtol = 1e-14
    @test em.x[2] ≈ phys2.(el2.nodes) rtol = 1e-14
    @test em.jac[1] ≈ jac1.(el1.nodes) rtol = 1e-14
    @test em.jac[2] ≈ jac2.(el2.nodes) rtol = 1e-14

    Mref, MD1ref, MD2ref = ref_volume(Float64, el1, el2, jac1, jac2, Nx, Ny)
    # mass-lumping theorem: the EXACT integral matrix is diagonal for affine maps
    @test maximum(abs, Mref - Diagonal(diag(Mref))) <= 1e-13 * maximum(abs, Mref)
    @test relerr(Matrix(em.M), Mref) <= 1e-13
    @test relerr(em.MD[1], MD1ref) <= 1e-13
    @test relerr(em.MD[2], MD2ref) <= 1e-13

    for dim in 1:2, s in (-1, 1)
        f = face_matrices(em, dim, s)
        @test f.dim == dim && f.side == s
        eln, elt = dim == 1 ? (el1, el2) : (el2, el1)
        jacn, jact = dim == 1 ? (jac1, jac2) : (jac2, jac1)
        Iref, MLref, MLDTref = ref_face(Float64, dim, s, eln, elt, jacn, jact, Nx, Ny)
        nt = dim == 1 ? Ny : Nx
        @test size(f.I) == (Nx * Ny, nt)
        @test f.njac ≈ jacn(float(s)) rtol = 1e-14
        @test relerr(f.I, Iref) <= 1e-13
        @test relerr(f.ML, MLref) <= 1e-13
        @test relerr(f.MLDT, MLDTref) <= 1e-13
    end
end

@testset "compaction element: definitional + analytic references" begin
    blk, id = compaction_blk(), compaction_id()
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))
    el1, el2 = refelement(Float64, Nx), refelement(Float64, Ny)
    phys1, jac1 = el_geometry(Float64, blockmap(blk, 1), id, 1)
    phys2, jac2 = el_geometry(Float64, blockmap(blk, 2), id, 2)

    # M vs the pointwise definition, independently recomputed geometry (main.tex:649)
    Mdef = [jac1(el1.nodes[i]) * jac2(el2.nodes[j]) * el1.weights[i] * el2.weights[j]
            for j in 1:Ny, i in 1:Nx]   # (j outer, i inner) — then flattened i-fastest below
    @test relerr(em.M.diag, vec(permutedims(Mdef))) <= 1e-13

    # MD_1 (radial) vs the EXACT analytic fine-quadrature reference: the 1/j1 chain-rule
    # factor cancels J's compaction part, leaving a polynomial integrand (j2 constant here)
    _, MD1ref, _ = ref_volume(Float64, el1, el2, jac1, jac2, Nx, Ny)
    @test relerr(em.MD[1], MD1ref) <= 1e-13
    # MD_2 vs the definitional collocation form with independent geometry (main.tex:650)
    MD2def = zeros(Float64, Nx * Ny, Nx * Ny)
    for q2 in 1:Ny, p2 in 1:Ny, i in 1:Nx
        MD2def[lin(i, p2, Nx), lin(i, q2, Nx)] =
            el1.weights[i] * el2.weights[p2] * jac1(el1.nodes[i]) * el2.D[p2, q2]
    end
    @test relerr(em.MD[2], MD2def) <= 1e-13

    # dim-1 faces: transverse Jacobian is the (constant) angular one → fine rule exact;
    # dim-2 faces: transverse Jacobian is the compaction one → definitional reference
    for s in (-1, 1)
        f1 = face_matrices(em, 1, s)
        Iref, MLref, MLDTref = ref_face(Float64, 1, s, el1, el2, jac1, jac2, Nx, Ny)
        @test relerr(f1.I, Iref) <= 1e-13
        @test relerr(f1.ML, MLref) <= 1e-13
        @test relerr(f1.MLDT, MLDTref) <= 1e-13
        f2 = face_matrices(em, 2, s)
        Iref, MLref, MLDTref = ref_face_defn(Float64, 2, s, el2, el1, jac2, jac1, Nx, Ny)
        @test relerr(f2.I, Iref) <= 1e-13
        @test relerr(f2.ML, MLref) <= 1e-13
        @test relerr(f2.MLDT, MLDTref) <= 1e-13
    end
end

sbp_fixtures = (("affine", affine_blk(), affine_id()),
                ("compaction", compaction_blk(), compaction_id()))
@testset "2D Jacobian-weighted generalized SBP + M·Dphys consistency [$name]" for (name, blk, id) in sbp_fixtures
    for (Nx, Ny) in ((2, 3), (5, 4), (1, 4))
        em = element_matrices(Float64, blk, id, (Nx, Ny))
        el1, el2 = refelement(Float64, Nx), refelement(Float64, Ny)
        for dim in 1:2
            lo, hi = face_matrices(em, dim, -1), face_matrices(em, dim, +1)
            lhs = em.MD[dim] + transpose(em.MD[dim])
            rhs = hi.ML * transpose(hi.I) - lo.ML * transpose(lo.I)
            scale = max(1.0, maximum(abs, rhs))
            @test maximum(abs, lhs - rhs) <= 1e-13 * scale
            # falsifier: the identity has teeth — wrong face orientation breaks it badly
            # (skip Nx==1 or Ny==1: a single Gauss node sits at ξ=0 where both lifts coincide)
            if Nx > 1 && Ny > 1
                wrong = lo.ML * transpose(lo.I) - hi.ML * transpose(hi.I)
                @test maximum(abs, lhs - wrong) > 1e-3 * scale
            end
            # MD_i == M · Dphys_i (the GLOBAL_DAG row-10 reading of the identity)
            D1d = dim == 1 ? el1.D : el2.D
            n = Nx * Ny
            Dphys = zeros(Float64, n, n)
            for t in 1:(dim == 1 ? Ny : Nx), a in 1:size(D1d, 1), b in 1:size(D1d, 2)
                p = dim == 1 ? lin(a, t, Nx) : lin(t, a, Nx)
                q = dim == 1 ? lin(b, t, Nx) : lin(t, b, Nx)
                Dphys[p, q] = D1d[a, b] / em.jac[dim][a]
            end
            @test relerr(em.MD[dim], Matrix(em.M) * Dphys) <= 1e-13
            # D annihilates constants ⇒ MD_i · 1 == 0 (discrete constant-flux divergence)
            @test maximum(abs, em.MD[dim] * ones(n)) <=
                  1e-13 * max(1.0, maximum(abs, em.MD[dim]))
        end
    end
end

@testset "polynomial-flux divergence integrated exactly (affine, physical coords)" begin
    rng = DGMRNG(0xd6_0b5e_55f0_2026_0a)
    blk, id = affine_blk(), affine_id()
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))
    el1, el2 = refelement(Float64, Nx), refelement(Float64, Ny)
    phys1, jac1 = el_geometry(Float64, blockmap(blk, 1), id, 1)
    phys2, jac2 = el_geometry(Float64, blockmap(blk, 2), id, 2)
    c1, c2 = randmat!(rng, Nx, Ny), randmat!(rng, Nx, Ny)   # degrees ≤ (Nx−1, Ny−1)
    x1, x2 = phys1.(el1.nodes), phys2.(el2.nodes)
    F1 = [polyval2(c1, x1[i], x2[j]) for j in 1:Ny for i in 1:Nx]
    F2 = [polyval2(c2, x1[i], x2[j]) for j in 1:Ny for i in 1:Nx]
    lhs = em.MD[1] * F1 + em.MD[2] * F2

    # reference: ∫ ψ_p (∂_xF¹ + ∂_yF²) J dξ² on the fine rule (polynomial → exact)
    NQ = 2 * max(Nx, Ny) + 4
    ξf, Wf = gl_nodes_weights(Float64, NQ)
    L1, _ = basis_tables(el1, ξf)
    L2, _ = basis_tables(el2, ξf)
    ref = zeros(Float64, Nx * Ny)
    for p2 in 1:Ny, p1 in 1:Nx
        acc = 0.0
        for k2 in 1:NQ, k1 in 1:NQ
            xk, yk = phys1(ξf[k1]), phys2(ξf[k2])
            divF = polydx2(c1, xk, yk) + polydy2(c2, xk, yk)
            acc += Wf[k1] * Wf[k2] * jac1(ξf[k1]) * jac2(ξf[k2]) *
                   L1[k1, p1] * L2[k2, p2] * divF
        end
        ref[lin(p1, p2, Nx)] = acc
    end
    @test relerr(lhs, ref) <= 1e-13

    # discrete divergence theorem: grid sum == exact surface integral ∮ F·n Ĵ
    surf = 0.0
    for k in 1:NQ
        surf += Wf[k] * jac2(ξf[k]) * (polyval2(c1, phys1(1.0), phys2(ξf[k])) -
                                       polyval2(c1, phys1(-1.0), phys2(ξf[k])))
        surf += Wf[k] * jac1(ξf[k]) * (polyval2(c2, phys1(ξf[k]), phys2(1.0)) -
                                       polyval2(c2, phys1(ξf[k]), phys2(-1.0)))
    end
    @test abs(sum(lhs) - surf) <= 1e-13 * max(1.0, abs(surf))
end

@testset "polynomial-flux divergence integrated exactly (compaction, radial/logical)" begin
    rng = DGMRNG(0x5eed_face_d6b2_0001)
    blk, id = compaction_blk(), compaction_id()
    Nx, Ny = 5, 4
    em = element_matrices(Float64, blk, id, (Nx, Ny))
    el1, el2 = refelement(Float64, Nx), refelement(Float64, Ny)
    _, jac1 = el_geometry(Float64, blockmap(blk, 1), id, 1)
    _, jac2 = el_geometry(Float64, blockmap(blk, 2), id, 2)
    c = randmat!(rng, Nx, Ny)                                # F(ξ,η), degrees ≤ (Nx−1, Ny−1)
    F = [polyval2(c, el1.nodes[i], el2.nodes[j]) for j in 1:Ny for i in 1:Nx]
    lhs = em.MD[1] * F

    # (MD_1 F)_p == ∫ ℓ_p1 ℓ_p2 F_ξ(ξ,η) j2(η) dξ dη — polynomial (j2 constant) → exact
    NQ = 2 * max(Nx, Ny) + 4
    ξf, Wf = gl_nodes_weights(Float64, NQ)
    L1, _ = basis_tables(el1, ξf)
    L2, _ = basis_tables(el2, ξf)
    ref = zeros(Float64, Nx * Ny)
    for p2 in 1:Ny, p1 in 1:Nx
        acc = 0.0
        for k2 in 1:NQ, k1 in 1:NQ
            acc += Wf[k1] * Wf[k2] * jac2(ξf[k2]) * L1[k1, p1] * L2[k2, p2] *
                   polydx2(c, ξf[k1], ξf[k2])
        end
        ref[lin(p1, p2, Nx)] = acc
    end
    @test relerr(lhs, ref) <= 1e-13

    # radial discrete divergence theorem at the collocation face nodes (definitional rule)
    surf = 0.0
    for j in 1:Ny
        surf += jac2(el2.nodes[j]) * el2.weights[j] *
                (polyval2(c, 1.0, el2.nodes[j]) - polyval2(c, -1.0, el2.nodes[j]))
    end
    @test abs(sum(lhs) - surf) <= 1e-13 * max(1.0, abs(surf))
end

bf_blocks = (("affine", () -> Block(3, LinearMap1D(big"2.8", big"5.2"),
                                    LinearMap1D(big"0.0", big"0.25"))),
             ("compaction", () -> Block(3, InverseRadialMap1D(big"25.3", big"1e4"),
                                        LinearMap1D(big"0.0", big"0.25"))))
@testset "BigFloat(256) pass-through [$name]" for (name, mkblk) in bf_blocks
    setprecision(BigFloat, 256) do
        blk = mkblk()
        id = ElementId(3, SegmentId(2, 3), SegmentId(1, 1))
        Nx, Ny = 4, 3
        em = element_matrices(BigFloat, blk, id, (Nx, Ny))
        @test eltype(em.M.diag) === BigFloat
        @test precision(em.M.diag[1]) == 256
        @test precision(em.MD[1][1, 2]) == 256
        el1 = refelement(BigFloat, Nx)
        _, jac1 = el_geometry(BigFloat, blockmap(blk, 1), id, 1)
        _, jac2 = el_geometry(BigFloat, blockmap(blk, 2), id, 2)
        # pointwise-definition M with independent geometry
        el2 = refelement(BigFloat, Ny)
        for j in 1:Ny, i in 1:Nx
            mij = jac1(el1.nodes[i]) * jac2(el2.nodes[j]) * el1.weights[i] * el2.weights[j]
            @test abs(em.M.diag[lin(i, j, Nx)] - mij) <= big"1e-70" * abs(mij)
        end
        # the 2D generalized-SBP identity at 256-bit accuracy (Float64 cannot pass this)
        for dim in 1:2
            lo, hi = face_matrices(em, dim, -1), face_matrices(em, dim, +1)
            lhs = em.MD[dim] + transpose(em.MD[dim])
            rhs = hi.ML * transpose(hi.I) - lo.ML * transpose(lo.I)
            @test maximum(abs, lhs - rhs) <= big"1e-70" * max(big"1.0", maximum(abs, rhs))
        end
    end
end

@testset "cache: 1D constituents drawn from the (RT, N)-keyed refelement cache" begin
    EllipticCore.refelement_cache_reset!()
    blkA, blkB = affine_blk(), compaction_blk()
    ids = Any[(blkA, ElementId(3, SegmentId(l, k), SegmentId(1, 0))) for (l, k) in
              ((0, 0), (1, 1), (2, 2))]
    append!(ids, Any[(blkB, ElementId(5, SegmentId(l, k), SegmentId(1, 1))) for (l, k) in
                     ((1, 0), (2, 3), (3, 5))])
    ems = [element_matrices(Float64, b, i, (5, 4)) for (b, i) in ids]
    # first build misses (Float64,5) and (Float64,4); the other 5 builds hit both
    @test EllipticCore.refelement_cache_stats() == (hits = 10, misses = 2)
    # ComplexF64 request normalizes to the Float64 entries (2 more hits, 0 misses) and
    # produces identical REAL-typed matrices (geometry lives in RT = real(T))
    emc = element_matrices(ComplexF64, blkA, ids[1][2], (5, 4))
    @test EllipticCore.refelement_cache_stats() == (hits = 12, misses = 2)
    @test eltype(emc.M.diag) === Float64
    @test emc.M == ems[1].M && emc.MD[1] == ems[1].MD[1] && emc.MD[2] == ems[1].MD[2]
    # distinct elements on the same block genuinely differ (Jacobian scaling by level)
    @test ems[1].M.diag != ems[2].M.diag
end

@testset "argument validation" begin
    blk, id = affine_blk(), affine_id()
    em = element_matrices(Float64, blk, id, (3, 2))
    @test_throws ArgumentError element_matrices(Float64, blk, ElementId(4, SegmentId(0, 0),
                                                                        SegmentId(0, 0)), (3, 2))
    @test_throws ArgumentError element_matrices(Float64, blk, id, (0, 2))
    @test_throws ArgumentError face_matrices(em, 3, 1)
    @test_throws ArgumentError face_matrices(em, 1, 0)
    @test face_matrices(em, 1, -1) === em.faces[1][1]
    @test face_matrices(em, 2, +1) === em.faces[2][2]
end

end # testset dg_matrices
