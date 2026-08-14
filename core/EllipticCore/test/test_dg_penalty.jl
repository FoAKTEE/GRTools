# test_dg_penalty.jl — grounded gates for src/dg_penalty.jl (DAG node: dg-penalty-flux,
# elliptic-ssf GLOBAL_DAG row 12 + TOPOLOGICAL_SCHEDULE W4-B, CRITICAL PATH). Written BEFORE
# the implementation (TDD red run captured in evidence/dg-penalty-flux/pkg_test_red.log).
#
# Target: the internal penalty numerical flux of CRP eq:numflux (main.tex:583-595),
#
#   Ψ*    = {{Ψ}} = (ψ_int + ψ_ext)/2,
#   (nF)* = {{nF}} − σ n_i 𝔄^{ij} n_j (ψ_int − ψ_ext),
#
# with the penalty law from Vu:2024 (ref-paper/arXiv-2405.06120v2/paper.tex:594-619):
#
#   σ = C (max(p_int, p_ext) + 1)² / min(h_int, h_ext),   C = 1.5 (CRP main.tex:595),
#   h = 2·𝒥/𝒥^surf = the element size perpendicular to the face,
#
# evaluated at the face (for the tensor-product diagonal maps h == 2·njac, the composed
# normal Jacobian at the face — exactly the physical element extent along the face-normal
# dimension on affine maps); plus the first external-boundary data instance: the
# zero-Dirichlet ghost u_ext = 2u_D − u_int|_{u_D=0}, F_ext = F_int (Vu:2024 eq:dirichlet,
# audit-speccore §"ubc ghost"), in the SENDER convention of dg_operator.jl N3:
# zero_dirichlet_ghost(fs) = FaceState(−ψ_int, −φ_int).
#
# Oracle gates (worker node contract, constants calibrated PRE-RED in scratchpad
# calibrate_dg_penalty.jl + calibrate_dg_penalty2.jl from the sealed wave-2/3/4A layers
# plus a throwaway flux prototype):
#   1. HEADLINE — first end-to-end solve: manufactured Poisson (𝔄=I, zero Dirichlet) on a
#      4×4-element affine mesh, massive system R·u = M·S solved with the wave-3 complex
#      GMRES as Complex{Float64}: L2 error decays EXPONENTIALLY over N = 4, 6, 8
#      (measured 6.79e-5 → 8.01e-8 → 5.31e-11; ratios ≤ 1e-2 gate vs ≤ 1.2e-3 measured;
#      LSQ log10 slope ≤ −1.0 per point of N vs −1.53 measured), err(N=8) < 1e-8
#      (margin ×188); GMRES converged to abstol 1e-12 (12/39/86 iterations).
#   2. complex Helmholtz variant β = −k² − iε (k=5, ε=0.4) on the same mesh, N=6: solved
#      to abstol 1e-11 in < 200 iterations (measured 37), L2 error < 1e-6 (measured 8.1e-8).
#   3. SIPG symmetry probe: global matrix column-by-column via unit vectors, 2×1 elements,
#      N=3, real Poisson: ‖A−Aᵀ‖∞/‖A‖∞ < 1e-12 (measured 1.1e-16); C=1.5 positive definite
#      (min eig of sym part = 1.336 > 0.5 gate) while C=0 is indefinite (−20.5 < −1 gate).
#   4. penalty-scaling structure: A(C) affine-linear in C — A(3.0)−A(1.5) == A(4.5)−A(3.0)
#      to 1e-12 rel (measured 6.5e-16) and equal to the independently assembled penalty
#      block structure Δσ·ML·diag(𝔄face)·Iᵀ (own/neighbor columns, ghost faces doubled)
#      to 1e-12 rel (measured 2.7e-16).
#   5. nonuniform-mesh spectral convergence (decomp-crp dg-penalty-flux gate iii): same
#      Poisson on an inverse-radial×linear block (element-width spread ×4, shared-face
#      njac bitwise-equal): errors 8.31e-4 → 6.38e-5 → 8.97e-7 over N = 4, 6, 8 (ratio
#      gate ≤ 0.25 vs ≤ 0.077 measured), err(N=8) < 1e-5; sym-part min eig > 0.01
#      (measured 0.092). NOTE the curved-map matrix is NOT symmetric to machine precision
#      (collocation aliasing, measured rel 1.4e-4) — machine-symmetry is an AFFINE gate.
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test files,
# ---- distinct name/seed) ---------------------------------------------------------------
mutable struct DGPRNG
    s::UInt64
end
dgp_unit!(r::DGPRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                        Int64(r.s >> 11) / 2.0^53 * 2 - 1)
dgp_vec!(r::DGPRNG, n) = [dgp_unit!(r) for _ in 1:n]
dgp_cvec!(r::DGPRNG, n) = [complex(dgp_unit!(r), dgp_unit!(r)) for _ in 1:n]

const DGP_SIDES = (-1, 1)

# ---- conforming one-block mesh harness (mesh-level orchestration is deliberately
# ---- test-side here, as in test_dg_operator.jl: connectivity containers belong to the
# ---- dg-mortars / matrix-free-operator nodes) ------------------------------------------

# nelx × nely single-block mesh (dyadic levels log2(nel)), N×N points per element
function dgp_setup(blk, nelx, nely, N, A, β, γ)
    lx, ly = Int(log2(nelx)), Int(log2(nely))
    ids = vec([ElementId(blk.id, SegmentId(lx, ix), SegmentId(ly, iy))
               for ix in 0:nelx-1, iy in 0:nely-1])
    ems = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in ids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, A, β, γ) for id in ids)
    off = Dict(id => (k - 1) * N * N for (k, id) in enumerate(ids))
    return ids, ems, cs, off
end

# global matrix-free apply of the massive residual y = R(u): two-pass sweep — every element
# builds its state, then residuals consume neighbor states via the exterior_state plumbing;
# block-boundary faces are closed with the zero-Dirichlet ghost (the node's BC instance)
function dgp_make_apply(ids, ems, cs, off, N, flux)
    npe = N * N
    return function (y, u)
        sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+npe))
                   for id in ids)
        for id in ids
            ext = [[begin
                        nb = neighbor_id(id, d, DGP_SIDES[sidx])
                        nb === nothing ? zero_dirichlet_ghost(sts[id].faces[d][sidx]) :
                            exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1])
                    end for sidx in 1:2] for d in 1:2]
            y[off[id]+1:off[id]+npe] .= strong_residual(ems[id], cs[id], sts[id], ext, flux)
        end
        return y
    end
end

# massive right-hand side b = M·S per element
function dgp_rhs(ids, ems, off, N, S)
    b = zeros(ComplexF64, length(ids) * N * N)
    for id in ids
        em = ems[id]
        b[off[id]+1:off[id]+N*N] .=
            em.M * vec([S(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
    end
    return b
end

# discrete (mass-weighted) L2 error of a global solution vector vs the exact solution
function dgp_l2err(ids, ems, off, N, u, uex)
    s = 0.0
    for id in ids
        em = ems[id]
        e = u[off[id]+1:off[id]+N*N] .-
            vec([uex(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
        s += real(dot(e, em.M * e))
    end
    return sqrt(s)
end

# dense global matrix column-by-column through the matrix-free apply (unit-vector probes)
function dgp_build_matrix(apply!, n)
    A = zeros(ComplexF64, n, n)
    y = zeros(ComplexF64, n)
    u = zeros(ComplexF64, n)
    for j in 1:n
        u[j] = 1
        apply!(y, u)
        A[:, j] = y
        u[j] = 0
    end
    return A
end

# least-squares slope of log10(err) against the points-per-dimension ladder
function dgp_lsq_slope(Ns, errs)
    x = float.(collect(Ns))
    y = log10.(errs)
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

@testset "dg_penalty (dg-penalty-flux)" begin

# constant / zero closures + the manufactured Poisson problem used throughout
dgp_one2(x, y) = 1.0
dgp_zero2(x, y) = 0.0
dgp_uex(x, y) = sinpi(x) * sinpi(y)                 # zero Dirichlet on [0,1]²
dgp_S(x, y) = 2 * pi^2 * dgp_uex(x, y)              # −∇²u = S
dgp_unit_blk() = Block(31, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))

@testset "penalty σ law + perpendicular size (Vu:2024 paper.tex:594-619)" begin
    # σ = C (max(p_int, p_ext) + 1)² / min(h_int, h_ext) — exact small-number arithmetic
    @test penalty_sigma(1.5, 5, 5, 0.25, 0.25) == 1.5 * 36 / 0.25
    @test penalty_sigma(1.0, 3, 3, 0.5, 0.5) == 32.0
    # max/min pairing + exchange symmetry (σ must be identical from both sides of a face)
    @test penalty_sigma(1.5, 3, 5, 0.5, 0.25) == 1.5 * 36 / 0.25
    @test penalty_sigma(1.5, 5, 3, 0.25, 0.5) == penalty_sigma(1.5, 3, 5, 0.5, 0.25)
    # linear in C
    @test penalty_sigma(3.0, 4, 4, 0.125, 0.125) == 2 * penalty_sigma(1.5, 4, 4, 0.125, 0.125)
    # BigFloat pass-through (precision inherited, DESIGN §4)
    σb = penalty_sigma(big"1.5", 5, 5, big"0.25", big"0.25")
    @test σb isa BigFloat && σb == 216

    # h = 2·𝒥/𝒥^surf == 2·njac: on affine maps this is EXACTLY the element's physical
    # extent along the face-normal dimension (block [0,2]×[0,1], level (1,1) elements)
    blk = Block(32, LinearMap1D(0.0, 2.0), LinearMap1D(0.0, 1.0))
    id = ElementId(32, SegmentId(1, 0), SegmentId(1, 1))
    em = element_matrices(Float64, blk, id, (4, 5))
    @test perpendicular_size(face_matrices(em, 1, -1)) == 1.0     # Δx = 2/2
    @test perpendicular_size(face_matrices(em, 1, +1)) == 1.0
    @test perpendicular_size(face_matrices(em, 2, -1)) == 0.5     # Δy = 1/2
    @test perpendicular_size(face_matrices(em, 2, +1)) == 0.5
    @test perpendicular_size(face_matrices(em, 1, -1)) == 2 * face_matrices(em, 1, -1).njac

    # shared-face h and p are BITWISE equal between conforming same-level neighbors on
    # affine AND inverse-radial maps (same dyadic logical endpoint through the same block
    # map) — the receiver-side σ evaluation IS the pairwise law on conforming interfaces
    for blk2 in (dgp_unit_blk(), Block(33, InverseRadialMap1D(1.0, 3.0), LinearMap1D(0.0, 1.0)))
        ids, ems, _, _ = dgp_setup(blk2, 4, 4, 5, (dgp_one2, dgp_one2), dgp_zero2,
                                   (dgp_zero2, dgp_zero2))
        for id2 in ids, d in 1:2, sidx in 1:2
            nb = neighbor_id(id2, d, DGP_SIDES[sidx])
            nb === nothing && continue
            fint = face_matrices(ems[id2], d, DGP_SIDES[sidx])
            fext = face_matrices(ems[nb], d, -DGP_SIDES[sidx])
            @test perpendicular_size(fint) == perpendicular_size(fext)   # bitwise
            @test size(fint.I) == size(fext.I)                           # same p
        end
    end
end

@testset "PenaltyFlux formula, consistency, central limit" begin
    @test PenaltyFlux() isa NumericalFlux
    @test PenaltyFlux().C == 1.5                     # CRP main.tex:595 default
    @test PenaltyFlux(3.0).C == 3.0

    blk = dgp_unit_blk()
    id = ElementId(31, SegmentId(1, 1), SegmentId(1, 0))
    N = 5
    em = element_matrices(Float64, blk, id, (N, N))
    rng = DGPRNG(0x00000000dead0031)
    for d in 1:2, sidx in 1:2
        fm = face_matrices(em, d, DGP_SIDES[sidx])
        nf = size(fm.I, 2)
        Af = dgp_cvec!(rng, nf)                      # complex 𝔄 at the face (m-mode case)
        int = FaceState(dgp_cvec!(rng, nf), dgp_cvec!(rng, nf))
        ext = FaceState(dgp_cvec!(rng, nf), dgp_cvec!(rng, nf))
        C = 1.5
        σ = penalty_sigma(C, N - 1, N - 1, 2 * fm.njac, 2 * fm.njac)
        ψs, φs = numerical_flux(PenaltyFlux(C), fm, Af, int, ext)
        # hand formula: Ψ* = {{Ψ}}; (nF)* = {{nF}} − σ 𝔄face (ψint − ψext) (receiver normal;
        # the sender-convention ext.nflux enters with the N3/N4 minus sign)
        @test ψs == (int.psi .+ ext.psi) ./ 2
        φref = (int.nflux .- ext.nflux) ./ 2 .- σ .* (Af .* (int.psi .- ext.psi))
        @test maximum(abs, φs .- φref) ≤ 1e-13 * maximum(abs, φref)

        # zero jump ⇒ penalty term OFF: identical ψ traces give exactly the central flux
        ext0 = FaceState(copy(int.psi), dgp_cvec!(rng, nf))
        ψc, φc = numerical_flux(CentralFlux(), fm, Af, int, ext0)
        ψp, φp = numerical_flux(PenaltyFlux(C), fm, Af, int, ext0)
        @test ψp == ψc && φp == φc

        # C → 0 limit IS the central flux (bitwise: σ = 0 kills the penalty term)
        ψ0, φ0 = numerical_flux(PenaltyFlux(0.0), fm, Af, int, ext)
        ψc2, φc2 = numerical_flux(CentralFlux(), fm, Af, int, ext)
        @test ψ0 == ψc2 && φ0 == φc2
    end

    # exact-solution flux consistency (decomp-crp dg-penalty-flux gate i): a polynomial the
    # basis represents EXACTLY, vanishing on ∂[0,1]² — all interface AND ghost-boundary
    # jumps are pure roundoff, so the penalty residual equals both the central residual and
    # the analytic M·S at machine precision (calibrated 1.6e-14 rel, σ-amplified roundoff)
    up(x, y) = x * (1 - x) * y * (1 - y)
    Sp(x, y) = 2 * x * (1 - x) + 2 * y * (1 - y)     # −∇²u, exact
    N2 = 5
    ids, ems, cs, off = dgp_setup(blk, 2, 2, N2, (dgp_one2, dgp_one2), dgp_zero2,
                                  (dgp_zero2, dgp_zero2))
    n = length(ids) * N2 * N2
    u = zeros(ComplexF64, n)
    bref = zeros(ComplexF64, n)
    for id2 in ids
        em2 = ems[id2]
        u[off[id2]+1:off[id2]+N2*N2] .=
            vec([up(em2.x[1][i], em2.x[2][j]) for i in 1:N2, j in 1:N2])
        bref[off[id2]+1:off[id2]+N2*N2] .=
            em2.M * vec([Sp(em2.x[1][i], em2.x[2][j]) for i in 1:N2, j in 1:N2])
    end
    yp = zeros(ComplexF64, n)
    yc = zeros(ComplexF64, n)
    dgp_make_apply(ids, ems, cs, off, N2, PenaltyFlux())(yp, u)
    dgp_make_apply(ids, ems, cs, off, N2, CentralFlux())(yc, u)
    scale = maximum(abs, bref)
    @test maximum(abs, yp .- bref) ≤ 1e-11 * scale
    @test maximum(abs, yp .- yc) ≤ 1e-11 * scale
end

@testset "zero-Dirichlet ghost (Vu:2024 eq:dirichlet, sender convention N3)" begin
    rng = DGPRNG(0x00000000dead0032)
    ψ = dgp_cvec!(rng, 5)
    φ = dgp_cvec!(rng, 5)
    g = zero_dirichlet_ghost(FaceState(ψ, φ))
    @test g isa FaceState
    @test g.psi == -ψ                    # u_ext = 2·0 − u_int
    @test g.nflux == -φ                  # F_ext = F_int ⇒ ghost's OWN normal flips the sign

    # through the flux: Dirichlet-on-average Ψ* == 0 exactly; {{nF}} == φ_int; the penalty
    # part rides on the doubled jump 2ψ_int
    blk = dgp_unit_blk()
    em = element_matrices(Float64, blk, ElementId(31, SegmentId(0, 0), SegmentId(0, 0)), (5, 5))
    fm = face_matrices(em, 1, -1)
    Af = dgp_cvec!(rng, 5)
    int = FaceState(ψ, φ)
    ψc, φc = numerical_flux(CentralFlux(), fm, Af, int, zero_dirichlet_ghost(int))
    @test all(iszero, ψc)
    @test φc == φ
    σ = penalty_sigma(1.5, 4, 4, 2 * fm.njac, 2 * fm.njac)
    ψp, φp = numerical_flux(PenaltyFlux(), fm, Af, int, zero_dirichlet_ghost(int))
    @test all(iszero, ψp)
    φpref = φ .- σ .* (Af .* (2 .* ψ))
    @test maximum(abs, φp .- φpref) ≤ 1e-13 * maximum(abs, φpref)

    # real input stays real through the ghost (no spurious complexification)
    gr = zero_dirichlet_ghost(FaceState(dgp_vec!(rng, 4), dgp_vec!(rng, 4)))
    @test eltype(gr.psi) == Float64
end

@testset "SIPG symmetry probe: 2×1 elements, N=3, real Poisson" begin
    blk = dgp_unit_blk()
    N = 3
    lids = [ElementId(31, SegmentId(1, ix), SegmentId(0, 0)) for ix in 0:1]
    ems = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in lids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, (dgp_one2, dgp_one2), dgp_zero2,
                                         (dgp_zero2, dgp_zero2)) for id in lids)
    off = Dict(id => (k - 1) * N * N for (k, id) in enumerate(lids))
    probe(flux) = dgp_build_matrix(dgp_make_apply(lids, ems, cs, off, N, flux), 2 * N * N)

    A15c = probe(PenaltyFlux())
    @test maximum(abs, imag.(A15c)) == 0            # real problem ⇒ exactly real columns
    A15 = real.(A15c)
    # SIPG symmetry (the headline symmetry gate; calibrated 1.1e-16)
    @test opnorm(A15 - transpose(A15), Inf) / opnorm(A15, Inf) < 1e-12
    # C = 1.5 penalty-stabilized: positive definite (calibrated min eig 1.336)
    @test minimum(eigvals(Symmetric((A15 + transpose(A15)) / 2))) > 0.5
    # C = 0 (central flux): symmetric but INDEFINITE (calibrated min eig −20.5)
    A0 = real.(probe(CentralFlux()))
    @test opnorm(A0 - transpose(A0), Inf) / opnorm(A0, Inf) < 1e-12
    @test minimum(eigvals(Symmetric((A0 + transpose(A0)) / 2))) < -1

    # ---- penalty-scaling structure: A(C) is affine-linear in C, and the C-slope is
    # ---- EXACTLY the lifted penalty blocks σ(ΔC)·ML·diag(𝔄face)·Iᵀ ----------------------
    A30 = real.(probe(PenaltyFlux(3.0)))
    A45 = real.(probe(PenaltyFlux(4.5)))
    D1 = A30 - A15
    D2 = A45 - A30
    @test opnorm(D1 - D2, Inf) / opnorm(D1, Inf) < 1e-12          # calibrated 6.5e-16

    # independent penalty-difference assembly (ΔC = 1.5): for every face the penalty term
    # lifts −σ(ψ_int − ψ_ext) through +ML·diag(𝔄face): own columns +Δσ·ML·diag(Af)·I_intᵀ,
    # neighbor columns −Δσ·ML·diag(Af)·I_nbᵀ; zero-Dirichlet ghost faces double the own
    # term (jump = 2ψ_int). The MLDT/{{·}} parts of the residual carry NO C-dependence.
    ΔC = 1.5
    nglob = 2 * N * N
    Dexp = zeros(nglob, nglob)
    for id in lids
        em = ems[id]
        for d in 1:2, sidx in 1:2
            fm = face_matrices(em, d, DGP_SIDES[sidx])
            nvol, nf = size(fm.I)
            Δσ = penalty_sigma(ΔC, nvol ÷ nf - 1, nvol ÷ nf - 1,
                               perpendicular_size(fm), perpendicular_size(fm))
            B = fm.ML * Diagonal(real.(cs[id].Aface[d][sidx]))
            rows = off[id]+1:off[id]+N*N
            nb = neighbor_id(id, d, DGP_SIDES[sidx])
            if nb === nothing
                Dexp[rows, rows] .+= 2 * Δσ .* (B * transpose(fm.I))
            else
                fmnb = face_matrices(ems[nb], d, -DGP_SIDES[sidx])
                Dexp[rows, rows] .+= Δσ .* (B * transpose(fm.I))
                Dexp[rows, off[nb]+1:off[nb]+N*N] .-= Δσ .* (B * transpose(fmnb.I))
            end
        end
    end
    @test opnorm(D1 - Dexp, Inf) / opnorm(D1, Inf) < 1e-12        # calibrated 2.7e-16
    @test opnorm(D1 - transpose(D1), Inf) / opnorm(D1, Inf) < 1e-12
end

@testset "HEADLINE: end-to-end Poisson, 4×4 affine mesh, zero Dirichlet, GMRES" begin
    Ns = (4, 6, 8)
    errs = Float64[]
    iters = Int[]
    for N in Ns
        ids, ems, cs, off = dgp_setup(dgp_unit_blk(), 4, 4, N, (dgp_one2, dgp_one2),
                                      dgp_zero2, (dgp_zero2, dgp_zero2))
        ap = dgp_make_apply(ids, ems, cs, off, N, PenaltyFlux())
        b = dgp_rhs(ids, ems, off, N, dgp_S)
        res = gmres(ap, b; abstol = 1e-12, reltol = 0.0, maxiter = 2000, restart = 2000)
        @test res.converged                          # ‖b − R·u‖ ≤ 1e-12 (true residual)
        push!(errs, dgp_l2err(ids, ems, off, N, res.x, dgp_uex))
        push!(iters, res.iterations)
    end
    # exponential (spectral) decay of the discrete L2 error (calibrated
    # 6.79e-5 → 8.01e-8 → 5.31e-11; ratios 1.2e-3, 6.6e-4; LSQ slope −1.53/point)
    @test errs[1] < 3e-4
    @test errs[2] / errs[1] ≤ 1e-2
    @test errs[3] / errs[2] ≤ 1e-2
    @test dgp_lsq_slope(Ns, errs) ≤ -1.0
    @test errs[3] < 1e-8                             # NODE HEADLINE GATE (calibrated 5.3e-11)
    @test all(iters .≤ 500)                          # calibrated 12/39/86
end

@testset "end-to-end complex Helmholtz (β = −k² − iε), < 200 iterations" begin
    k, ε = 5.0, 0.4
    βH = -k^2 - im * ε
    SH(x, y) = (2 * pi^2 + βH) * dgp_uex(x, y)       # (−∇² + β)u = S with u = uex exactly
    N = 6
    ids, ems, cs, off = dgp_setup(dgp_unit_blk(), 4, 4, N, (dgp_one2, dgp_one2),
                                  (x, y) -> βH, (dgp_zero2, dgp_zero2))
    ap = dgp_make_apply(ids, ems, cs, off, N, PenaltyFlux())
    b = dgp_rhs(ids, ems, off, N, SH)
    @test eltype(b) == ComplexF64                    # genuinely complex system
    res = gmres(ap, b; abstol = 1e-11, reltol = 0.0, maxiter = 2000, restart = 2000)
    @test res.converged
    @test res.iterations < 200                       # calibrated 37
    @test dgp_l2err(ids, ems, off, N, res.x, dgp_uex) < 1e-6      # calibrated 8.1e-8
end

@testset "nonuniform inverse-radial mesh: spectral convergence + positivity" begin
    # inverse-radial x ∈ [1,3] (element widths 0.2/0.3/0.5/1.0 at level 2) × linear y ∈ [0,1]
    blk = Block(34, InverseRadialMap1D(1.0, 3.0), LinearMap1D(0.0, 1.0))
    uexn(x, y) = sinpi((x - 1) / 2) * sinpi(y)       # zero Dirichlet on [1,3]×[0,1]
    Sn(x, y) = (5 * pi^2 / 4) * uexn(x, y)
    Ns = (4, 6, 8)
    errs = Float64[]
    for N in Ns
        ids, ems, cs, off = dgp_setup(blk, 4, 4, N, (dgp_one2, dgp_one2), dgp_zero2,
                                      (dgp_zero2, dgp_zero2))
        if N == 4
            # the mesh is genuinely nonuniform: face-h spread across the block is ×4
            hs = sort!([perpendicular_size(ems[id].faces[1][1]) for id in ids])
            @test hs[end] / hs[1] > 3.9
        end
        ap = dgp_make_apply(ids, ems, cs, off, N, PenaltyFlux())
        b = dgp_rhs(ids, ems, off, N, Sn)
        res = gmres(ap, b; abstol = 1e-12, reltol = 0.0, maxiter = 4000, restart = 4000)
        @test res.converged
        push!(errs, dgp_l2err(ids, ems, off, N, res.x, uexn))
    end
    # calibrated 8.31e-4 → 6.38e-5 → 8.97e-7 (ratios 0.077, 0.014; LSQ slope −0.74/point)
    @test errs[2] / errs[1] ≤ 0.25
    @test errs[3] / errs[2] ≤ 0.25
    @test dgp_lsq_slope(Ns, errs) ≤ -0.5
    @test errs[3] < 1e-5

    # penalty keeps the curved-map operator positive definite (calibrated min eig 0.092);
    # NOTE: machine-precision symmetry is an AFFINE property — Gauss collocation aliasing
    # on the curved map leaves a small adjoint inconsistency (calibrated rel 1.4e-4), so
    # the curved gate is positivity + smallness, not 1e-12 symmetry
    N = 4
    ids, ems, cs, off = dgp_setup(blk, 4, 4, N, (dgp_one2, dgp_one2), dgp_zero2,
                                  (dgp_zero2, dgp_zero2))
    An = real.(dgp_build_matrix(dgp_make_apply(ids, ems, cs, off, N, PenaltyFlux()),
                                length(ids) * N * N))
    @test opnorm(An - transpose(An), Inf) / opnorm(An, Inf) < 1e-2
    @test minimum(eigvals(Symmetric((An + transpose(An)) / 2))) > 0.01
end

@testset "BigFloat pass-through" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        blk = Block(35, LinearMap1D(big"0.0", big"1.0"), LinearMap1D(big"0.0", big"1.0"))
        em = element_matrices(BigFloat, blk, ElementId(35, SegmentId(0, 0), SegmentId(0, 0)),
                              (4, 4))
        fm = face_matrices(em, 2, +1)
        @test perpendicular_size(fm) isa BigFloat
        rng = DGPRNG(0x00000000dead0035)
        ψ = BigFloat.(dgp_vec!(rng, 4))
        φ = BigFloat.(dgp_vec!(rng, 4))
        Af = ones(BigFloat, 4)
        ψs, φs = numerical_flux(PenaltyFlux(), fm, Af, FaceState(ψ, φ),
                                zero_dirichlet_ghost(FaceState(ψ, φ)))
        @test eltype(φs) == BigFloat && precision(real(eltype(φs))) == 256
        σ = penalty_sigma(1.5, 3, 3, perpendicular_size(fm), perpendicular_size(fm))
        @test φs ≈ φ .- σ .* (2 .* ψ) rtol = big"1e-70"
    finally
        setprecision(BigFloat, old)
    end
end

end # testset dg_penalty
