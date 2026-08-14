# test_mortars.jl — grounded gates for src/mortars.jl (DAG node: dg-mortars, elliptic-ssf
# GLOBAL_DAG row 13 + TOPOLOGICAL_SCHEDULE W5-A, CRITICAL PATH). Written BEFORE the
# implementation (TDD red run captured in evidence/dg-mortars/pkg_test_red.log).
#
# Target: nonconforming interface coupling per CRP main.tex:678-686 — mortars covering the
# full element face with the LARGER of the two abutting resolutions; both neighbors project
# their face data to the mortar (Lagrange interpolation P), compute the numerical flux
# there, and restrict back with the mass-conservative adjoint R. 1D h/p projection building
# blocks PORTED from specjl SpecCore spectral/mortars.jl (battery-verified there), adapted
# to the Legendre-GAUSS refelement cache: on Gauss nodes the diagonal Gauss mass IS the
# exact L2 Gram matrix (ℓ_aℓ_b has degree 2N−2 ≤ 2N−1), so R = jac·W_f⁻¹·Pᵀ·W_m with
# jac = 1 (full) or 1/2 (half-interval mortar measure). Wired into the wave-4 flux exchange
# through the untouched dg_operator.jl seam (conventions N1-N7 inherited verbatim): the
# nonconforming penalty σ is the PAIRWISE Vu:2024 law σ = C(max(p_int,p_ext)+1)²/
# min(h_int,h_ext) on mortar-paired (p, h) — the W4-B handoff — via the exported
# penalty_sigma.
#
# Oracle gates (worker node contract, constants calibrated PRE-RED in scratchpad
# calibrate_dg_mortars.jl from the sealed wave-2/3/4 layers only):
#   1. R is the mass-conservative adjoint of P w.r.t. the Gauss mass matrices to 1e-13
#      (measured ≤ 2.8e-17), INDEPENDENTLY re-verified as the integral identity
#      ∫_face (Rg)·u = jac·∫_mortar g·(u∘map) with an external fine Gauss rule to 1e-13
#      (measured ≤ 3.0e-15); the Gauss diagonal mass equals the fine-quadrature Gram to
#      1e-14 (measured ≤ 7.8e-16).
#   2. R∘P = I on the coarse space to 1e-14: full mortars (measured ≤ 6.7e-16) and
#      half-pair sums R_lo·P_lo + R_up·P_up (measured ≤ 1.2e-15); P exact on polynomials
#      of degree ≤ N_face−1 to 1e-14 (measured ≤ 4.5e-16).
#   3. polynomial/conforming pass-through: on a conforming face P and R are BITWISE
#      identity matrices, the mortar flux is BITWISE equal to PenaltyFlux, and a whole-mesh
#      apply with mortar coupling on every interior face is BITWISE equal to the direct
#      exchange (all measured true).
#   4. MMS Poisson on a 2:1-refined mesh (one element h-refined once in both dims,
#      neighbor not; 5 elements) with the full penalty flux: spectral L2 decay over
#      N = 4, 6, 8 (measured 8.22e-3 → 1.72e-4 → 1.96e-6, ratios ≤ 0.021, LSQ log10 slope
#      −0.906/point), and NO O(1) interface error: max nodal error on the mortar-interface
#      layer stays within ×4 of the interior max (measured ratio ≤ 1.42) and decays at the
#      same spectral rate (interface ratios ≤ 0.025). Conservation across the 2:1
#      interface telescopes to 1e-12 (measured 1.9e-16). SIPG symmetry probe (N=4 dense):
#      sym rel 1e-12 (measured 2.0e-16), min sym eig > 0.05 (measured 0.217), and the
#      RECEIVER-side-σ falsifier breaks symmetry > 1e-2 (measured 0.63) — the pairwise
#      pairing is load-bearing.
#   5. p-nonconforming variant: ladder (Nw,Ne) = (4,6), (6,8), (8,10) spectral (measured
#      8.23e-3 → 1.71e-4 → 1.95e-6) + the explicit N=5-vs-N=7 case (measured L2 3.05e-4,
#      interface ratio 0.99); same interface gates.
#   6. complex smoke: mortar flux vs hand formula on random complex data; real data stays
#      real-typed; 2:1 complex Helmholtz β = −k²−iε (k=5, ε=0.4), N=6: converged in 123
#      iterations (< 300 gate), L2 1.9e-4 (< 2e-3 gate).
#   7. BigFloat(256) projection exactness smoke (measured ≤ 1.5e-76 ≤ 1e-70 gates).
#
# Helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (same Knuth-MMIX LCG convention as sibling test files,
# ---- distinct name/seed) ---------------------------------------------------------------
mutable struct MORRNG
    s::UInt64
end
mor_unit!(r::MORRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                        Int64(r.s >> 11) / 2.0^53 * 2 - 1)
mor_vec!(r::MORRNG, n) = [mor_unit!(r) for _ in 1:n]
mor_cvec!(r::MORRNG, n) = [complex(mor_unit!(r), mor_unit!(r)) for _ in 1:n]

const MOR_SIDES = (-1, 1)
const MOR_C = 1.5                                   # CRP main.tex:595 penalty parameter

mor_one2(x, y) = 1.0
mor_zero2(x, y) = 0.0
mor_uex(x, y) = sinpi(x) * sinpi(y)                 # zero Dirichlet on [0,1]²
mor_S(x, y) = 2 * pi^2 * mor_uex(x, y)              # −∇²u = S

# least-squares slope of log10(err) against a resolution ladder
function mor_lsq_slope(Ns, errs)
    x = float.(collect(Ns))
    y = log10.(errs)
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

# map a mortar-logical point into the face coordinate (MORTAR_FULL/UPPER/LOWER)
mor_map(sz, x) = sz == MORTAR_FULL ? copy(x) :
                 sz == MORTAR_UPPER ? (x .+ 1) ./ 2 : (x .- 1) ./ 2

# ---- 2:1 h-nonconforming single-block mesh: W = [0,.5]×[0,1]; east half refined once in
# ---- BOTH dims (proper h-refinement) → 4 elements of 0.25×0.5. Mesh orchestration is
# ---- deliberately test-side (predecessor boundary: containers land with
# ---- matrix-free-operator); mortar objects come from the module under test. ------------
mor_blk21() = Block(71, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
mor_ids21() = (W  = ElementId(71, SegmentId(1, 0), SegmentId(0, 0)),
               E11 = ElementId(71, SegmentId(2, 2), SegmentId(1, 0)),
               E12 = ElementId(71, SegmentId(2, 2), SegmentId(1, 1)),
               E21 = ElementId(71, SegmentId(2, 3), SegmentId(1, 0)),
               E22 = ElementId(71, SegmentId(2, 3), SegmentId(1, 1)))

# conforming interior pairs (dim, lo-el, hi-el) of the 2:1 mesh
mor_conf21(t) = [(1, t.E11, t.E21), (1, t.E12, t.E22), (2, t.E11, t.E12), (2, t.E21, t.E22)]

# setup: element matrices/coefficients/offsets + mortar-coupled per-face fluxes.
# σ on the mortars is the PAIRWISE law from both sides' (p, h) — the node's W4-B handoff.
function mor_setup21(N, A, β; C = MOR_C, sigma_override = nothing)
    blk = mor_blk21()
    t = mor_ids21()
    ids = collect(t)
    ems = Dict(id => element_matrices(Float64, blk, id, (N, N)) for id in ids)
    cs = Dict(id => element_coefficients(ems[id], blk, id, (A, A), β,
                                         (mor_zero2, mor_zero2)) for id in ids)
    off = Dict(id => (k - 1) * N * N for (k, id) in enumerate(ids))
    fmW = face_matrices(ems[t.W], 1, +1)
    fmE = face_matrices(ems[t.E11], 1, -1)          # E12's west face is geometry-identical
    Am = ones(N)                                     # 𝔄 = A ≡ 1 at the mortar points
    mWlo = Mortar(fmW, fmE, MORTAR_LOWER, Am; C = C)
    mWup = Mortar(fmW, fmE, MORTAR_UPPER, Am; C = C)
    mE = Mortar(fmE, fmW, MORTAR_FULL, Am; C = C)
    if sigma_override !== nothing        # falsifier seam: RECEIVER-side σ — pairing each
        CW, CE = sigma_override          # side with its OWN face degenerates the pairwise
        mWlo = Mortar(fmW, fmW, MORTAR_LOWER, Am; C = CW)   # law to h_int (W4-B shortcut)
        mWup = Mortar(fmW, fmW, MORTAR_UPPER, Am; C = CW)
        mE = Mortar(fmE, fmE, MORTAR_FULL, Am; C = CE)
    end
    pen = PenaltyFlux(C)
    flux = Dict{ElementId,NumericalFlux}(
        t.W   => PerFaceFlux(((pen, MortarFlux([mWlo, mWup])), (pen, pen))),
        t.E11 => PerFaceFlux(((MortarFlux([mE]), pen), (pen, pen))),
        t.E12 => PerFaceFlux(((MortarFlux([mE]), pen), (pen, pen))),
        t.E21 => PerFaceFlux(((pen, pen), (pen, pen))),
        t.E22 => PerFaceFlux(((pen, pen), (pen, pen))))
    return t, ids, ems, cs, off, flux, (mWlo, mWup, mE)
end

# exterior FaceStates of element id (sender convention; mortar data pre-projected by the
# SENDER through project_to_mortar — "both neighbors project their face data to the mortar")
function mor_ext21(t, sts, id, mortars)
    mWlo, mWup, mE = mortars
    conf = mor_conf21(t)
    function confnb(d, sidx)
        for (dd, lo, hi) in conf
            dd == d || continue
            id == lo && sidx == 2 && return sts[hi].faces[d][1]
            id == hi && sidx == 1 && return sts[lo].faces[d][2]
        end
        return nothing
    end
    return [[begin
                 if id == t.W && d == 1 && sidx == 2
                     lo = project_to_mortar(mE.P, sts[t.E11].faces[1][1])
                     up = project_to_mortar(mE.P, sts[t.E12].faces[1][1])
                     FaceState(vcat(lo.psi, up.psi), vcat(lo.nflux, up.nflux))
                 elseif (id == t.E11 || id == t.E12) && d == 1 && sidx == 1
                     P = id == t.E11 ? mWlo.P : mWup.P
                     project_to_mortar(P, sts[t.W].faces[1][2])
                 else
                     nb = confnb(d, sidx)
                     nb === nothing ? zero_dirichlet_ghost(sts[id].faces[d][sidx]) :
                         exterior_state(nb)
                 end
             end for sidx in 1:2] for d in 1:2]
end

function mor_apply21(t, ids, ems, cs, off, flux, mortars, N)
    npe = N * N
    return function (y, u)
        sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+npe))
                   for id in ids)
        for id in ids
            y[off[id]+1:off[id]+npe] .= strong_residual(ems[id], cs[id], sts[id],
                                                        mor_ext21(t, sts, id, mortars),
                                                        flux[id])
        end
        return y
    end
end

function mor_rhs(ids, ems, off, Ns, S)
    b = zeros(ComplexF64, sum(N * N for N in values(Ns)))
    for id in ids
        em = ems[id]
        N = Ns[id]
        b[off[id]+1:off[id]+N*N] .=
            em.M * vec([S(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
    end
    return b
end

function mor_l2err(ids, ems, off, Ns, u, ue)
    s = 0.0
    for id in ids
        em = ems[id]
        N = Ns[id]
        e = u[off[id]+1:off[id]+N*N] .-
            vec([ue(em.x[1][i], em.x[2][j]) for i in 1:N, j in 1:N])
        s += real(dot(e, em.M * e))
    end
    return sqrt(s)
end

# max abs nodal error split: mortar-interface node layer vs everything else
function mor_iferr(ids, ems, off, Ns, u, ue, layerof)
    ifm, oth = 0.0, 0.0
    for id in ids
        em = ems[id]
        N = Ns[id]
        for j in 1:N, i in 1:N
            e = abs(u[off[id] + tensor_lin(i, j, N)] - ue(em.x[1][i], em.x[2][j]))
            layerof(id, i, N) ? (ifm = max(ifm, e)) : (oth = max(oth, e))
        end
    end
    return ifm, oth
end

# dense global matrix via unit-vector probes through the matrix-free apply
function mor_build_matrix(apply!, n)
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

@testset "mortars (dg-mortars)" begin

@testset "R = mass-conservative adjoint of P (Gauss masses + independent quadrature)" begin
    combos = [(4, 4, MORTAR_FULL), (4, 6, MORTAR_FULL), (5, 7, MORTAR_FULL),
              (5, 9, MORTAR_FULL), (4, 4, MORTAR_LOWER), (4, 4, MORTAR_UPPER),
              (3, 5, MORTAR_UPPER), (6, 8, MORTAR_LOWER), (6, 6, MORTAR_UPPER)]
    rng = MORRNG(0x00000000dead0071)
    for (Nf, Nm, sz) in combos
        P = projection_prolong(Float64, Nf, Nm, sz)
        R = projection_restrict(Float64, Nf, Nm, sz)
        @test size(P) == (Nm, Nf) && size(R) == (Nf, Nm)
        wf = EllipticCore.refelement(Float64, Nf).weights
        wm = EllipticCore.refelement(Float64, Nm).weights
        jac = sz == MORTAR_FULL ? 1.0 : 0.5
        # THE node gate: W_f·R == jac·Pᵀ·W_m (mass-conservative adjoint w.r.t. the Gauss
        # mass matrices; calibrated ≤ 2.8e-17)
        @test maximum(abs, Diagonal(wf) * R .- jac .* (transpose(P) * Diagonal(wm))) ≤ 1e-13
        # independent check with an EXTERNAL fine Gauss rule: the Gauss diagonal mass is
        # the exact L2 Gram (degree 2N−2 ≤ 2N−1; calibrated ≤ 7.8e-16), and the adjoint
        # identity holds as the integral statement ∫_face (Rg)·u = jac·∫_mortar g·(u∘map)
        # (calibrated ≤ 3.0e-15 relative)
        Nq = Nf + Nm + 3
        elq = EllipticCore.refelement(Float64, Nq)
        elf = EllipticCore.refelement(Float64, Nf)
        elm = EllipticCore.refelement(Float64, Nm)
        Lf = interpolation_matrix(elf, elq.nodes)
        @test maximum(abs, transpose(Lf) * Diagonal(elq.weights) * Lf .- Diagonal(wf)) ≤ 1e-14
        Lfm = interpolation_matrix(elf, mor_map(sz, elq.nodes))
        Lm = interpolation_matrix(elm, elq.nodes)
        g = mor_vec!(rng, Nm)
        u = mor_vec!(rng, Nf)
        D = sum(elq.weights .* (Lf * (R * g)) .* (Lf * u))
        Dm = sum(elq.weights .* (Lm * g) .* (Lfm * u))
        @test abs(D - jac * Dm) ≤ 1e-13 * max(abs(D), abs(jac * Dm), 1.0)
        # P is exact on every polynomial the face represents (degree ≤ Nf−1; ≤ 4.5e-16)
        tm = mor_map(sz, elm.nodes)
        for k in 0:Nf-1
            @test maximum(abs, P * (elf.nodes .^ k) .- tm .^ k) ≤ 1e-14
        end
    end
    # losslessness contract: a mortar coarser than the face is a construction error
    @test_throws ArgumentError projection_prolong(Float64, 5, 4, MORTAR_FULL)
    @test_throws ArgumentError projection_restrict(Float64, 5, 3, MORTAR_LOWER)
end

@testset "R∘P = I on the coarse space (full + half-pair sums)" begin
    for (Nf, Nm) in ((4, 4), (4, 6), (5, 7), (5, 9), (6, 8))
        P = projection_prolong(Float64, Nf, Nm, MORTAR_FULL)
        R = projection_restrict(Float64, Nf, Nm, MORTAR_FULL)
        @test maximum(abs, R * P - I) ≤ 1e-14                     # calibrated ≤ 6.7e-16
        # 2:1 split face: the two half-mortar round trips SUM to the identity
        Plo = projection_prolong(Float64, Nf, Nm, MORTAR_LOWER)
        Rlo = projection_restrict(Float64, Nf, Nm, MORTAR_LOWER)
        Pup = projection_prolong(Float64, Nf, Nm, MORTAR_UPPER)
        Rup = projection_restrict(Float64, Nf, Nm, MORTAR_UPPER)
        @test maximum(abs, Rlo * Plo + Rup * Pup - I) ≤ 1e-14     # calibrated ≤ 1.2e-15
    end
end

@testset "conforming pass-through: mortar exchange == direct exchange (bitwise)" begin
    N = 5
    # P and R are BITWISE identity matrices on a conforming (same-N, same-level) face
    Pid = projection_prolong(Float64, N, N, MORTAR_FULL)
    Rid = projection_restrict(Float64, N, N, MORTAR_FULL)
    @test Pid == Matrix{Float64}(I, N, N)
    @test Rid == Matrix{Float64}(I, N, N)

    blk = Block(72, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    id = ElementId(72, SegmentId(1, 0), SegmentId(0, 0))
    em = element_matrices(Float64, blk, id, (N, N))
    fm = face_matrices(em, 1, +1)
    rng = MORRNG(0x00000000dead0072)
    Af = mor_cvec!(rng, N)
    int = FaceState(mor_cvec!(rng, N), mor_cvec!(rng, N))
    ext = FaceState(mor_cvec!(rng, N), mor_cvec!(rng, N))
    # conforming Mortar: both sides ARE this face's geometry ⇒ σ is bitwise the
    # receiver-side PenaltyFlux σ, and the flux is bitwise PenaltyFlux
    m = Mortar(fm, fm, MORTAR_FULL, Af; C = MOR_C)
    ψm, φm = numerical_flux(MortarFlux([m]), fm, Af, int, ext)
    ψp, φp = numerical_flux(PenaltyFlux(MOR_C), fm, Af, int, ext)
    @test ψm == ψp && φm == φp
    # project_to_mortar through the identity P is the identity on the FaceState data
    fs2 = project_to_mortar(m.P, int)
    @test fs2.psi == int.psi && fs2.nflux == int.nflux

    # whole-mesh falsifier: mortar coupling on EVERY interior face of a conforming 2×2
    # mesh reproduces the direct exchange BITWISE through strong_residual
    N2 = 4
    ids = vec([ElementId(72, SegmentId(1, ix), SegmentId(1, iy)) for ix in 0:1, iy in 0:1])
    ems = Dict(i2 => element_matrices(Float64, blk, i2, (N2, N2)) for i2 in ids)
    cs = Dict(i2 => element_coefficients(ems[i2], blk, i2, (mor_one2, mor_one2), mor_zero2,
                                         (mor_zero2, mor_zero2)) for i2 in ids)
    off = Dict(i2 => (k - 1) * N2 * N2 for (k, i2) in enumerate(ids))
    pen = PenaltyFlux(MOR_C)
    mortartable(i2) = PerFaceFlux(ntuple(d -> ntuple(sidx -> begin
            nb = neighbor_id(i2, d, MOR_SIDES[sidx])
            if nb === nothing
                pen
            else
                fmi = face_matrices(ems[i2], d, MOR_SIDES[sidx])
                fme = face_matrices(ems[nb], d, -MOR_SIDES[sidx])
                MortarFlux([Mortar(fmi, fme, MORTAR_FULL, ones(N2); C = MOR_C)])
            end
        end, 2), 2))
    function mkapply(fluxof)
        return function (y, u)
            sts = Dict(i2 => element_state(ems[i2], cs[i2], view(u, off[i2]+1:off[i2]+N2*N2))
                       for i2 in ids)
            for i2 in ids
                ext2 = [[begin
                             nb = neighbor_id(i2, d, MOR_SIDES[sidx])
                             nb === nothing ?
                                 zero_dirichlet_ghost(sts[i2].faces[d][sidx]) :
                                 exterior_state(sts[nb].faces[d][sidx == 1 ? 2 : 1])
                         end for sidx in 1:2] for d in 1:2]
                y[off[i2]+1:off[i2]+N2*N2] .= strong_residual(ems[i2], cs[i2], sts[i2],
                                                              ext2, fluxof(i2))
            end
            return y
        end
    end
    rng2 = MORRNG(0x00000000dead0073)
    u = mor_cvec!(rng2, 4 * N2 * N2)
    y1 = zeros(ComplexF64, 4 * N2 * N2)
    y2 = zeros(ComplexF64, 4 * N2 * N2)
    mkapply(i2 -> pen)(y1, u)
    mkapply(mortartable)(y2, u)
    @test y1 == y2                                   # bitwise (calibrated true)
end

@testset "2:1 h-nonconforming MMS Poisson: spectral, no O(1) interface error" begin
    Ns = (4, 6, 8)
    errs = Float64[]
    iferrs = Float64[]
    for N in Ns
        t, ids, ems, cs, off, flux, mortars = mor_setup21(N, mor_one2, mor_zero2)
        # mortar geometry sanity: two half-mortars on the big face, full on the small
        @test size(mortars[1].P) == (N, N) && size(mortars[3].P) == (N, N)
        ap = mor_apply21(t, ids, ems, cs, off, flux, mortars, N)
        NsD = Dict(id => N for id in ids)
        b = mor_rhs(ids, ems, off, NsD, mor_S)
        res = gmres(ap, b; abstol = 1e-12, reltol = 0.0, maxiter = 2000, restart = 2000)
        @test res.converged
        @test res.iterations ≤ 500                   # calibrated 51/125/200
        push!(errs, mor_l2err(ids, ems, off, NsD, res.x, mor_uex))
        layer(id, i, Nl) = (id == t.W && i == Nl) || ((id == t.E11 || id == t.E12) && i == 1)
        fe, oe = mor_iferr(ids, ems, off, NsD, res.x, mor_uex, layer)
        push!(iferrs, fe)
        # NO O(1) interface error: the mortar-interface node layer is the same order as
        # the interior (calibrated ratio ≤ 1.42)
        @test fe ≤ 4 * oe
    end
    # spectral decay (calibrated 8.22e-3 → 1.72e-4 → 1.96e-6, ratios ≤ 0.021)
    @test errs[1] < 3e-2
    @test errs[2] / errs[1] ≤ 0.1
    @test errs[3] / errs[2] ≤ 0.1
    @test mor_lsq_slope(Ns, errs) ≤ -0.7             # calibrated −0.906/point
    @test errs[3] < 2e-5                             # NODE HEADLINE (calibrated 1.96e-6)
    # the interface layer decays at the same spectral rate (calibrated ratios ≤ 0.025)
    @test iferrs[2] / iferrs[1] ≤ 0.1
    @test iferrs[3] / iferrs[2] ≤ 0.1
end

@testset "2:1 conservation + SIPG symmetry + pairwise-σ falsifier" begin
    # conservation: Σ over the nonconforming interface of the lifted numerical flux
    # telescopes — ∫_Wface (nF)* + Σ_c ∫_cface (nF)* = 0 (mass-conservative adjoint +
    # exact antisymmetry of the sender-convention mortar algebra; calibrated 1.9e-16)
    N = 5
    t, ids, ems, cs, off, flux, mortars = mor_setup21(N, mor_one2, mor_zero2)
    rng = MORRNG(0x00000000dead0074)
    u = mor_cvec!(rng, length(ids) * N * N)
    sts = Dict(id => element_state(ems[id], cs[id], view(u, off[id]+1:off[id]+N*N))
               for id in ids)
    tot = 0.0 + 0im
    scale = 0.0
    for id in (t.W, t.E11, t.E12)
        d, sidx = 1, id == t.W ? 2 : 1
        fm = face_matrices(ems[id], d, MOR_SIDES[sidx])
        _, φstar = numerical_flux(flux[id].faces[d][sidx], fm, cs[id].Aface[d][sidx],
                                  sts[id].faces[d][sidx],
                                  mor_ext21(t, sts, id, mortars)[d][sidx])
        v = sum(fm.ML * φstar)                       # 1ᵀ(ML φ*) = ∫_face (nF)*
        tot += v
        scale = max(scale, abs(v))
    end
    @test abs(tot) ≤ 1e-12 * scale

    # SIPG symmetry probe on the 2:1 mesh (affine, 𝔄 = 1): machine-symmetric and
    # penalty-stabilized (calibrated sym 2.0e-16, min sym eig 0.217)
    N4 = 4
    t4, ids4, ems4, cs4, off4, flux4, mortars4 = mor_setup21(N4, mor_one2, mor_zero2)
    n = length(ids4) * N4 * N4
    A = mor_build_matrix(mor_apply21(t4, ids4, ems4, cs4, off4, flux4, mortars4, N4), n)
    @test maximum(abs, imag.(A)) == 0
    Ar = real.(A)
    @test opnorm(Ar - transpose(Ar), Inf) / opnorm(Ar, Inf) < 1e-12
    @test minimum(eigvals(Symmetric((Ar + transpose(Ar)) / 2))) > 0.05
    # FALSIFIER: receiver-side σ (each side using only its own h, the W4-B conforming
    # shortcut) breaks SIPG symmetry — the PAIRWISE mortar σ is load-bearing
    # (calibrated asymmetry 0.63)
    tb, idsb, emsb, csb, offb, fluxb, mortarsb =
        mor_setup21(N4, mor_one2, mor_zero2; sigma_override = (MOR_C, MOR_C))
    Ab = real.(mor_build_matrix(mor_apply21(tb, idsb, emsb, csb, offb, fluxb, mortarsb, N4), n))
    @test opnorm(Ab - transpose(Ab), Inf) / opnorm(Ab, Inf) > 1e-2
end

@testset "p-nonconforming MMS (incl. N=5 vs N=7): spectral, no O(1) interface error" begin
    blk = Block(73, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    PW = ElementId(73, SegmentId(1, 0), SegmentId(0, 0))
    PE = ElementId(73, SegmentId(1, 1), SegmentId(0, 0))
    function solvep(Nw, Ne)
        ems = Dict(PW => element_matrices(Float64, blk, PW, (Nw, Nw)),
                   PE => element_matrices(Float64, blk, PE, (Ne, Ne)))
        cs = Dict(id => element_coefficients(ems[id], blk, id, (mor_one2, mor_one2),
                                             mor_zero2, (mor_zero2, mor_zero2))
                  for id in (PW, PE))
        off = Dict(PW => 0, PE => Nw * Nw)
        NsD = Dict(PW => Nw, PE => Ne)
        n = Nw * Nw + Ne * Ne
        fmW = face_matrices(ems[PW], 1, +1)
        fmE = face_matrices(ems[PE], 1, -1)
        Nm = max(Nw, Ne)
        Am = ones(Nm)
        mW = Mortar(fmW, fmE, MORTAR_FULL, Am; C = MOR_C)   # mortar carries the LARGER N
        mE = Mortar(fmE, fmW, MORTAR_FULL, Am; C = MOR_C)
        @test size(mW.P) == (Nm, Nw) && size(mE.P) == (Nm, Ne)
        # both sides computed the SAME pairwise σ (max-p / min-h law)
        @test mW.sigma == mE.sigma
        @test mW.sigma == penalty_sigma(MOR_C, Nw - 1, Ne - 1,
                                        perpendicular_size(fmW), perpendicular_size(fmE))
        pen = PenaltyFlux(MOR_C)
        fluxW = PerFaceFlux(((pen, MortarFlux([mW])), (pen, pen)))
        fluxE = PerFaceFlux(((MortarFlux([mE]), pen), (pen, pen)))
        function ap(y, u)
            stW = element_state(ems[PW], cs[PW], view(u, 1:Nw*Nw))
            stE = element_state(ems[PE], cs[PE], view(u, Nw*Nw+1:n))
            gW(d, s) = zero_dirichlet_ghost(stW.faces[d][s])
            gE(d, s) = zero_dirichlet_ghost(stE.faces[d][s])
            extW = [[gW(1, 1), project_to_mortar(mE.P, stE.faces[1][1])],
                    [gW(2, 1), gW(2, 2)]]
            extE = [[project_to_mortar(mW.P, stW.faces[1][2]), gE(1, 2)],
                    [gE(2, 1), gE(2, 2)]]
            y[1:Nw*Nw] .= strong_residual(ems[PW], cs[PW], stW, extW, fluxW)
            y[Nw*Nw+1:n] .= strong_residual(ems[PE], cs[PE], stE, extE, fluxE)
            return y
        end
        b = mor_rhs((PW, PE), ems, off, NsD, mor_S)
        res = gmres(ap, b; abstol = 1e-12, reltol = 0.0, maxiter = 2000, restart = 2000)
        layer(id, i, Nl) = (id == PW && i == Nl) || (id == PE && i == 1)
        fe, oe = mor_iferr((PW, PE), ems, off, NsD, res.x, mor_uex, layer)
        return mor_l2err((PW, PE), ems, off, NsD, res.x, mor_uex), fe, oe,
               res.iterations, res.converged
    end
    # ladder (calibrated 8.23e-3 → 1.71e-4 → 1.95e-6, ratios ≤ 0.021; iters 34/70/111)
    perrs = Float64[]
    for (Nw, Ne) in ((4, 6), (6, 8), (8, 10))
        e, fe, oe, it, cv = solvep(Nw, Ne)
        @test cv
        @test it ≤ 500
        @test fe ≤ 4 * oe                            # calibrated ratio ≤ 1.28
        push!(perrs, e)
    end
    @test perrs[2] / perrs[1] ≤ 0.1
    @test perrs[3] / perrs[2] ≤ 0.1
    @test perrs[3] < 2e-5
    # the node contract's explicit case (calibrated L2 3.05e-4, interface ratio 0.99)
    e57, fe57, oe57, it57, cv57 = solvep(5, 7)
    @test cv57
    @test e57 < 3e-3
    @test fe57 ≤ 4 * oe57
end

@testset "complex smoke: mortar flux algebra + 2:1 Helmholtz" begin
    # hand-formula check on random complex data through a genuinely nonconforming mortar
    blk = Block(74, LinearMap1D(0.0, 1.0), LinearMap1D(0.0, 1.0))
    PW = ElementId(74, SegmentId(1, 0), SegmentId(0, 0))
    PE = ElementId(74, SegmentId(1, 1), SegmentId(0, 0))
    Nw, Ne = 4, 6
    emW = element_matrices(Float64, blk, PW, (Nw, Nw))
    emE = element_matrices(Float64, blk, PE, (Ne, Ne))
    fmW = face_matrices(emW, 1, +1)
    fmE = face_matrices(emE, 1, -1)
    rng = MORRNG(0x00000000dead0075)
    Am = mor_cvec!(rng, Ne)                          # complex 𝔄 at the mortar (m-mode case)
    mW = Mortar(fmW, fmE, MORTAR_FULL, Am; C = MOR_C)
    int = FaceState(mor_cvec!(rng, Nw), mor_cvec!(rng, Nw))
    ext = FaceState(mor_cvec!(rng, Ne), mor_cvec!(rng, Ne))   # already mortar data
    ψs, φs = numerical_flux(MortarFlux([mW]), fmW, ones(ComplexF64, Nw), int, ext)
    ψi = mW.P * int.psi
    φi = mW.P * int.nflux
    ψref = mW.R * ((ψi .+ ext.psi) ./ 2)
    φref = mW.R * ((φi .- ext.nflux) ./ 2 .- mW.sigma .* (Am .* (ψi .- ext.psi)))
    @test maximum(abs, ψs .- ψref) ≤ 1e-13 * max(maximum(abs, ψref), 1.0)
    @test maximum(abs, φs .- φref) ≤ 1e-13 * max(maximum(abs, φref), 1.0)
    @test length(ψs) == Nw && length(φs) == Nw       # restricted back to the RECEIVER face
    # real data stays real-typed through a real-coefficient mortar
    mR = Mortar(fmW, fmE, MORTAR_FULL, ones(Ne); C = MOR_C)
    ψr, φr = numerical_flux(MortarFlux([mR]), fmW, ones(Nw),
                            FaceState(mor_vec!(rng, Nw), mor_vec!(rng, Nw)),
                            FaceState(mor_vec!(rng, Ne), mor_vec!(rng, Ne)))
    @test eltype(ψr) == Float64 && eltype(φr) == Float64
    # exterior data not matching the mortar layout is a hard error
    @test_throws ArgumentError numerical_flux(MortarFlux([mW]), fmW, ones(ComplexF64, Nw),
                                              int, FaceState(mor_cvec!(rng, Ne + 1),
                                                             mor_cvec!(rng, Ne + 1)))

    # 2:1 complex Helmholtz β = −k² − iε (calibrated 123 iterations, L2 1.9e-4)
    k, ε = 5.0, 0.4
    βH = -k^2 - im * ε
    SH(x, y) = (2 * pi^2 + βH) * mor_uex(x, y)
    N = 6
    t, ids, ems, cs, off, flux, mortars = mor_setup21(N, mor_one2, (x, y) -> βH)
    ap = mor_apply21(t, ids, ems, cs, off, flux, mortars, N)
    NsD = Dict(id => N for id in ids)
    b = mor_rhs(ids, ems, off, NsD, SH)
    @test eltype(b) == ComplexF64
    res = gmres(ap, b; abstol = 1e-11, reltol = 0.0, maxiter = 2000, restart = 2000)
    @test res.converged
    @test res.iterations < 300
    @test mor_l2err(ids, ems, off, NsD, res.x, mor_uex) < 2e-3
end

@testset "BigFloat(256) projection exactness smoke" begin
    old = precision(BigFloat)
    try
        setprecision(BigFloat, 256)
        Nf, Nm = 5, 7
        P = projection_prolong(BigFloat, Nf, Nm, MORTAR_FULL)
        R = projection_restrict(BigFloat, Nf, Nm, MORTAR_FULL)
        @test eltype(P) == BigFloat && precision(eltype(P)) == 256
        xf = EllipticCore.refelement(BigFloat, Nf).nodes
        xm = EllipticCore.refelement(BigFloat, Nm).nodes
        for k in 0:Nf-1                              # calibrated ≤ 2.6e-77
            @test maximum(abs, P * (xf .^ k) .- xm .^ k) ≤ big"1e-70"
        end
        @test maximum(abs, R * P - I) ≤ big"1e-70"   # calibrated ≤ 1.5e-76
        Plo = projection_prolong(BigFloat, 4, 4, MORTAR_LOWER)
        Rlo = projection_restrict(BigFloat, 4, 4, MORTAR_LOWER)
        Pup = projection_prolong(BigFloat, 4, 4, MORTAR_UPPER)
        Rup = projection_restrict(BigFloat, 4, 4, MORTAR_UPPER)
        @test maximum(abs, Rlo * Plo + Rup * Pup - I) ≤ big"1e-70"   # calibrated ≤ 3.5e-77
    finally
        setprecision(BigFloat, old)
    end
end

end # testset mortars
