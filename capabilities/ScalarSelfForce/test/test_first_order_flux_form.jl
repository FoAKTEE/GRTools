# test_first_order_flux_form.jl — ESSF node `first-order-flux-form` (schedule W6-E,
# GLOBAL_DAG row 30, CRIT-join closure), written BEFORE the implementation (TDD). The node
# delivers the ASSEMBLED first-order flux form of arXiv:2606.04998 eq:flux/eq:first-order
# (main.tex:365-383) — `fluxform_callables(orbit, m, v::VtuSlicing, region)`: the coefficient
# CALLABLES (A = (A^rr, A^zz), beta, gamma = (γ^r, γ^z)) over the discretization chart
# q = (r, z) with the CONSTANT per-region vtu boost H ∈ {−1, 0, +1} (H′ ≡ 0, main.tex:316-330
# + 347), in exactly the shape EllipticCore's DG residual layer consumes
# (`element_coefficients(em, blk, id, fc.A, fc.beta, fc.gamma)`).
#
# PREDECESSOR SCOPE (sealed; cited, never re-derived here):
#   - elliptic-operator-tsliced: the coefficient ALGEBRA vs eq:triangle + replacement rules
#     (pointwise-FD twin of gate (c), 1.4e-16), β_m horizon finiteness at H = −1, conjugation,
#     and the structural pins A^rr(r₊) = 0 / A^zz(±1) = 0 AT the exact points;
#   - horizon-penetrating-compactification (W6-D): the pointwise bundle
#     `fluxform_coefficients` + r-chart ≡ r_*-chart closure + deep-horizon finiteness;
#   - vtu-slicing-jumps: the piecewise H instance, domain classification, junction DATA
#     (the [∂_{r_*}Ψ] = −imΩΨ hooks are node dg-jump-conditions, NOT here);
#   - hyperboloidal-transformation: r_*(r) + slicing kinematics.
# THIS node's remaining gap (row-30 oracles not yet gated):
#  (a) assembled-callables contract — bitwise the sealed pointwise composition at interior
#      points; the PER-REGION constant H ON the interface faces (where the pointwise vtu
#      boost takes the closed-t value — per-element classification, vtu.jl contract); the
#      t region == tsliced (0, 0) defaults; unknown region refused; RunConfig slicing chain.
#  (b) flux-vanishing closure (main.tex:385): F^r(r₊) = 0 EXACT through the callables (Δ
#      factored off cached horizons) — gated AT r₊ and at its float neighbors (one-ulp
#      probes bounded by eps(r₊)·(r₊−r₋)/(r₊²+a²)); F^z(z = ±1) = 0 exact + neighbors.
#      No explicit horizon/pole boundary conditions are ever needed.
#  (c) THE row-30 identity: −∂_i F^i + β Ψ + γ^i ∂_i Ψ == △_m Ψ on 10 random smooth complex
#      Ψ at 1e-10, per region H ∈ {−1, 0, +1}, via TEST-SIDE SPECTRAL differentiation —
#      Legendre-Gauss grids + differentiation matrices, the nested flux divergence
#      ∂_r(A^rr ∂_r Ψ) evaluated exactly as a nodal discretization evaluates it. Reference =
#      INDEPENDENT transcription of eq:triangle (main.tex:247-251) + replacement rules
#      (main.tex:306-308) + ∂_{r_*} = Δ/(r²+a²)∂_r (main.tex:356), raw Δ = r² − 2Mr + a²
#      (never the implementation's factored horizons), region-constant H from the printed
#      vtu boost table (main.tex:316-322). The pointwise-FD twin of this transcription is
#      the SEALED tsliced gate (a); this file re-transcribes the same PRINTED equations for
#      the spectral instrument (matrix form) — the assembled-callable path is what's new.
#  (d) DG readiness — EllipticCore element_coefficients/element_state/strong_residual consume
#      the callables directly on a v element whose upper radial face IS the vt interface
#      (𝔄 is evaluated AT that face, convention N2 — the point where per-region H matters).
#  (e) BigFloat(256) — precision inheritance through every callable + the spectral identity
#      at the BigFloat truncation floor (N = 48; coefficient-aliasing tail ρ^{-N} ≲ 1e-30 at
#      these box margins, far below the 1e-25 gate; 256-bit roundoff ~1e-70 is invisible).
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore: refelement, Block, LinearMap1D, ElementId, SegmentId,
                                    element_matrices, element_coefficients, element_state,
                                    strong_residual, FaceState, ElementCoefficients

# ---- deterministic pseudo-random complex polynomial fields (no Random dep) ---------------
# PCG-style 64-bit LCG (fixed constants) → uniform-ish complex coefficients in [−1,1)².
mutable struct Lcg64
    s::UInt64
end
_unit!(g::Lcg64) = ((g.s = 0x5851f42d4c957f2d * g.s + 0x14057b7ef767814f;
                     (g.s >> 11) * 2.0^-52 - 1.0))
_crand!(g::Lcg64) = complex(_unit!(g), _unit!(g))

# one random smooth complex field = bivariate polynomial of degree (DX, DY) in the LOGICAL
# box coordinates (x, y) ∈ [−1,1]² — smooth, exactly representable on N ≥ DX+1 Gauss nodes
const DX, DY = 8, 5
const FIELDS10 = let g = Lcg64(0x2606049980000001)
    [[_crand!(g) for _ in 0:DX, _ in 0:DY] for _ in 1:10]
end

vander(x::AbstractVector, d::Int) = [xi^k for xi in x, k in 0:d]
fieldgrid(C, x) = vander(x, DX) * C * transpose(vander(x, DY))   # U[i,j] = Ψ(x_i, y_j)

# ---- spectral instrument: Legendre-Gauss grid + differentiation on a physical box ---------
function specgrid(::Type{T}, N, ra, rb, za, zb) where {T}
    re = refelement(T, N)
    x = re.nodes
    r = @. (T(ra) + T(rb)) / 2 + (T(rb) - T(ra)) / 2 * x
    z = @. (T(za) + T(zb)) / 2 + (T(zb) - T(za)) / 2 * x
    (x = x, r = r, z = z,
     Dr = re.D .* (2 / (T(rb) - T(ra))), Dz = re.D .* (2 / (T(zb) - T(za))))
end

cgrid(f, r, z) = [f(ri, zj) for ri in r, zj in z]

# LHS: the assembled first-order form through the node's callables, ALL differentiation
# spectral — including the nested flux divergence ∂_i(A^ii ∂_i Ψ), exactly the structure a
# nodal DG discretization evaluates (collocation product then differentiation matrix)
function apply_first_order(fc, g, U)
    dUr = g.Dr * U
    dUz = U * transpose(g.Dz)
    -(g.Dr * (cgrid(fc.A[1], g.r, g.z) .* dUr)) .-
        ((cgrid(fc.A[2], g.r, g.z) .* dUz) * transpose(g.Dz)) .+
        cgrid(fc.beta, g.r, g.z) .* U .+
        cgrid(fc.gamma[1], g.r, g.z) .* dUr .+
        cgrid(fc.gamma[2], g.r, g.z) .* dUz
end

# INDEPENDENT reference: eq:triangle (main.tex:247-251) + replacement rules (main.tex:306-308)
# with the region-constant boost (H′ ≡ 0, main.tex:347) and ∂_{r_*} = Δ/(r²+a²)∂_r
# (main.tex:356). Raw Δ = r² − 2Mr + a², NOT the implementation's factored horizons.
# Transcribed AS PRINTED (m ≥ 0, the production branch; conjugation is the sealed tsliced
# gate (d)). Same spectral instrument; matrix twin of the sealed pointwise-FD transcription
# in test_operator_tsliced.jl.
function apply_triangle_spectral(orb, m, H, g, U)
    a, M, Om = orb.a, orb.M, orb.Omega
    rr, zt = g.r, transpose(g.z)
    Dl = @. rr^2 - 2M * rr + a^2
    w = @. rr^2 + a^2
    lam = Dl ./ w                                    # dr/dr_* — the chart binding
    dUr = g.Dr * U
    dUz = U * transpose(g.Dz)
    D1 = lam .* dUr                                  # ∂Ψ/∂r_*
    D2 = lam .* (g.Dr * (lam .* dUr))                # ∂²Ψ/∂r_*²
    D1h = D1 .+ (im * m * Om * H) .* U
    D2h = D2 .+ (2im * m * Om * H) .* D1 .- (m^2 * Om^2 * H^2) .* U    # H′ ≡ 0
    Zang = ((1 .- zt .^ 2) .* dUz) * transpose(g.Dz) .- 2m .* (zt .* dUz) .-
           m * (m + 1) .* U
    Sig = @. rr^2 + a^2 * zt^2
    @. -(w / Dl) * D2h + (2a^2 / (rr * w) - 2im * m * a / Dl) * D1h -
       (Zang + m^2 * Om * (2M * rr * (Om * w - 2a) + Om * Dl * Sig) / Dl * U +
        2 * (a^2 - im * m * a * rr - M * rr) / rr^2 * U) / w
end

# per-region physical boxes, strongly inside each slab with wide Bernstein-ellipse margins to
# the coefficient poles (r ∈ {0, r₋, r₊} and the complex w = 0 pair): the N = 32 spectral
# truncation of the rational coefficients sits many orders below the 1e-10 gate. Type-generic
# (BigFloat boxes from a BigFloat orbit/slicing).
function region_box(orb, v, region)
    region === :v ? (orb.rp + (v.r_v - orb.rp) / 2, v.r_v) :
    region === :t ? (orb.r0 - 3, orb.r0 + 3) : (v.r_u, v.r_u + 20)
end
const ZBOX = (-0.9, 0.7)                             # asymmetric: no equatorial specialness

# the printed vtu boost table (main.tex:316-322) — the reference's OWN region → H map
vtu_H_table(region) = region === :v ? -1 : region === :t ? 0 : 1

@testset "MModeOperator assembled first-order flux form (first-order-flux-form)" begin

    @testset "exports" begin
        @test Base.isexported(ScalarSelfForce, :fluxform_callables)
    end

    @testset "(a) assembled callables: bitwise composition + per-region H at the faces" begin
        for (a, r0) in ((0.5, 10.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            v = vtu_slicing(orb, 2.8, 24.4)
            hf = vtu_height_function(v)
            pts = Dict(:v => (orb.rp + 0.4, 0.33), :t => (10.7, -0.71), :u => (31.0, 0.52))
            for region in (:v, :t, :u), m in (0, 1, 5)
                fc = fluxform_callables(orb, m, v, region)
                r, z = pts[region]
                # strictly interior ⇒ the pointwise W6-D bundle (sealed) must be reproduced
                # BITWISE by the region callables
                c = fluxform_coefficients(orb, m, hf, r, z)
                @test fc.A[1](r, z) == c.Arr
                @test fc.A[2](r, z) == c.Azz
                @test fc.beta(r, z) == c.beta
                @test fc.gamma[1](r, z) == c.gamma_r
                @test fc.gamma[2](r, z) == c.gamma_z
            end
            # per-region H ON the interface faces (per-element classification, the vtu.jl
            # contract): the pointwise vtu boost takes the closed-t value 0 AT both
            # interfaces — a v/u element's callables must carry the slab constant there
            z, m = 0.44, 3
            @test slicing_kinematics(hf, orb, v.r_v).H == 0      # the pointwise trap
            @test slicing_kinematics(hf, orb, v.r_u).H == 0
            fcv = fluxform_callables(orb, m, v, :v)
            fcu = fluxform_callables(orb, m, v, :u)
            @test fcv.beta(v.r_v, z) == beta_m(orb, m, v.r_v, z, -1, 0)
            @test fcv.gamma[1](v.r_v, z) == gamma_r_m(orb, m, v.r_v, -1)
            @test fcu.beta(v.r_u, z) == beta_m(orb, m, v.r_u, z, 1, 0)
            @test fcu.gamma[1](v.r_u, z) == gamma_r_m(orb, m, v.r_u, 1)
            @test fcv.beta(v.r_v, z) != beta_m(orb, m, v.r_v, z, 0, 0)   # not the t value
            # t-region callables == the tsliced (H, H′) = (0, 0) defaults, bitwise
            fct = fluxform_callables(orb, m, v, :t)
            @test fct.beta(9.1, z) == beta_m(orb, m, 9.1, z)
            @test fct.gamma[1](9.1, z) == gamma_r_m(orb, m, 9.1)
            # unknown region refused
            @test_throws ArgumentError fluxform_callables(orb, m, v, :w)
        end
        # RunConfig slicing chain: the production transition data feeds the constructor
        orb = CircularOrbit(0.5, 10.0)
        cfg = crp_paper_config(spin = 0.5)
        vc = vtu_slicing(orb, cfg)
        fc = fluxform_callables(orb, cfg.orbit.m_mode_number, vc, :t)
        @test fc.beta(10.3, 0.2) == beta_m(orb, cfg.orbit.m_mode_number, 10.3, 0.2)
    end

    @testset "(b) flux vanishing: F^r(r₊) = 0 exact + neighbors; F^z(±1) = 0 exact + neighbors" begin
        dpsi = 1.3 - 0.7im       # any finite ∂Ψ — the flux zero is the coefficient's
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            rv = a == 0.998 ? 2.0 : 2.8                  # r₊ < r_v < r0 for every spin
            v = vtu_slicing(orb, rv, r0 + 14.4)
            for region in (:v, :t, :u)
                fc = fluxform_callables(orb, 2, v, region)
                for z in (-0.52, 0.37)
                    # F^r(r₊) = 0 EXACTLY (main.tex:385): Δ factored off the cached horizons
                    # (tsliced pin) ⇒ the callable is zero AT r₊, hence the flux with ANY
                    # finite ∂_rΨ — no explicit horizon boundary condition
                    @test fc.A[1](orb.rp, z) == 0
                    @test iszero(fc.A[1](orb.rp, z) * dpsi)
                    # float neighbors: one-ulp probes bounded by eps(r₊)·(r₊−r₋)/(r₊²+a²)
                    @test abs(fc.A[1](nextfloat(orb.rp), z)) <= 4 * eps(orb.rp)
                    @test abs(fc.A[1](prevfloat(orb.rp), z)) <= 4 * eps(orb.rp)
                end
                for r in (orb.rp + 0.4, r0, 88.0)
                    # F^z(z = ±1) = 0 EXACTLY: pole-accurate product form (1−z)(1+z)
                    @test fc.A[2](r, 1.0) == 0
                    @test fc.A[2](r, -1.0) == 0
                    @test iszero(fc.A[2](r, -1.0) * dpsi)
                    @test abs(fc.A[2](r, prevfloat(1.0))) <= 4 * eps()
                    @test abs(fc.A[2](r, nextfloat(-1.0))) <= 4 * eps()
                end
            end
        end
        # BigFloat horizon: exact zero at the 256-bit cached horizon too
        setprecision(BigFloat, 256) do
            orbB = CircularOrbit(big"0.5", big"10.0")
            vB = vtu_slicing(orbB, big"2.8", big"24.4")
            fcB = fluxform_callables(orbB, 2, vB, :v)
            @test fcB.A[1](orbB.rp, big"0.3") == 0
            @test fcB.A[2](big"7.7", big"1.0") == 0
        end
    end

    @testset "(c) −∂F + βΨ + γ∂Ψ == △_m Ψ: spectral, 10 random Ψ, per region H (1e-10)" begin
        N = 32
        worst = 0.0
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            v = vtu_slicing(orb, 2.8, 24.4)
            for region in (:v, :t, :u)
                ra, rb = region_box(orb, v, region)
                g = specgrid(Float64, N, ra, rb, ZBOX...)
                H = vtu_H_table(region)
                for m in (0, 2, 5)
                    fc = fluxform_callables(orb, m, v, region)
                    for C in FIELDS10
                        U = fieldgrid(C, g.x)
                        lhs = apply_first_order(fc, g, U)
                        ref = apply_triangle_spectral(orb, m, H, g, U)
                        dev = maximum(abs.(lhs .- ref)) / max(1.0, maximum(abs.(ref)))
                        worst = max(worst, dev)
                        @test dev <= 1e-10
                    end
                end
            end
        end
        println("first-order-flux-form gate (c) worst spectral identity rel dev: ", worst)
    end

    @testset "(d) DG readiness: element_coefficients/element_state/strong_residual" begin
        orb = CircularOrbit(0.5, 10.0)
        v = vtu_slicing(orb, 2.8, 24.4)
        m = 3
        fc = fluxform_callables(orb, m, v, :v)
        # innermost v element: physical block [r₊ + 0.4, r_v] × z ∈ [0, 0.5]; its UPPER radial
        # face IS the vt interface — element_coefficients evaluates 𝔄 there (convention N2),
        # exactly where the per-region H contract matters
        blk = Block(7, LinearMap1D(orb.rp + 0.4, v.r_v), LinearMap1D(0.0, 0.5))
        id = ElementId(7, SegmentId(0, 0), SegmentId(0, 0))
        em = element_matrices(Float64, blk, id, (6, 5))
        c = element_coefficients(em, blk, id, fc.A, fc.beta, fc.gamma)
        @test c isa ElementCoefficients{ComplexF64}
        x1, x2 = em.x
        @test c.A[1] == vec([fc.A[1](ξ, η) for ξ in x1, η in x2])
        @test c.A[2] == vec([fc.A[2](ξ, η) for ξ in x1, η in x2])
        @test c.beta == vec([fc.beta(ξ, η) for ξ in x1, η in x2])
        @test c.gamma[1] == vec([fc.gamma[1](ξ, η) for ξ in x1, η in x2])
        @test c.gamma[2] == vec([fc.gamma[2](ξ, η) for ξ in x1, η in x2])
        # 𝔄^{rr} at the vt-interface face points: finite, positive (r_v > r₊)
        @test c.Aface[1][2] == [fc.A[1](v.r_v, η) for η in x2]
        @test all(>(0), real.(c.Aface[1][2]))
        # apply: element_state + N7 transmissive exterior ⇒ finite complex residual
        ψ = vec([complex(sin(0.7ξ) + 0.2η, 0.1ξ * η) for ξ in x1, η in x2])
        st = element_state(em, c, ψ)
        ext = ntuple(d -> ntuple(s -> FaceState(st.faces[d][s].psi,
                                                -st.faces[d][s].nflux), 2), 2)
        R = strong_residual(em, c, st, ext)
        @test eltype(R) == ComplexF64
        @test all(isfinite, R)
    end

    @testset "(e) BigFloat(256): precision inheritance + spectral identity (1e-25)" begin
        setprecision(BigFloat, 256) do
            orb = CircularOrbit(big"0.5", big"10.0")
            v = vtu_slicing(orb, big"2.8", big"24.4")
            for region in (:v, :t, :u)
                fc = fluxform_callables(orb, 3, v, region)
                b = fc.beta(big"7.1", big"0.44")
                @test b isa Complex{BigFloat} && precision(real(b)) == 256
                @test fc.A[1](big"7.1", big"0.44") isa BigFloat
                @test fc.A[2](big"7.1", big"0.44") isa BigFloat
                @test fc.gamma[1](big"7.1", big"0.44") isa Complex{BigFloat}
                @test fc.gamma[2](big"7.1", big"0.44") isa BigFloat
            end
            # spectral identity at the BigFloat truncation floor on the two BOOSTED slabs
            # (N = 48: coefficient-aliasing tail ρ^{-N} ≲ 1e-30 at these box margins)
            N, m = 48, 3
            for region in (:v, :u)
                ra, rb = region_box(orb, v, region)
                g = specgrid(BigFloat, N, ra, rb, ZBOX...)
                fc = fluxform_callables(orb, m, v, region)
                U = fieldgrid(FIELDS10[1], g.x)
                lhs = apply_first_order(fc, g, U)
                ref = apply_triangle_spectral(orb, m, vtu_H_table(region), g, U)
                dev = Float64(maximum(abs.(lhs .- ref)) /
                              max(one(BigFloat), maximum(abs.(ref))))
                @test dev <= 1e-25
            end
        end
    end
end
