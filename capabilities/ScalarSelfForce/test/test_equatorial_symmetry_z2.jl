# test_equatorial_symmetry_z2.jl — ESSF node `equatorial-symmetry-z2` (schedule W6-F,
# GLOBAL_DAG row 31, closure), written BEFORE the implementation (TDD). The node delivers
# the z² equatorial-symmetry chart of the ASSEMBLED first-order flux form
# (arXiv:2606.04998 main.tex:389-394): `fluxform_callables_z2(orbit, m, v, region)` — the
# coefficient callables over the HALF-domain chart q = (r, w), w = z² ∈ [0, 1], with the
# printed chain-rule factors 𝓕^{z²} = 4z²𝓕^z and γ^{z²} = 2zγ^z + 2(1−z²)/(r²+a²), in the
# same bundle shape EllipticCore's DG layer consumes — selectable via the RunConfig
# equatorial-symmetry flag (`Background.CircularOrbit.ImposeEquatorialSymmetry`, yaml:36):
# `fluxform_callables(orbit, config::RunConfig, region)`. The production lattice is (r, z²)
# (`block_bounds_z2 = [0, dz², 1]`), so this bundle is what the domain/solve nodes evaluate
# at quadrature nodes when the flag is on.
#
# PREDECESSOR SCOPE (sealed; cited, never re-derived here):
#   - elliptic-operator-tsliced: the z² coefficient FUNCTIONS flux_coeff_z2/gamma_z2_m with
#     factors read from main.tex:389-394, gated pointwise (w-FD divergence) vs the z form at
#     4.3e-42, incl. A^{z²} = 4z²A^zz bitwise + zeros at w ∈ {0, 1} AT the points;
#   - first-order-flux-form (W6-E): the assembled z-chart `fluxform_callables` — region → H
#     slab map, bitwise-composition, flux-vanishing and the spectral operator identity vs an
#     independent eq:triangle transcription (5.28e-14); that sealed bundle is THE reference
#     of gate (a) below;
#   - runtime-config-fixture: RunConfig/OrbitConfig incl. impose_equatorial_symmetry (yaml:36);
#   - vtu-slicing-jumps: VtuSlicing + the RunConfig → vtu_slicing chain.
# THIS node's remaining gap (row-31 oracles not yet gated):
#  (a) THE chain-rule chart identity at the ASSEMBLED level: the (r, w = z²)-chart operator
#      applied to even-z fields == the (r, z)-chart operator at 1e-12, with TEST-SIDE
#      SPECTRAL differentiation in BOTH charts (Legendre-Gauss grids + differentiation
#      matrices per chart, nested collocation flux divergence exactly as a nodal DG evaluates
#      it; the two nodal results compared at common physical points via barycentric
#      interpolation in the angular direction — exact here, the angular dependence of every
#      coefficient is polynomial). Even-z fields are polynomials in w = z², so ONE field
#      evaluates on both grids; z-chart nodes take BOTH signs of z (evenness is part of the
#      identity, γ^z is odd in z). INSTRUMENT ARITHMETIC IS EXACT (BigFloat 256, the sealed
#      tsliced-gate-(c) precedent) while the callables under test are evaluated in PRODUCTION
#      Float64 — the two charts' differentiation matrices share no roundoff (unlike the
#      same-grid W6-E comparison), and a pure-Float64 instrument floors at ~3.5e-12 of its
#      own uncorrelated matmul noise, above the gate; exact instrument arithmetic makes the
#      1e-12 gate measure the shipped callables' chart identity (measured ~6e-15).
#  (b) polar-boundary flux: F^{z²} vanishes at w = 1 EXACTLY through the callable (+ float
#      neighbors), and at the equator w = 0 EXACTLY too — the even sector needs no explicit
#      boundary data at either end of the half-domain (main.tex:385 + 399-400 via the chart).
#  (c) printed factors bitwise (main.tex:390-393) + equator regularity: no 1/z = w^{−1/2}
#      term survives for even fields — every angular dependence is EXACTLY the printed
#      polynomial in w (A^{z²} = 4w(1−w)/(r²+a²), γ^{z²} = (4|m|w + 2(1−w))/(r²+a²), β affine
#      in w), finite at w = 0.
#  (d) RunConfig selection: flag ON (the yaml:36 production default) ⇒ the z² bundle with
#      (m, slicing) from the config; flag OFF ⇒ the sealed z-chart bundle unchanged; unknown
#      region refused through both entries.
#  (e) BigFloat(256): precision inheritance through every callable + gate (a) at N = 48 on
#      both boosted slabs at 1e-30 (the r-direction truncation is IDENTICAL in the two charts
#      and cancels — the residual is pure roundoff, invisible at 256 bits).
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore: refelement

# ---- deterministic pseudo-random complex polynomial fields (no Random dep) ---------------
# PCG-style 64-bit LCG (fixed constants, distinct seed from the sibling W6-E instrument).
mutable struct Z2Lcg
    s::UInt64
end
z2_unit!(g::Z2Lcg) = ((g.s = 0x5851f42d4c957f2d * g.s + 0x14057b7ef767814f;
                       (g.s >> 11) * 2.0^-52 - 1.0))
z2_crand!(g::Z2Lcg) = complex(z2_unit!(g), z2_unit!(g))

# one random smooth EVEN field = bivariate polynomial of degree (DX, DW) in (logical-r,
# logical-w): as a function of z it is a polynomial in z² — even by construction, exactly
# representable on either chart's Gauss grid (degree 2·DW in z ≤ N − 1)
const Z2DX, Z2DW = 8, 5
const Z2FIELDS = let g = Z2Lcg(0x2606049980000002)
    [[z2_crand!(g) for _ in 0:Z2DX, _ in 0:Z2DW] for _ in 1:10]
end

z2_vander(x, d) = [xi^k for xi in x, k in 0:d]
# evaluate the field on a chart grid: pass that chart's logical-r and logical-w coordinates
z2_field(C, xr, yw) = z2_vander(xr, Z2DX) * C * transpose(z2_vander(yw, Z2DW))

# ---- spectral instrument: 1D Legendre-Gauss pieces per chart ------------------------------
function z2_grid1d(::Type{T}, N, qa, qb) where {T}
    re = refelement(T, N)
    (x = re.nodes,
     q = (T(qa) + T(qb)) / 2 .+ (T(qb) - T(qa)) / 2 .* re.nodes,
     D = re.D .* (2 / (T(qb) - T(qa))))
end

z2_logical(q, qa, qb) = @. (2q - (qa + qb)) / (qb - qa)

# barycentric interpolation matrix: nodal values at x ↦ values at query points xq
# (exact on polynomials of degree ≤ length(x) − 1 — gated below before first use)
function z2_interpmat(x::AbstractVector, xq)
    bw = [inv(prod(x[j] - x[k] for k in eachindex(x) if k != j)) for j in eachindex(x)]
    P = [bw[j] / (q - x[j]) for q in xq, j in eachindex(x)]
    P ./ sum(P, dims = 2)
end

z2_cgrid(f, r, q) = [f(ri, qj) for ri in r, qj in q]

# assembled-bundle application, chart-agnostic: nested collocation flux divergence
# (collocation product then differentiation matrix), exactly the sealed W6-E instrument
function z2_apply(fc, gr, gq, U)
    dUr = gr.D * U
    dUq = U * transpose(gq.D)
    -(gr.D * (z2_cgrid(fc.A[1], gr.q, gq.q) .* dUr)) .-
        ((z2_cgrid(fc.A[2], gr.q, gq.q) .* dUq) * transpose(gq.D)) .+
        z2_cgrid(fc.beta, gr.q, gq.q) .* U .+
        z2_cgrid(fc.gamma[1], gr.q, gq.q) .* dUr .+
        z2_cgrid(fc.gamma[2], gr.q, gq.q) .* dUq
end

# gate-(a) twin of z2_apply: same instrument, but the callables under test are evaluated in
# PRODUCTION Float64 (at the Float64 image of each node) inside exact BigFloat grids/matrices
z2_cgrid64(f, r, q) = [f(Float64(ri), Float64(qj)) for ri in r, qj in q]
function z2_apply64(fc, gr, gq, U)
    dUr = gr.D * U
    dUq = U * transpose(gq.D)
    -(gr.D * (z2_cgrid64(fc.A[1], gr.q, gq.q) .* dUr)) .-
        ((z2_cgrid64(fc.A[2], gr.q, gq.q) .* dUq) * transpose(gq.D)) .+
        z2_cgrid64(fc.beta, gr.q, gq.q) .* U .+
        z2_cgrid64(fc.gamma[1], gr.q, gq.q) .* dUr .+
        z2_cgrid64(fc.gamma[2], gr.q, gq.q) .* dUq
end

# per-region radial boxes — the sealed W6-E choice (wide Bernstein-ellipse margins to the
# coefficient poles; the angular direction has NO poles, both charts are polynomial there)
z2_region_box(orb, v, region) =
    region === :v ? (orb.rp + (v.r_v - orb.rp) / 2, v.r_v) :
    region === :t ? (orb.r0 - 3, orb.r0 + 3) : (v.r_u, v.r_u + 20)

const Z2_ZBOX = (-0.85, 0.9)   # z chart: asymmetric, BOTH signs of z (no equatorial specialness)
const Z2_WBOX = (0.0, 0.81)    # w chart: covers every z-node image w = z² (0.9² = 0.81)

@testset "MModeOperator z² equatorial-symmetry chart (equatorial-symmetry-z2)" begin

    @testset "exports" begin
        @test Base.isexported(ScalarSelfForce, :fluxform_callables_z2)
    end

    @testset "(a) chain-rule chart identity: (r, w)-operator == (r, z)-operator, 1e-12" begin
        # exact-arithmetic instrument, production Float64 callables (see the file header):
        # grids, differentiation matrices, fields and the interpolation are BigFloat(256);
        # the ONLY Float64-rounded quantities are the shipped callables' values
        setprecision(BigFloat, 256) do
            N = 32
            gw = z2_grid1d(BigFloat, N, Z2_WBOX...)
            gz = z2_grid1d(BigFloat, N, Z2_ZBOX...)
            P = z2_interpmat(gw.q, gz.q .^ 2)
            # instrument self-check: the interpolation is exact on a degree-7 polynomial in
            # w, so a gate failure below indicts the bundle, never the interpolant
            pcoef = [0.3 + 0.7k for k in 0:7]
            pw = [sum(pcoef[k+1] * wj^k for k in 0:7) for wj in gw.q]
            pz = [sum(pcoef[k+1] * zk^(2k) for k in 0:7) for zk in gz.q]
            @test maximum(abs.(P * pw .- pz)) <= 1e-50
            yw_w = gw.x
            yw_z = z2_logical(gz.q .^ 2, Z2_WBOX...)
            worst = zero(BigFloat)
            for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (-0.7, 10.0))
                orb = CircularOrbit(a, r0)
                v = vtu_slicing(orb, 2.8, 24.4)
                for region in (:v, :t, :u)
                    ra, rb = z2_region_box(orb, v, region)
                    gr = z2_grid1d(BigFloat, N, ra, rb)
                    for m in (0, 2, 5)
                        fcz = fluxform_callables(orb, m, v, region)    # sealed z-chart ref
                        fc2 = fluxform_callables_z2(orb, m, v, region) # the node
                        for C in Z2FIELDS
                            Uw = z2_field(C, gr.x, yw_w)  # ψ(r, w) at the (r, w) nodes
                            Uz = z2_field(C, gr.x, yw_z)  # SAME even field at the (r, z) nodes
                            Lz = z2_apply64(fcz, gr, gz, Uz)
                            dev = maximum(abs.(z2_apply64(fc2, gr, gw, Uw) * transpose(P) .-
                                               Lz)) / max(one(BigFloat), maximum(abs.(Lz)))
                            worst = max(worst, dev)
                            @test dev <= 1e-12
                        end
                    end
                end
            end
            println("equatorial-symmetry-z2 gate (a) worst chart-identity rel dev: ",
                    Float64(worst))
        end
    end

    @testset "(b) F^{z²} = 0 at the pole w = 1 AND the equator w = 0: exact + neighbors" begin
        dpsi = 1.3 - 0.7im     # any finite ∂_wΨ — the flux zero is the coefficient's
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            rv = a == 0.998 ? 2.0 : 2.8                  # r₊ < r_v < r0 for every spin
            v = vtu_slicing(orb, rv, r0 + 14.4)
            for region in (:v, :t, :u)
                fc2 = fluxform_callables_z2(orb, 2, v, region)
                for r in (orb.rp + 0.4, r0, 88.0)
                    # pole w = 1 (main.tex:385 through the chart): exact zero coefficient
                    @test fc2.A[2](r, 1.0) == 0
                    @test iszero(fc2.A[2](r, 1.0) * dpsi)
                    @test abs(fc2.A[2](r, prevfloat(1.0))) <= 8 * eps()
                    # equator w = 0: the even half-domain needs no boundary data either
                    @test fc2.A[2](r, 0.0) == 0
                    @test iszero(fc2.A[2](r, 0.0) * dpsi)
                    @test abs(fc2.A[2](r, nextfloat(0.0))) <= 8 * eps()
                end
                # the radial sector passes through the chart untouched: F^r(r₊) = 0 exact
                # (sealed W6-E gate, re-pinned through the z² bundle)
                @test fc2.A[1](orb.rp, 0.37) == 0
            end
        end
    end

    @testset "(c) printed factors (main.tex:390-393) bitwise + equator regularity (no 1/z)" begin
        orb = CircularOrbit(0.5, 10.0)
        v = vtu_slicing(orb, 2.8, 24.4)
        for region in (:v, :t, :u), m in (0, 3)
            fcz = fluxform_callables(orb, m, v, region)
            fc2 = fluxform_callables_z2(orb, m, v, region)
            for (r, w) in ((5.2, 0.04), (10.0, 0.36), (30.0, 0.7225))
                z = sqrt(w)
                s = r^2 + orb.a^2
                # 𝓕^{z²} = 4z²𝓕^z at coefficient level (main.tex:390-391) — bitwise
                @test fc2.A[2](r, w) == 4z^2 * fcz.A[2](r, z)
                # γ^{z²} = 2zγ^z + 2(1−z²)/(r²+a²) (main.tex:392-393) — bitwise
                @test fc2.gamma[2](r, w) == 2z * fcz.gamma[2](r, z) +
                                            2 * (1 - z) * (1 + z) / s
                # the r sector and β are the z-chart values at z = √w — bitwise
                @test fc2.A[1](r, w) == fcz.A[1](r, z)
                @test fc2.beta(r, w) == fcz.beta(r, z)
                @test fc2.gamma[1](r, w) == fcz.gamma[1](r, z)
                # equator regularity: no 1/z = w^{−1/2} survives in the even sector — every
                # angular dependence is EXACTLY the printed polynomial in w (a w^{−1/2} term
                # could not agree with these to ~eps at generic w, and diverges at w = 0)
                @test abs(fc2.A[2](r, w) - 4 * w * (1 - w) / s) <= 8 * eps()
                @test abs(fc2.gamma[2](r, w) - (4 * abs(m) * w + 2 * (1 - w)) / s) <=
                      8 * eps() * (1 + 4m)
                b0, b1 = fc2.beta(r, 0.0), fc2.beta(r, 1.0)   # β affine in w (1 − z² = 1 − w)
                @test abs(fc2.beta(r, w) - (b0 + (b1 - b0) * w)) <= 1e-13 * max(1.0, abs(b0))
            end
            # finite AT the equator (the m ≥ 1 sector of the z chart has γ^z ∝ z odd — the
            # even-field chart absorbs it into the regular γ^{z²}): exact values at w = 0
            for r in (5.2, 10.0, 30.0)
                @test fc2.gamma[2](r, 0.0) == 2 / (r^2 + orb.a^2)
                @test isfinite(fc2.beta(r, 0.0))
                @test isfinite(fc2.gamma[2](r, 1.0))
            end
        end
    end

    @testset "(d) RunConfig equatorial-symmetry flag selects the chart" begin
        orb = CircularOrbit(0.5, 10.0)
        cfg = crp_paper_config(spin = 0.5)
        @test cfg.orbit.impose_equatorial_symmetry === true    # yaml:36 production default
        vc = vtu_slicing(orb, cfg)
        m = cfg.orbit.m_mode_number
        fc = fluxform_callables(orb, cfg, :t)
        ref2 = fluxform_callables_z2(orb, m, vc, :t)
        for (r, q) in ((10.3, 0.36), (5.2, 0.04))
            @test fc.A[1](r, q) == ref2.A[1](r, q)
            @test fc.A[2](r, q) == ref2.A[2](r, q)
            @test fc.beta(r, q) == ref2.beta(r, q)
            @test fc.gamma[1](r, q) == ref2.gamma[1](r, q)
            @test fc.gamma[2](r, q) == ref2.gamma[2](r, q)
        end
        # the second coordinate IS w = z², not z — and m comes from the config (m_mode = 2)
        @test fc.A[2](10.3, 0.36) == flux_coeff_z2(orb, 10.3, sqrt(0.36))
        @test fc.A[2](10.3, 0.36) != flux_coeff_z(orb, 10.3, 0.36)
        @test fc.beta(10.3, 0.36) == beta_m(orb, m, 10.3, sqrt(0.36))
        # flag OFF ⇒ the sealed z-chart bundle, unchanged
        oc = cfg.orbit
        oc_off = ScalarSelfForce.OrbitConfig{Float64}(
            oc.black_hole_mass, oc.black_hole_spin, oc.orbital_radius, oc.m_mode_number,
            oc.hyperboloidal_slicing_transitions, oc.penetrating_horizon, false)
        cfg_off = RunConfig{Float64}(
            oc_off, cfg.null_slicing_blocks, cfg.initial_guess_is_background,
            cfg.randomize_initial_guess, cfg.domain, cfg.amr, cfg.dg, cfg.observers,
            cfg.gmres, cfg.multigrid, cfg.schwarz, cfg.build_matrix_subfile_name,
            cfg.build_matrix_verbosity, cfg.enable_direct_solve, cfg.element_distribution,
            cfg.avoid_global_proc0, cfg.singletons)
        fcoff = fluxform_callables(orb, cfg_off, :t)
        refz = fluxform_callables(orb, m, vc, :t)
        for (r, z) in ((10.3, 0.6), (5.2, -0.3))
            @test fcoff.A[2](r, z) == refz.A[2](r, z)
            @test fcoff.beta(r, z) == refz.beta(r, z)
            @test fcoff.gamma[1](r, z) == refz.gamma[1](r, z)
            @test fcoff.gamma[2](r, z) == refz.gamma[2](r, z)
        end
        # unknown region refused through both entries (the sealed region → H map validates)
        @test_throws ArgumentError fluxform_callables(orb, cfg, :x)
        @test_throws ArgumentError fluxform_callables_z2(orb, 2, vc, :w)
    end

    @testset "(e) BigFloat(256): precision inheritance + chart identity (1e-30)" begin
        setprecision(BigFloat, 256) do
            orbB = CircularOrbit(big"0.5", big"10.0")
            vB = vtu_slicing(orbB, big"2.8", big"24.4")
            for region in (:v, :t, :u)
                fc2 = fluxform_callables_z2(orbB, 3, vB, region)
                rB, wB = big"7.1", big"0.44"
                @test fc2.A[1](rB, wB) isa BigFloat
                a2 = fc2.A[2](rB, wB)
                @test a2 isa BigFloat && precision(a2) == 256
                b = fc2.beta(rB, wB)
                @test b isa Complex{BigFloat} && precision(real(b)) == 256
                @test fc2.gamma[1](rB, wB) isa Complex{BigFloat}
                g2 = fc2.gamma[2](rB, wB)
                @test g2 isa BigFloat && precision(g2) == 256
                # exact zeros at the 256-bit pole AND equator
                @test fc2.A[2](rB, big"1.0") == 0
                @test fc2.A[2](rB, big"0.0") == 0
            end
            # chart identity on the two boosted slabs at N = 48: the r-direction truncation
            # is IDENTICAL in both charts (the radial sector never touches the angular chart)
            # and cancels; the angular direction is exactly polynomial in each chart — the
            # residual is pure roundoff, ≲ 1e-60 at 256 bits (gate 1e-30 with huge margin)
            N = 48
            gw = z2_grid1d(BigFloat, N, Z2_WBOX...)
            gz = z2_grid1d(BigFloat, N, Z2_ZBOX...)
            P = z2_interpmat(gw.q, gz.q .^ 2)
            yw_w = gw.x
            yw_z = z2_logical(gz.q .^ 2, big(Z2_WBOX[1]), big(Z2_WBOX[2]))
            for region in (:v, :u)
                ra, rb = z2_region_box(orbB, vB, region)
                gr = z2_grid1d(BigFloat, N, ra, rb)
                fcz = fluxform_callables(orbB, 3, vB, region)
                fc2 = fluxform_callables_z2(orbB, 3, vB, region)
                Uw = z2_field(Z2FIELDS[1], gr.x, yw_w)
                Uz = z2_field(Z2FIELDS[1], gr.x, yw_z)
                Lz = z2_apply(fcz, gr, gz, Uz)
                dev = Float64(maximum(abs.(z2_apply(fc2, gr, gw, Uw) * transpose(P) .- Lz)) /
                              max(one(BigFloat), maximum(abs.(Lz))))
                @test dev <= 1e-30
            end
        end
    end
end
