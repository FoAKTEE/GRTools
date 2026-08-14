# test_fluxform.jl — ESSF node `horizon-penetrating-compactification` (schedule W6-D), written
# BEFORE the implementation (TDD). The node is the CLOSURE of the horizon-penetrating
# compactification (arXiv:2606.04998 Sec. II.B.3, main.tex:354-360): the DG solver
# discretizes and differentiates in Boyer–Lindquist r, NOT r_*, and the binding of the
# operator to that DISCRETIZATION chart is ∂_{r_*} = Δ/(r²+a²) ∂_r. The r-chart coefficient
# substance was already delivered by predecessors (tsliced.jl coefficients are r-chart
# functions with general (H, H′); hyperboloidal.jl owns r_*(r); vtu.jl owns the piecewise H);
# this node adds ONLY the thin chart-converted coefficient bundle `fluxform_coefficients`
# (the pointwise product the domain/solve nodes consume at quadrature nodes) and the
# horizon-penetration closure gates that no predecessor ran.
#
# Gates (node contract, GLOBAL_DAG row 29):
#  (a) r-form ≡ r_*-form on smooth fields, 1e-10 — the SAME physical operator applied two
#      ways at the SAME physical points, away from the horizon: the r-chart side uses the
#      bundle with ALL radial differentiation performed numerically ON THE r CHART
#      (central FD in r of physical samples, including the nested divergence ∂_r(A^rr ∂_r Ψ)
#      exactly as a discretization sees it); the r_*-chart side uses an INDEPENDENT
#      transcription of the printed r_*-form radial operator (eq:triangle radial part +
#      replacement rules, main.tex:247-251 + 306-308: principal coefficient w/Δ, first-order
#      coefficient 2a²/(rw) − 2ima/Δ − 2imΩH·w/Δ, raw Δ = r² − 2Mr + a²) with ALL radial
#      differentiation ON THE r_* CHART (the field is DEFINED as a function of r_*). β and
#      the angular part are chart-shared and enter both sides identically — their algebra
#      is the sealed predecessor claim (ESSF::elliptic-operator-tsliced gate (a)), NOT
#      re-verified here; what this gate adds is precisely the chart-conversion content.
#  (b) horizon-penetration closure — ALL five bundle coefficients finite (no NaN/Inf) at
#      r = r₊(1 + 1e-12) with the vtu boost H = −1 (the DAG row-29 acceptance point), and
#      down the whole collocation-reachable approach: r = r₊(1 + 10^−k), k = 2..12, plus the
#      Legendre-Gauss nodes (N = 32) of an innermost element [r₊, r_v]. β accuracy pinned
#      against the SAME-IEEE-POINT BigFloat(256) bundle (pure O(eps/Δ) Float64 cancellation
#      noise of the 2Mr/w + H term, measured ≲ 1e-4 at k = 12). The smooth H = −1 limit
#      curve (main.tex:386) is gated with a FULL-BIGFLOAT orbit — the implementation
#      factors Δ off the orbit's CACHED horizons, so a Float64 orbit's β carries a residual
#      pole at fl(r₊) with residue ~eps (probing offsets ≲ 1e-16 from the true horizon is
#      outside its faithful domain; the mixed-precision boundary is separately bounded at
#      the k = 12 acceptance point, measured ~1e-4 ≤ 1e-2). Discriminator: t-slicing
#      (H = 0) blows up ~1/Δ at the same point (gate not vacuous).
#  (c) u-domain interplay smoke — all five coefficients finite AND smooth (spectral tail of
#      the Legendre expansion ≤ 1e-8 relative) across N = 24 Gauss samples of the wavezone
#      1/r-compaction map (EllipticCore InverseRadialMap1D, eq:compaction main.tex:546-549,
#      paper values r_u = 24.4, R = 1e4) with the vtu boost H = +1 — the coefficients are
#      analytic in 1/r at infinity, so composed with the affine-in-1/r map they must be
#      spectrally resolvable in ξ; this is the property the u-domain discretization needs.
#  (d) a = 0 degeneracy — the bundle collapses to an independent test-side Schwarzschild
#      transcription (Δ = r² − 2Mr, all a-terms dropped) at rtol 1e-13, incl. H′ ≠ 0.
#  (e) BigFloat(256) smoke — full-BigFloat orbit + slicing ⇒ Complex{BigFloat}/BigFloat
#      bundle at precision 256; chart identity (a) at the BigFloat noise floor (gate 1e-20).
#  (f) thin-layer contract — the bundle is BITWISE the slicing_kinematics ∘ tsliced.jl
#      composition (no algebraic drift is possible), and the t_slicing() bundle reproduces
#      the tsliced.jl (H, H′) = (0, 0) defaults bitwise.
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore: InverseRadialMap1D, map_forward
using ScalarSelfForce.NumericsCore: gl_nodes_weights

# ---- BigFloat central FD ----------------------------------------------------------------
cdiff(f, x, h) = (f(x + h) - f(x - h)) / (2h)

# ---- smooth complex test fields DEFINED ON THE r_* CHART ---------------------------------
# Ψ(r_*, z) = P(z)·g₁(r_*) + Q(z)·g₂(r_*): complex z-polynomials times two Gaussians in r_*
# (analytic everywhere on the real r_* line, decaying both ways). Analytic FIRST derivatives
# only (division-consistent d/ds to keep BigFloat residuals algebraically clean); the second
# r_* derivative in gate (a) enters through cdiff of drs — differentiation ON the r_* chart.
function rsfield(pc, qc, c1, s1, c2, s2)
    dpc = ntuple(i -> i * pc[i + 1], length(pc) - 1)
    dqc = ntuple(i -> i * qc[i + 1], length(qc) - 1)
    g1(s) = exp(-((s - c1) / s1)^2)
    dg1(s) = -2 * ((s - c1) / s1) / s1 * g1(s)
    g2(s) = exp(-((s - c2) / s2)^2)
    dg2(s) = -2 * ((s - c2) / s2) / s2 * g2(s)
    (f   = (s, z) -> evalpoly(z, pc) * g1(s) + evalpoly(z, qc) * g2(s),
     drs = (s, z) -> evalpoly(z, pc) * dg1(s) + evalpoly(z, qc) * dg2(s),
     dz  = (s, z) -> evalpoly(z, dpc) * g1(s) + evalpoly(z, dqc) * g2(s))
end

# deterministic complex coefficients (no Random dep), r_* centers/widths sized for the
# sample range r_* ∈ ~[−4, 45] of the gate-(a) points
const FFFIELDS = (
    rsfield((complex(0.83, -0.21), complex(0.40, 0.90), complex(-0.37, 0.62),
             complex(0.155, 0.201)),
            (complex(-0.058, -0.113), complex(0.91, -0.44), complex(0.21, 0.077)),
            12.0, 9.0, -2.0, 14.0),
    rsfield((complex(-0.50, 0.70), complex(0.31, -0.11), complex(0.09, 0.44),
             complex(-0.23, 0.05)),
            (complex(0.66, 0.13), complex(-0.40, 0.52), complex(0.12, -0.31)),
            4.0, 6.5, 20.0, 11.0),
)

# ---- height functions ---------------------------------------------------------------------
# smooth tanh boost with a CONSISTENT (h, H, H′) triple (exact antiderivative/derivative):
# h = log(cosh(k(r_*−c)))/k, H = tanh(k(r_*−c)), H′ = k/cosh²(k(r_*−c)). Type-generic.
tanh_boost(k, c) = HeightFunction(
    s -> log(cosh(k * (s - c))) / k,
    s -> tanh(k * (s - c)),
    s -> k / cosh(k * (s - c))^2,
)

# ---- gate (a) operator applications --------------------------------------------------------
# Chart-shared part (identical calls on both sides): β·Ψ + the angular divergence + γ^z ∂_z Ψ,
# with the angular divergence by cdiff in z of A^zz(r, ·)·∂_zΨ (A^zz from the bundle).
function apply_shared(orb, m, hf, fld, r, z, rs, h)
    c = fluxform_coefficients(orb, m, hf, r, z)
    -cdiff(x -> fluxform_coefficients(orb, m, hf, r, x).Azz * fld.dz(rs, x), z, h) +
        c.beta * fld.f(rs, z) + c.gamma_z * fld.dz(rs, z)
end

# r-chart radial part: ALL radial differentiation is central FD IN r of PHYSICAL samples
# g(x) = Ψ(r_*(x), z) — including the nested flux divergence ∂_r(A^rr ∂_r g), exactly the
# structure a nodal discretization in r evaluates. Uses ONLY the exported bundle + rstar.
function apply_radial_r_chart(orb, m, hf, fld, r, z, h)
    g(x) = fld.f(rstar(orb, x), z)
    drg(x) = cdiff(g, x, h)
    F(x) = fluxform_coefficients(orb, m, hf, x, z).Arr * drg(x)
    c = fluxform_coefficients(orb, m, hf, r, z)
    -cdiff(F, r, h) + c.gamma_r * drg(r)
end

# r_*-chart radial part: INDEPENDENT transcription of the printed r_*-form radial operator —
# eq:triangle radial part (main.tex:247-251) with the replacement-rule boost terms collected
# on the ∂_{r_*} slots (main.tex:306-308; the Ψ-proportional remainder is β's, chart-shared):
#   −(w/Δ)·∂²_{r_*}Ψ + (2a²/(rw) − 2ima/Δ − 2imΩH·w/Δ)·∂_{r_*}Ψ ,  Δ = r² − 2Mr + a² RAW
# (not the implementation's factored horizons). All radial differentiation ON the r_* chart:
# analytic ∂_{r_*}Ψ + cdiff of it in r_* steps. H from the exported kinematics.
function apply_radial_rstar_chart(orb, m, hf, fld, r, z, h)
    a, M, Om = orb.a, orb.M, orb.Omega
    D = r^2 - 2M * r + a^2
    w = r^2 + a^2
    rs = rstar(orb, r)
    H = slicing_kinematics(hf, orb, r).H
    d2rs = cdiff(s -> fld.drs(s, z), rs, h)
    -(w / D) * d2rs +
        (2a^2 / (r * w) - 2im * m * a / D - 2im * m * Om * H * w / D) * fld.drs(rs, z)
end

# ---- gate (c) helper: Legendre spectral tail of Gauss samples ------------------------------
legP(l, x) = l == 0 ? one(x) :
             l == 1 ? x :
             begin
                 pm, p = one(x), x
                 for n in 1:(l - 1)
                     pm, p = p, ((2n + 1) * x * p - n * pm) / (n + 1)
                 end
                 p
             end
function legendre_tail_ratio(vals, xs, ws)
    n = length(xs)
    cs = [abs((2l + 1) / 2 * sum(ws[i] * legP(l, xs[i]) * vals[i] for i in 1:n))
          for l in 0:(n - 1)]
    maximum(cs[(end - 2):end]) / max(maximum(cs), 1e-300)
end

@testset "MModeOperator horizon-penetrating chart bundle (horizon-penetrating-compactification)" begin

    @testset "exports" begin
        @test Base.isexported(ScalarSelfForce, :fluxform_coefficients)
    end

    @testset "(f) thin-layer contract: bitwise composition + t-slicing defaults" begin
        for (a, r0) in ((0.5, 10.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            vtu = vtu_slicing(orb, orb.rp + 1.0, r0 + 5.0)
            for hf in (t_slicing(), tanh_boost(0.35, 8.0), vtu_height_function(vtu)),
                (r, z, m) in ((orb.rp + 0.6, 0.33, 2), (7.3, -0.71, 5), (r0 + 11.0, 0.52, 1))

                c = fluxform_coefficients(orb, m, hf, r, z)
                kin = slicing_kinematics(hf, orb, r)
                @test c.Arr == flux_coeff_r(orb, r)
                @test c.Azz == flux_coeff_z(orb, r, z)
                @test c.beta == beta_m(orb, m, r, z, kin.H, kin.Hprime)
                @test c.gamma_r == gamma_r_m(orb, m, r, kin.H)
                @test c.gamma_z == gamma_z_m(orb, m, r, z)
            end
            # t_slicing() bundle == the tsliced.jl (H, H′) = (0, 0) defaults, bitwise
            ct = fluxform_coefficients(orb, 2, t_slicing(), 7.0, 0.62)
            @test ct.beta == beta_m(orb, 2, 7.0, 0.62)
            @test ct.gamma_r == gamma_r_m(orb, 2, 7.0)
        end
    end

    @testset "(a) r-form == r_*-form on smooth fields (1e-10, per-chart differentiation)" begin
        worst = 0.0
        setprecision(BigFloat, 256) do
            h = big"1e-13"   # nested-FD noise floor: smooth inner truncation ~h² = 1e-26
            for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
                orb = CircularOrbit(a, r0)
                vtu = vtu_slicing(orb, orb.rp + 1.0, r0 + 5.0)
                # one point per vtu domain (v | t | u), all away from the horizon
                pts = ((orb.rp + 0.6, 0.33), ((orb.rp + 1.0 + r0 + 5.0) / 2, -0.71),
                       (r0 + 11.0, 0.52))
                for hf in (t_slicing(), tanh_boost(0.35, 8.0), vtu_height_function(vtu)),
                    m in (0, 1, 2, 5), fld in FFFIELDS, (rr, zz) in pts

                    r, z = BigFloat(rr), BigFloat(zz)
                    rs = rstar(orb, r)
                    shared = apply_shared(orb, m, hf, fld, r, z, rs, h)
                    lhs = apply_radial_r_chart(orb, m, hf, fld, r, z, h) + shared
                    ref = apply_radial_rstar_chart(orb, m, hf, fld, r, z, h) + shared
                    dev = Float64(abs(lhs - ref) / max(one(BigFloat), abs(ref)))
                    worst = max(worst, dev)
                    @test dev <= 1e-10
                end
            end
        end
        println("fluxform gate (a) worst r-chart vs r*-chart rel dev: ", worst)
    end

    @testset "(b) horizon-penetration closure: all coefficients finite at r₊(1+1e-12), H = -1" begin
        worst_f64 = 0.0
        worst_lim = 0.0
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            r_v = a == 0.998 ? 2.0 : 2.8          # keep r₊ < r_v < r0 for every spin
            vtu = vtu_slicing(orb, r_v, r0 + 14.4)
            hf = vtu_height_function(vtu)
            z = 0.37
            # full-BigFloat twin: horizon cache exact at working precision, so the H = −1
            # limit curve of the COEFFICIENT LAYER is probed without the fl(r₊) residual
            # pole of a Float64 orbit's factored Δ (see the header note)
            orbB, vtuB = setprecision(BigFloat, 256) do
                oB = CircularOrbit(BigFloat(a), BigFloat(r0))
                oB, vtu_slicing(oB, BigFloat(r_v), BigFloat(r0 + 14.4))
            end
            hfB = vtu_height_function(vtuB)
            for m in (1, 2, 5)
                # BigFloat(256) H = −1 horizon limit reference (offset 1e-30 off the
                # BigFloat-orbit horizon: cancellation noise eps₂₅₆/1e-30 ~ 1e-47,
                # genuinely ON the smooth limit curve)
                blim = setprecision(BigFloat, 256) do
                    fluxform_coefficients(orbB, m, hfB, orbB.rp * (1 + big"1e-30"),
                                          BigFloat(z)).beta
                end
                # the k = 12 acceptance point ON the consistent-precision limit curve
                devlim = setprecision(BigFloat, 256) do
                    b12 = fluxform_coefficients(orbB, m, hfB, orbB.rp * (1 + big"1e-12"),
                                                BigFloat(z)).beta
                    Float64(abs(b12 - blim) / (1 + abs(blim)))
                end
                worst_lim = max(worst_lim, devlim)
                @test devlim <= 1e-9
                # geometric approach: the DAG row-29 acceptance point k = 12 and shallower
                for k in 2:12
                    r = orb.rp * (1 + 10.0^-k)
                    c = fluxform_coefficients(orb, m, hf, r, z)
                    @test slicing_kinematics(hf, orb, r).H == -1
                    @test isfinite(c.Arr) && isfinite(c.Azz) && isfinite(c.beta) &&
                          isfinite(c.gamma_r) && isfinite(c.gamma_z)
                    # β against the SAME-IEEE-POINT BigFloat bundle (BigFloat(r) is exact;
                    # same Float64 orbit, so the factored-Δ horizons agree bitwise): the
                    # Float64 evaluation carries only the bounded O(eps/Δ) cancellation
                    # noise of the 2Mr/w + H term
                    bref = setprecision(BigFloat, 256) do
                        fluxform_coefficients(orb, m, hf, BigFloat(r), BigFloat(z)).beta
                    end
                    dev = Float64(abs(c.beta - bref) / (1 + abs(bref)))
                    worst_f64 = max(worst_f64, dev)
                    @test dev <= 3e-2
                    # mixed-precision boundary at the acceptance point: the Float64-orbit
                    # value sits within the fl(r₊) residual-pole bound ~eps/offset ~ 1e-4
                    # of the true (BigFloat-orbit) limit
                    if k == 12
                        @test Float64(abs(bref - blim) / (1 + abs(blim))) <= 1e-2
                    end
                end
                # collocation-reachable approach: Legendre-Gauss nodes (N = 32) of an
                # innermost element [r₊, r_v] — strictly interior, exactly what a DG
                # discretization can ever sample (main.tex:556)
                xi, _ = gl_nodes_weights(Float64, 32)
                for x in xi
                    r = orb.rp + (r_v - orb.rp) * (1 + x) / 2
                    c = fluxform_coefficients(orb, m, hf, r, z)
                    @test isfinite(c.Arr) && isfinite(c.Azz) && isfinite(c.beta) &&
                          isfinite(c.gamma_r) && isfinite(c.gamma_z)
                end
                # discriminator: WITHOUT horizon penetration of the slicing (t-chart H = 0)
                # β blows up ~1/Δ at the acceptance point — the H = −1 gate is not vacuous
                @test abs(fluxform_coefficients(orb, m, t_slicing(),
                                                orb.rp * (1 + 1e-12), z).beta) > 1e6
            end
        end
        println("fluxform gate (b) worst Float64-vs-BigFloat beta rel dev: ", worst_f64,
                "; worst k=12 limit-curve dev: ", worst_lim)
    end

    @testset "(c) u-domain smoke: finite + spectrally smooth across the 1/r-compaction map, H = +1" begin
        orb = CircularOrbit(0.5, 10.0)
        vtu = vtu_slicing(orb, 2.8, 24.4)              # the crp production transitions
        hf = vtu_height_function(vtu)
        cmap = InverseRadialMap1D(24.4, 1e4)           # eq:compaction paper values
        xi, wq = gl_nodes_weights(Float64, 24)         # strictly interior ⇒ r > r_u ⇒ H = +1
        z = 0.4
        worst = 0.0
        for m in (0, 2, 5)
            rs_pts = [map_forward(cmap, x) for x in xi]
            @test all(r -> slicing_kinematics(hf, orb, r).H == 1, rs_pts)
            cs = [fluxform_coefficients(orb, m, hf, r, z) for r in rs_pts]
            for name in (:Arr, :Azz, :beta, :gamma_r, :gamma_z)
                vals = [getfield(c, name) for c in cs]
                @test all(v -> isfinite(v), vals)
                tail = legendre_tail_ratio(vals, xi, wq)
                worst = max(worst, tail)
                @test tail <= 1e-8
            end
        end
        println("fluxform gate (c) worst compaction-map spectral tail ratio: ", worst)
    end

    @testset "(d) a = 0 degeneracy: independent Schwarzschild transcription (rtol 1e-13)" begin
        orb = CircularOrbit(0.0, 10.0)
        vtu = vtu_slicing(orb, 3.5, 15.0)
        Om = orb.Omega
        for hf in (t_slicing(), tanh_boost(0.35, 8.0), vtu_height_function(vtu)),
            m in (0, 1, 5), (r, z) in ((2.6, 0.33), (7.3, -0.71), (26.0, 0.52))

            kin = slicing_kinematics(hf, orb, r)
            H, Hp = kin.H, kin.Hprime
            am = abs(m)
            c = fluxform_coefficients(orb, m, hf, r, z)
            @test isapprox(c.Arr, 1 - 2 / r; rtol = 1e-13)
            @test isapprox(c.Azz, (1 - z^2) / r^2; rtol = 1e-13)
            @test isapprox(c.beta,
                           (am * (am + 1) + 2 / r) / r^2 +
                           (r^2 / (r^2 - 2r)) * (m^2 * Om^2 * (H^2 - 1) - im * m * Om * Hp);
                           rtol = 1e-13, atol = 1e-15)
            @test isapprox(c.gamma_r, -2im * m * Om * H; rtol = 1e-13, atol = 1e-15)
            @test isapprox(c.gamma_z, 2 * am * z / r^2; rtol = 1e-13, atol = 1e-15)
        end
    end

    @testset "(e) BigFloat(256) smoke: precision inheritance + chart identity" begin
        setprecision(BigFloat, 256) do
            orb = CircularOrbit(big"0.5", big"10.0")
            @test orb isa CircularOrbit{BigFloat}
            vtu = vtu_slicing(orb, big"2.8", big"24.4")
            @test vtu isa VtuSlicing{BigFloat}
            for hf in (tanh_boost(0.35, 8.0), vtu_height_function(vtu))
                c = fluxform_coefficients(orb, 3, hf, big"7.1", big"0.44")
                @test c.Arr isa BigFloat && precision(c.Arr) == 256
                @test c.Azz isa BigFloat
                @test c.beta isa Complex{BigFloat} && precision(real(c.beta)) == 256
                @test c.gamma_r isa Complex{BigFloat}
                @test c.gamma_z isa BigFloat
            end
            # chart identity at the BigFloat noise floor (nested-FD smooth truncation ~h²)
            h = big"1e-13"
            hf = tanh_boost(0.35, 8.0)
            for fld in FFFIELDS, m in (0, 3)
                r, z = big"7.1", big"0.44"
                rs = rstar(orb, r)
                shared = apply_shared(orb, m, hf, fld, r, z, rs, h)
                lhs = apply_radial_r_chart(orb, m, hf, fld, r, z, h) + shared
                ref = apply_radial_rstar_chart(orb, m, hf, fld, r, z, h) + shared
                @test abs(lhs - ref) <= big"1e-20" * max(one(BigFloat), abs(ref))
            end
            # deep-horizon BigFloat: the acceptance point at full precision, H = −1
            hfv = vtu_height_function(vtu)
            cdeep = fluxform_coefficients(orb, 2, hfv, orb.rp * (1 + big"1e-12"), big"0.37")
            @test isfinite(cdeep.beta) && isfinite(cdeep.gamma_r)
        end
    end
end
