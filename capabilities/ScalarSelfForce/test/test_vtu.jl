# test_vtu.jl — ESSF node `vtu-slicing-jumps` (schedule W5-D, GLOBAL_DAG row 28), written
# BEFORE the implementation (TDD). The node delivers the vtu piecewise slicing instance
# (arXiv:2606.04998 Sec. II.B.2, main.tex:315-347): the piecewise-constant boost
# H(r_*) ∈ {−1, 0, +1} with transition radii (r^v_*, r^u_*), the continuous piecewise
# height function h(r_*), v/t/u domain classification of a Boyer–Lindquist radius, and the
# two junction conditions [∂_{r_*}Ψ_m] = −imΩ Ψ_m (eq:jump_vt main.tex:337-339, eq:jump_tu
# main.tex:341-343) as pure interface DATA for the DG jump hooks.
# Predecessor scope (NOT re-delivered here, only composed against): the HeightFunction
# carrier + rstar (hyperboloidal.jl, node hyperboloidal-transformation), the phase map
# Ψ_m = Ψ̃_m e^{−imΩh} (ansatz.jl psi_from_tilde, node mode-ansatz-decomposition), the
# (H, H′) coefficient slots incl. flux_coeff_r (tsliced.jl, node elliptic-operator-tsliced),
# and the yaml-twin RunConfig transition radii (config.jl, node runtime-config-fixture).
#
# Gates (node contract):
#  (a) construction + tortoise consistency — VtuSlicing carries the transition radii in
#      BOTH charts, rstar_v/rstar_u bitwise = rstar(orbit, r_v/r_u); the paper constraint
#      r^v_* < r_{*0} < r^u_* (main.tex:323) is enforced (worldline in the t domain), as
#      are r_v < r_u and r_v > r₊ (tortoise defined); violations throw ArgumentError.
#  (b) RunConfig round-trip — vtu_slicing reads crp_paper_config()'s
#      HyperboloidalSlicingTransitions (v→t, v→t, t→u, t→u) EXACTLY (the IEEE doubles the
#      rendered yaml carries), the transitions coincide with the radial block bounds
#      [2]/[5] (the DG double-sampling premise: interfaces are element faces), the paper
#      ordering r^v_* < r_{*0} < r^u_* holds, a second Jinja2 parameter point round-trips,
#      and orbit/config mismatches + non-degenerate (smooth-band) transitions throw.
#  (c) piecewise boost + height exactness (main.tex:316-330) — H = −1/0/+1 exact in the
#      open v/t/u domains and H(±far) = ±1 (hyperboloidal condition, main.tex:278);
#      h = (r_*^v − r_*) / 0 / (r_* − r_*^u) bitwise per domain, h(interfaces) = 0 exact
#      (continuity anchors) with prevfloat/nextfloat continuity at the ulp scale;
#      H′ ≡ 0 exact EVERYWHERE pointwise (the interface deltas live in the junction data,
#      main.tex:347), FD of h reproduces H and FD of H vanishes inside open domains;
#      evaluations are type-generic (BigFloat in ⇒ BigFloat out).
#  (d) domain classification — vtu_domain: :v strictly below r_v, :t on the closed
#      [r_v, r_u] (interface convention pinned; Legendre-Gauss volume nodes never sample
#      a face, so element classification only ever queries interior points), :u strictly
#      above; consistent with the boost value at rstar(r); the crp_paper_config block
#      columns classify v|t|t|t|u — matching the yaml NullSlicingBlocks (Block0/5 v-sliced,
#      Block4/9 u-sliced).
#  (e) JUNCTION CONDITIONS (the row's oracle) — per-domain analytic branches
#      Ψ_dom(r_*) = Ψ̃(r_*)·e^{−imΩ h_dom(r_*)} of a smooth t-sliced Ψ̃ (test-side
#      transcription of the main.tex:325-330 height branches through the exported
#      psi_from_tilde): Ψ continuous across both interfaces (1e-14; branch values agree AT
#      the interface and the global-height Ψ matches across prevfloat/nextfloat), and the
#      derivative jump [∂_{r_*}Ψ] = (∂ from above) − (∂ from below) equals −imΩ Ψ|_iface
#      to 1e-12 at BOTH interfaces — each one-sided limit computed by BigFloat(256)
#      central FD (h = 1e-25) of the branch that coincides with Ψ on that side; the sign
#      orientation (increasing r_*, eq:jump_vt/eq:jump_tu LHS order) is pinned by the
#      flipped-sign target failing; the VtuJunction DATA (coeff = −imΩ at r/rstar)
#      reproduces the measured jump; and [F^r] = coeff·Ψ (F^r = A^rr ∂_rΨ IS ∂_{r_*}Ψ)
#      binds the data to the flux language the DG FaceState hooks consume.
#  (f) m = 0 degeneracy — junction coefficients exactly zero, the phase map is the
#      identity (Ψ ≡ Ψ̃), and the measured FD jump vanishes.
#  (g) BigFloat inheritance — full-BigFloat orbit ⇒ VtuSlicing{BigFloat} with
#      precision-256 tortoise radii, BigFloat height/boost evaluations, Complex{BigFloat}
#      junction coefficients.
#
# All local names are vsj_-prefixed: sibling test files share the runtests.jl include
# scope, and this file must also run standalone as the node's verification command.
using Test
using ScalarSelfForce

# ---- BigFloat central FD ------------------------------------------------------------------
vsj_cdiff(f, x, h) = (f(x + h) - f(x - h)) / (2h)

# ---- smooth complex t-sliced test field Ψ̃(r_*): entire, O(1), nonzero at interfaces ------
# low-order complex polynomial × unit-modulus oscillation; type-generic (BigFloat inherited).
vsj_tld(rs) = (complex(0.83, -0.29) + complex(0.31, 0.57) * (rs / 40) +
               complex(-0.12, 0.23) * (rs / 40)^2) * cis(rs / 17)

# ---- per-domain height branches: test-side transcription of main.tex:325-330 --------------
# each branch is the domain's linear height law extended smoothly (entire), so central FD of
# the branch at the interface IS the one-sided derivative limit from that domain's side.
vsj_branch_h(dom, rsv, rsu) =
    dom === :v ? (rs -> rsv - rs) :
    dom === :u ? (rs -> rs - rsu) :
                 (rs -> zero(rs - rsv))

# per-domain Ψ branch through the EXPORTED phase map (eq:toTilde): Ψ_dom = Ψ̃ e^{−imΩ h_dom}
function vsj_branch_psi(dom, m, Om, rsv, rsu)
    hb = vsj_branch_h(dom, rsv, rsu)
    rs -> psi_from_tilde(m, vsj_tld(rs), Om, hb(rs))
end

# (a, r0, r_v, r_u) cases: the CRP-config geometry at a = 0 plus spinning + retrograde draws
const VSJ_CASES = ((0.0, 10.0, 2.8, 24.4), (0.5, 10.0, 3.0, 26.0), (-0.7, 10.0, 4.0, 30.0))

@testset "MModeOperator vtu slicing + junction data (vtu-slicing-jumps)" begin

    @testset "exports" begin
        for s in (:VtuSlicing, :vtu_slicing, :vtu_height_function, :vtu_domain,
                  :VtuJunction, :vtu_junctions)
            @test Base.isexported(ScalarSelfForce, s)
        end
    end

    @testset "(a) construction + tortoise consistency" begin
        for (a, r0, rv, ru) in VSJ_CASES
            orb = CircularOrbit(a, r0)
            v = vtu_slicing(orb, rv, ru)
            @test v isa VtuSlicing{Float64}
            @test v.r_v == rv && v.r_u == ru
            # tortoise chart bitwise-consistent with the kinematics layer
            @test v.rstar_v == rstar(orb, rv)
            @test v.rstar_u == rstar(orb, ru)
            # paper constraint main.tex:323: r^v_* < r_{*0} < r^u_*
            @test v.rstar_v < rstar(orb, r0) < v.rstar_u
        end
        orb = CircularOrbit(0.0, 10.0)
        @test_throws ArgumentError vtu_slicing(orb, 24.4, 2.8)    # r_v ≥ r_u
        @test_throws ArgumentError vtu_slicing(orb, 3.0, 8.0)     # worldline not in t (r0 ≥ r_u)
        @test_throws ArgumentError vtu_slicing(orb, 12.0, 24.0)   # worldline not in t (r0 ≤ r_v)
        @test_throws ArgumentError vtu_slicing(orb, 1.5, 24.4)    # r_v ≤ r₊ = 2 (tortoise undefined)
    end

    @testset "(b) RunConfig round-trip (crp_paper_config)" begin
        cfg = crp_paper_config()
        orb = CircularOrbit(cfg.orbit.black_hole_spin, cfg.orbit.orbital_radius)
        v = vtu_slicing(orb, cfg)
        tr = cfg.orbit.hyperboloidal_slicing_transitions
        # the yaml 4-tuple is the degenerate (sharp) form (v→t, v→t, t→u, t→u); the radii
        # round-trip as the EXACT IEEE doubles the rendered yaml carries (10 − 4.8·1.5 =
        # 2.8000000000000007 in binary64, runtime-config-fixture gate (b))
        @test v.r_v === tr[1] === tr[2]
        @test v.r_u === tr[3] === tr[4]
        @test v.r_v == 10.0 - 4.8 * 1.5 && v.r_u == 10.0 + 4.8 * 3
        # transitions ARE radial block interfaces (bounds {inner, v→t, r0−wr, r0+wr, t→u, R}):
        # the interfaces are element faces, sampled from both sides — the DG jump premise
        @test v.r_v == cfg.domain.block_bounds_r[2]
        @test v.r_u == cfg.domain.block_bounds_r[5]
        # paper ordering r^v_* < r_{*0} < r^u_* for the production config
        @test v.rstar_v < rstar(orb, cfg.orbit.orbital_radius) < v.rstar_u
        # second Jinja2 parameter point: v→t = r0 − 1.5wr, t→u = r0 + 3wr
        cfg2 = crp_paper_config(orbital_radius = 8.0, worldtube_radius = 2.0)
        v2 = vtu_slicing(CircularOrbit(0.0, 8.0), cfg2)
        @test v2.r_v == 8.0 - 2.0 * 1.5 && v2.r_u == 8.0 + 2.0 * 3
        # orbit/config mismatch is refused (spin, then radius)
        @test_throws ArgumentError vtu_slicing(CircularOrbit(0.5, 10.0), cfg)
        @test_throws ArgumentError vtu_slicing(CircularOrbit(0.0, 9.0), cfg)
        # non-degenerate (smooth transition band) yaml tuples are NOT vtu — refused
        oc = ScalarSelfForce.OrbitConfig{Float64}(1.0, 0.0, 10.0, 2,
                                                  (2.8, 3.0, 24.4, 24.4), true, true)
        @test_throws ArgumentError vtu_slicing(CircularOrbit(0.0, 10.0), oc)
    end

    @testset "(c) piecewise boost H + height h exact (main.tex:316-330)" begin
        for (a, r0, rv, ru) in VSJ_CASES
            orb = CircularOrbit(a, r0)
            v = vtu_slicing(orb, rv, ru)
            hf = vtu_height_function(v)
            @test hf isa HeightFunction
            rsv, rsu = v.rstar_v, v.rstar_u
            vpts = (rsv - 40.0, rsv - 3.7, prevfloat(rsv))
            tpts = (nextfloat(rsv), (rsv + rsu) / 2, prevfloat(rsu))
            upts = (nextfloat(rsu), rsu + 5.1, rsu + 1.0e6)
            for rs in vpts
                @test hf.H(rs) == -1
                @test hf.h(rs) == rsv - rs          # bitwise: the printed v-branch law
            end
            for rs in tpts
                @test hf.H(rs) == 0
                @test hf.h(rs) == 0
            end
            for rs in upts
                @test hf.H(rs) == 1
                @test hf.h(rs) == rs - rsu          # bitwise: the printed u-branch law
            end
            # H′ ≡ 0 pointwise EVERYWHERE (main.tex:347: the H′ deltas at the transition
            # radii are handled by the junction conditions, not the pointwise function)
            for rs in (vpts..., rsv, tpts..., rsu, upts...)
                @test hf.Hprime(rs) == 0
            end
            # interface pins: h continuity anchors (exact zeros) + the H convention (the
            # pointwise function takes the t-domain value on the closed [r^v_*, r^u_*];
            # per-side values at the double-valued faces come from the domain branches)
            @test hf.h(rsv) == 0 && hf.h(rsu) == 0
            @test hf.H(rsv) == 0 && hf.H(rsu) == 0
            # ulp-scale continuity of h across both interfaces
            @test hf.h(prevfloat(rsv)) == rsv - prevfloat(rsv)
            @test hf.h(nextfloat(rsu)) == nextfloat(rsu) - rsu
            # hyperboloidal condition H(±∞) = ±1 (main.tex:278) — owned by the assembly
            @test hf.H(rsv - 1.0e12) == -1 && hf.H(rsu + 1.0e12) == 1
            # boost/height relation inside each open domain: FD of (linear) h equals H,
            # FD of (constant) H vanishes — exact to roundoff
            for (rs, d) in ((rsv - 8.0, 1.5), ((rsv + rsu) / 2, (rsu - rsv) / 8),
                            (rsu + 9.0, 2.25))
                @test vsj_cdiff(hf.h, rs, d) ≈ hf.H(rs) atol = 1e-14
                @test vsj_cdiff(hf.H, rs, d) == 0
            end
        end
        # type-generic evaluation: BigFloat in ⇒ BigFloat out (Float64 slicing)
        v = vtu_slicing(CircularOrbit(0.0, 10.0), 2.8, 24.4)
        hf = vtu_height_function(v)
        for rs in (big(v.rstar_v) - 5, big(v.rstar_v) + 3, big(v.rstar_u) + 2)
            @test hf.h(rs) isa BigFloat && hf.H(rs) isa BigFloat && hf.Hprime(rs) isa BigFloat
        end
        @test hf.h(big(v.rstar_v) - 5) == 5
    end

    @testset "(d) domain classification vtu_domain" begin
        orb = CircularOrbit(0.0, 10.0)
        v = vtu_slicing(orb, 2.8, 24.4)
        hf = vtu_height_function(v)
        @test vtu_domain(v, 2.1) === :v
        @test vtu_domain(v, prevfloat(v.r_v)) === :v
        @test vtu_domain(v, orb.r0) === :t
        @test vtu_domain(v, nextfloat(v.r_u)) === :u
        @test vtu_domain(v, 1.0e4) === :u
        # interface convention pinned (matches the closed-t pointwise H; the double-valued
        # faces are the junction-data seam, and Gauss volume nodes never sample a face)
        @test vtu_domain(v, v.r_v) === :t && vtu_domain(v, v.r_u) === :t
        # classification ↔ boost consistency at interior radii
        for (r, dom, Hval) in ((2.2, :v, -1), (10.0, :t, 0), (500.0, :u, 1))
            @test vtu_domain(v, r) === dom
            @test hf.H(rstar(orb, r)) == Hval
        end
        # crp_paper_config block columns classify v|t|t|t|u — the yaml NullSlicingBlocks
        # (Block0/Block5 v-sliced, Block4/Block9 u-sliced, middle columns t-sliced)
        cfg = crp_paper_config()
        vc = vtu_slicing(CircularOrbit(cfg.orbit.black_hole_spin,
                                       cfg.orbit.orbital_radius), cfg)
        bb = cfg.domain.block_bounds_r
        mids = [(bb[i] + bb[i + 1]) / 2 for i in 1:5]
        @test [vtu_domain(vc, r) for r in mids] == [:v, :t, :t, :t, :u]
    end

    @testset "(e) junction conditions [∂Ψ] = −imΩΨ (eq:jump_vt/eq:jump_tu, 1e-12)" begin
        setprecision(BigFloat, 256) do
            for (a, r0, rv, ru) in VSJ_CASES, m in (1, 2, 5)
                orb = CircularOrbit(a, r0)
                v = vtu_slicing(orb, rv, ru)
                jc = vtu_junctions(v, orb, m)
                # the DATA: both interfaces carry coeff = −imΩ at their r/rstar locations
                @test jc.vt isa VtuJunction && jc.tu isa VtuJunction
                @test jc.vt.r == v.r_v && jc.vt.rstar == v.rstar_v
                @test jc.tu.r == v.r_u && jc.tu.rstar == v.rstar_u
                @test jc.vt.coeff == -im * m * orb.Omega
                @test jc.tu.coeff == jc.vt.coeff
                # continuity of the GLOBAL-height Ψ across each interface (1e-14)
                hf = vtu_height_function(v)
                psig(rs) = psi_from_tilde(m, vsj_tld(rs), orb.Omega, hf.h(rs))
                for rsx in (v.rstar_v, v.rstar_u)
                    @test abs(psig(prevfloat(rsx)) - psig(nextfloat(rsx))) <= 1e-14
                end
                # one-sided limits: central FD (h = 1e-25, truncation ~h² = 1e-50) of the
                # per-domain branch that coincides with Ψ on that side of the interface
                fdh = big"1e-25"
                for (iface, lo, hi) in ((jc.vt, :v, :t), (jc.tu, :t, :u))
                    rsB = big(iface.rstar)
                    psilo = vsj_branch_psi(lo, m, orb.Omega, v.rstar_v, v.rstar_u)
                    psihi = vsj_branch_psi(hi, m, orb.Omega, v.rstar_v, v.rstar_u)
                    # branch values agree AT the interface (h branches share the zero)
                    @test abs(psihi(rsB) - psilo(rsB)) <= 1e-14
                    jump = vsj_cdiff(psihi, rsB, fdh) - vsj_cdiff(psilo, rsB, fdh)
                    target = big(iface.coeff) * psihi(rsB)   # −imΩ Ψ|_iface
                    @test abs(target) > 1e-3                  # non-vacuous gate
                    @test abs(jump - target) <= 1e-12 * abs(target)
                    # eq:jump_vt/eq:jump_tu orientation (above − below) is pinned: the
                    # flipped-sign target misses by ~2|target|
                    @test abs(jump + target) > 1e6 * abs(jump - target)
                end
            end
            # flux binding, one draw: F^r = A^rr ∂_rΨ IS ∂_{r_*}Ψ (main.tex:354-358 chart +
            # tsliced A^rr = dr/dr_*), so the SAME data reads [F^r] = coeff·Ψ — the form the
            # DG FaceState (φ = n·F) jump hooks consume
            orb = CircularOrbit(0.5, 10.0)
            v = vtu_slicing(orb, 3.0, 26.0)
            m = 2
            jc = vtu_junctions(v, orb, m)
            rB = big(v.r_v)
            psilo = vsj_branch_psi(:v, m, orb.Omega, v.rstar_v, v.rstar_u)
            psihi = vsj_branch_psi(:t, m, orb.Omega, v.rstar_v, v.rstar_u)
            Flo = flux_coeff_r(orb, rB) * vsj_cdiff(x -> psilo(rstar(orb, x)), rB, big"1e-20")
            Fhi = flux_coeff_r(orb, rB) * vsj_cdiff(x -> psihi(rstar(orb, x)), rB, big"1e-20")
            target = big(jc.vt.coeff) * psihi(big(v.rstar_v))
            @test abs((Fhi - Flo) - target) <= 1e-12 * abs(target)
        end
    end

    @testset "(f) m = 0 degeneracy — no jumps" begin
        orb = CircularOrbit(0.0, 10.0)
        v = vtu_slicing(orb, 2.8, 24.4)
        jc0 = vtu_junctions(v, orb, 0)
        @test iszero(jc0.vt.coeff) && iszero(jc0.tu.coeff)
        # the phase map degenerates to the identity: Ψ ≡ Ψ̃ in every domain
        hf = vtu_height_function(v)
        for rs in (v.rstar_v - 2.0, v.rstar_v, (v.rstar_v + v.rstar_u) / 2, v.rstar_u + 3.3)
            @test psi_from_tilde(0, vsj_tld(rs), orb.Omega, hf.h(rs)) == vsj_tld(rs)
        end
        # measured FD jump vanishes (branches are the same function at m = 0)
        setprecision(BigFloat, 256) do
            rsB = big(v.rstar_v)
            psilo = vsj_branch_psi(:v, 0, orb.Omega, v.rstar_v, v.rstar_u)
            psihi = vsj_branch_psi(:t, 0, orb.Omega, v.rstar_v, v.rstar_u)
            @test abs(vsj_cdiff(psihi, rsB, big"1e-25") -
                      vsj_cdiff(psilo, rsB, big"1e-25")) <= 1e-40
        end
    end

    @testset "(g) BigFloat inheritance" begin
        setprecision(BigFloat, 256) do
            orbB = CircularOrbit(big(0.5), big(10.0))
            vB = vtu_slicing(orbB, big(3.0), big(26.0))
            @test vB isa VtuSlicing{BigFloat}
            @test precision(vB.rstar_v) == 256 && precision(vB.rstar_u) == 256
            @test vB.rstar_v == rstar(orbB, big(3.0))
            hfB = vtu_height_function(vB)
            @test hfB.h(vB.rstar_v - 1) == 1 && hfB.h(vB.rstar_v - 1) isa BigFloat
            @test hfB.H(vB.rstar_u + 2) == 1 && hfB.H(vB.rstar_u + 2) isa BigFloat
            jcB = vtu_junctions(vB, orbB, 3)
            @test jcB.vt.coeff isa Complex{BigFloat}
            @test jcB.vt.coeff == -im * 3 * orbB.Omega
            @test vtu_domain(vB, big(4.0)) === :t && vtu_domain(vB, big(2.9)) === :v
        end
    end
end
