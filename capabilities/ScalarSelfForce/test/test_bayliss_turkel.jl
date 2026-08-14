# test_bayliss_turkel.jl — ESSF node `bayliss-turkel-bc` (schedule W7-D, GLOBAL_DAG row 32),
# written BEFORE the implementation (TDD). The node delivers the second-order Bayliss–Turkel
# outer radiation condition of arXiv:2606.04998 (following Osburn:2022bby):
#
#   t-sliced (main.tex:403-407):  ∂²_{r_*}Ψ̃ − 2imΩ ∂_{r_*}Ψ̃ − m²Ω²Ψ̃ = 0
#                                 = (∂_{r_*} − imΩ)² Ψ̃,
#   u-sliced (main.tex:409-411, eq:bayliss_turkel_bc):  ∂²_{r_*}Ψ = 0,
#   DG realization (main.tex:614-618): (n·F)_N = −Δ/(r²+a²)·(β_m/γ^r_m)·Ψ with the ghost
#   data (n·F)_ext = 2(n·F)_N − (n·F)_int, Ψ_ext = Ψ_int (Vu:2024),
#
# as `MModeOperator/bayliss_turkel.jl`: `bt_tsliced_residual` (the printed t-form),
# `bt_robin_coefficient` (the pointwise (n·F)_N = robin·Ψ coefficient, all four branches of
# the spectre Sommerfeld oracle, Sommerfeld.cpp:34-69), `BaylissTurkelGhost` (the
# EllipticCore GhostBoundary policy — linear-homogeneous, so linearized == full), and
# `outer_boundary_ghost` (the RunConfig m = 0 no-outer-BC flag semantics, yaml:115-125).
#
# PREDECESSOR SCOPE (sealed; cited, never re-derived here):
#   - elliptic-operator-tsliced / first-order-flux-form: β_m, γ^r_m, A^rr = Δ/(r²+a²) and
#     the assembled per-region callables (row-30 spectral operator identity at 5.28e-14) —
#     gate (d) consumes fluxform_callables(:u) as THE coefficient source;
#   - hyperboloidal-transformation: the replacement rules ∂_{r_*} → ∂_{r_*} + imΩH,
#     ∂²_{r_*} → ∂²_{r_*} + 2imΩH∂_{r_*} + imΩH′ − m²Ω²H² (main.tex:306-308) and
#     rstar(orbit, r) — gate (c) applies them to the t-form COEFFICIENT TRIPLE;
#   - boundary-data-hooks (EllipticCore W5-B): the ghost seam ghost_state /
#     ghost_state_linearized + conventions N1-N7 — this node only ADDS a GhostBoundary
#     subtype, never re-plumbs;
#   - runtime-config-fixture (W4-E): SommerfeldBC + the m == 0 ⇒ outer_radial_bc === nothing
#     branch (gated there against the yaml; re-pinned here only as the input of gate (f)).
# LATER NODE (not here): dg-boundary-conditions (row 41) wires the ghost into the Mesh2D
# sweep (Solve/bc_glue.jl) and runs the reduced-ODE domain-contamination oracle.
#
# THIS node's gates (constants calibrated pre-red in scratchpad calibrate_bayliss_turkel.jl):
#  (a) exports.
#  (b) the t-form ANNIHILATES the outgoing asymptotics Ψ̃ = e^{imΩr_*}/r^n, n ∈ {1, 2}
#      (symbolic-in-test differentiation: closed-form ∂_{r_*}, ∂²_{r_*} via α = Δ/(r²+a²) =
#      dr/dr_*, hand-derived below; BigFloat(256) — the cancellation depth ~1e13 at r = 1e8
#      is invisible to Float64): relative residual |B|/(m²Ω²|Ψ̃|) ≤ 1e-12 at r = 1e8
#      (measured worst 6.19e-13, m = 1, n = 2, a = 0.5) falling as 1/r² (ratio 1e-4 per two
#      decades, measured 1.0000e-4); the EXACT kernel (c₁ + c₂r_*)e^{imΩr_*} is annihilated
#      to BigFloat roundoff ≤ 1e-60 (measured ≤ 2.2e-77; B = (∂_{r_*} − imΩ)² exactly);
#      FALSIFIER: the ingoing mode e^{−imΩr_*} gives relative residual 4 (measured 4.0).
#  (c) u-transform equivalence (main.tex:408-411): pushing the t-form coefficient triple
#      (probed from bt_tsliced_residual by unit vectors — the shipped equation is what is
#      checked) through the hyperboloidal replacement rules with the u-slab boost H = +1,
#      H′ = 0 yields EXACTLY (∂²_{r_*}, 0, 0) — eq:bayliss_turkel_bc; H = 0 is the identity
#      (t-domain, bitwise), H = −1 gives (∂_{r_*} − 2imΩ)² (v-domain ingoing annihilator,
#      the horizon-side counterpart). Float64 ≤ 4 eps scale + BigFloat(256) ≤ 1e-70.
#  (d) spectre HARD REGRESSION values (audit-spectre-ssf §4.4; Test_Sommerfeld.cpp:91-93 +
#      123-124, M = 1, χ = 0.5, r0 = 10, m = 2, Ψ = 1): no-slicing order 1 ⇒
#      n·F = i·0.06226111848298833 (= i·mΩ, Sommerfeld.cpp:59-60) at 1e-12; hyperboloidal
#      order 2 with (α, β, γ_r) = (1, 0.5, 2) ⇒ n·F = −0.25 (= −αβ/γ_r, Sommerfeld.cpp:45-46)
#      at 1e-12; hyperboloidal order 1 ⇒ 0 exactly (plain outflow, Sommerfeld.cpp:43);
#      no-slicing order 2 ⇒ (k² − β)/(γ_r − 2ik) (Sommerfeld.cpp:61-64, transcribed as-built);
#      order 3 refused in both slicings (spectre ERROR).
#  (e) DG-flux derivation + ghost-policy exactness:
#      (e1) main.tex:616 through the OPERATOR's own u-region coefficients: robin ==
#           −A^rr·β_m/γ^r_m with fluxform_callables(:u) at R = 1e4 (identical expression,
#           ==), and the derivation's "solved for ∂Ψ" step: the imposed slope
#           ∂_rΨ = robin·Ψ/A^rr makes β_mΨ + γ^r_m∂_rΨ vanish (≤ 1e-13 rel) — which IS
#           ∂_rF^r = 0 ⇔ ∂²_{r_*}Ψ = 0 substituted into the interior PDE.
#      (e2) the ghost on a MANUFACTURED u-sliced outgoing solution (a = 0 ⇒ exactly
#           z-independent; series Ψ = 1 + b₁/r + b₂/r² with b₁, b₂ derived below from the
#           radial ODE's asymptotic balance): Ψ_ext == Ψ_int (exact), (nF)* == robin·Ψ
#           (Neumann-on-average, ≤ 1e-13 rel, penalty inert — PenaltyFlux == CentralFlux
#           bitwise), and the boundary residual |(nF)* − (n·F)_int| sits AT the physical
#           reflection scale of the second-order BT condition, (m+1)/(2mΩ²R³) — within
#           [0.9, 1.1]× (measured 0.99995×), falling as 1/R³ (ratio 1.0000e-3 per decade);
#           FALSIFIER: the ingoing wave e^{−2imΩr_*} leaves residual ~2mΩ (measured 0.063,
#           m = 1) — ≥ 1e6× the outgoing one. Linearized ghost bitwise == full (the BC is
#           linear-homogeneous: it lives in the GMRES matrix, spectre Sommerfeld.cpp:71-81).
#  (f) m = 0 flag honored + config plumbing + BigFloat(256) precision inheritance.
using Test
using ScalarSelfForce
using ScalarSelfForce.EllipticCore: FaceState, GhostBoundary, ghost_state,
                                    ghost_state_linearized, PenaltyFlux, CentralFlux,
                                    numerical_flux, element_matrices, Block, ElementId,
                                    SegmentId, LinearMap1D, InverseRadialMap1D

@testset "MModeOperator Bayliss-Turkel outer BC (bayliss-turkel-bc)" begin

    @testset "(a) exports" begin
        for nm in (:bt_tsliced_residual, :bt_robin_coefficient, :BaylissTurkelGhost,
                   :outer_boundary_ghost)
            @test Base.isexported(ScalarSelfForce, nm)
        end
    end

    @testset "(b) t-form annihilates outgoing asymptotics e^{imΩr_*}/r^n (1e-12 at r=1e8)" begin
        setprecision(BigFloat, 256) do
            # symbolic-in-test differentiation of Ψ̃ = E·f, E = e^{imΩr_*}, f = r^{−n}:
            # ∂_{r_*}(E·g) = E·(imΩ g + α g′) with α = Δ/(r²+a²) = dr/dr_* (main.tex:356-358)
            #   ⇒ ∂_{r_*}Ψ̃  = E·(imΩ f + α f′)
            #   ⇒ ∂²_{r_*}Ψ̃ = E·(−m²Ω² f + 2imΩ α f′ + α α′ f′ + α² f″)
            # (α′ = dα/dr; the r_*-phase makes every m²Ω²-sized term cancel INSIDE the
            # t-form, leaving B = E·∂²_{r_*}f = O(1/r³) — the BT annihilation property)
            for (a, r0) in ((big"0.0", big"10.0"), (big"0.5", big"10.0"),
                            (big"-0.7", big"10.0"))
                orb = CircularOrbit(a, r0)
                Om = orb.Omega
                for m in (1, 2, 3), n in (1, 2)
                    rel = Dict{Float64,BigFloat}()
                    for r in (big"1.0e6", big"1.0e8")
                        w  = r^2 + orb.a^2
                        D  = (r - orb.rp) * (r - orb.rm)
                        al = D / w
                        alp = ((2r - 2orb.M) * w - 2r * D) / w^2
                        f, fp, fpp = r^(-n), -n * r^(-n - 1), n * (n + 1) * r^(-n - 2)
                        E = exp(im * m * Om * rstar(orb, r))
                        psi   = E * f
                        dpsi  = E * (im * m * Om * f + al * fp)
                        d2psi = E * (-(m * Om)^2 * f + 2im * m * Om * al * fp +
                                     al * alp * fp + al^2 * fpp)
                        B = bt_tsliced_residual(orb, m, psi, dpsi, d2psi)
                        rel[Float64(r)] = abs(B) / ((m * Om)^2 * abs(psi))
                    end
                    @test rel[1.0e8] <= 1e-12          # calibrated worst 6.19e-13
                    @test rel[1.0e6] <= 1e-8
                    # BT reflection-error decay ~1/r²: exactly 1e-4 per two decades
                    @test 0.9e-4 <= rel[1.0e8] / rel[1.0e6] <= 1.1e-4
                end
            end
            # EXACT kernel: B = (∂_{r_*} − imΩ)² annihilates (c₁ + c₂ r_*)e^{imΩr_*}
            # identically — closed-form derivatives, residual is pure BigFloat roundoff
            orb = CircularOrbit(big"0.5", big"10.0")
            Om = orb.Omega
            c1, c2 = big"0.7" - big"0.2" * im, big"0.3" + big"0.5" * im
            for m in (1, 3), r in (big"30.0", big"1.0e4")
                rs = rstar(orb, r)
                q  = c1 + c2 * rs
                E  = exp(im * m * Om * rs)
                B = bt_tsliced_residual(orb, m, E * q, E * (c2 + im * m * Om * q),
                                        E * (2im * m * Om * c2 - (m * Om)^2 * q))
                @test abs(B) / ((m * Om)^2 * abs(q)) <= 1e-60    # measured ≤ 2.2e-77
                # FALSIFIER: the ingoing mode e^{−imΩr_*} is NOT annihilated —
                # B = −4m²Ω²ψ (relative residual exactly 4)
                psi_i = exp(-im * m * Om * rs)
                B_i = bt_tsliced_residual(orb, m, psi_i, -im * m * Om * psi_i,
                                          -(m * Om)^2 * psi_i)
                @test abs(B_i) / ((m * Om)^2 * abs(psi_i)) >= 1
                @test abs(abs(B_i) / ((m * Om)^2 * abs(psi_i)) - 4) <= 1e-70
            end
        end
    end

    @testset "(c) u-transform: t-form ⇒ ∂²_{r_*}Ψ = 0 identically (eq:bayliss_turkel_bc)" begin
        # probe the SHIPPED t-form's coefficient triple (c₂, c₁, c₀) on (∂², ∂, 1) by unit
        # vectors (bt_tsliced_residual is linear), then push it through the hyperboloidal
        # replacement rules (main.tex:306-308, sealed hyperboloidal-transformation node):
        #   ∂_{r_*} → ∂_{r_*} + imΩH,   ∂²_{r_*} → ∂²_{r_*} + 2imΩH ∂_{r_*} + imΩH′ − m²Ω²H²
        # i.e. on the triple: n₂ = c₂, n₁ = c₁ + 2imΩH c₂,
        #                     n₀ = c₀ + (imΩH′ − m²Ω²H²)c₂ + imΩH c₁.
        subst(k, c2, c1, c0, H, Hp) = (c2, c1 + 2im * k * H * c2,
                                       c0 + (im * k * Hp - k^2 * H^2) * c2 + im * k * H * c1)
        z, o = 0.0 + 0im, 1.0 + 0im
        for (a, r0) in ((0.5, 10.0), (0.0, 10.0)), m in (1, 2, 5)
            orb = CircularOrbit(a, r0)
            k = m * orb.Omega
            c2 = bt_tsliced_residual(orb, m, z, z, o)
            c1 = bt_tsliced_residual(orb, m, z, o, z)
            c0 = bt_tsliced_residual(orb, m, o, z, z)
            @test c2 == 1 && c1 == -2im * k && c0 == -k^2   # the printed main.tex:405 triple
            # u domain (H = +1, H′ = 0; vtu boost table main.tex:316-322):
            # the t-form collapses to ∂²_{r_*}Ψ = 0 — eq:bayliss_turkel_bc, main.tex:409-411
            n2, n1, n0 = subst(k, c2, c1, c0, 1, 0)
            @test n2 == 1
            @test abs(n1) <= 4 * eps() * k
            @test abs(n0) <= 4 * eps() * k^2
            # t domain (H = 0): identity, bitwise
            @test subst(k, c2, c1, c0, 0, 0) == (c2, c1, c0)
            # v domain (H = −1): (∂_{r_*} − 2imΩ)² — the ingoing/horizon-side annihilator
            n2v, n1v, n0v = subst(k, c2, c1, c0, -1, 0)
            @test n2v == 1
            @test abs(n1v - (-4im * k)) <= 4 * eps() * 4k
            @test abs(n0v - (-4 * k^2)) <= 4 * eps() * 4k^2
        end
        # exact-arithmetic twin: BigFloat(256), residual pure roundoff
        setprecision(BigFloat, 256) do
            orb = CircularOrbit(big"0.5", big"10.0")
            m = 3
            k = m * orb.Omega
            zb, ob = zero(Complex{BigFloat}), one(Complex{BigFloat})
            c2 = bt_tsliced_residual(orb, m, zb, zb, ob)
            c1 = bt_tsliced_residual(orb, m, zb, ob, zb)
            c0 = bt_tsliced_residual(orb, m, ob, zb, zb)
            n2, n1, n0 = subst(k, c2, c1, c0, 1, 0)
            @test n2 == 1 && abs(n1) <= 1e-70 && abs(n0) <= 1e-70
        end
    end

    @testset "(d) spectre hard regression values (Test_Sommerfeld.cpp:91-93, 123-124)" begin
        # the spectre test point: M = 1, χ = 0.5, r0 = 10, m = 2, Ψ = 1 on a 5-point face —
        # pointwise here, the coefficient IS the n·F value at Ψ = 1
        # (1) no slicing, order 1: n·F = i·mΩ with Ω = 1/(a + √(r0³/M)) (Sommerfeld.cpp:54-60)
        bc1 = SommerfeldBC(1.0, 0.5, 10.0, 2, false, 1)
        c = bt_robin_coefficient(bc1, 1.0, 0.0 + 0im, 0.0 + 0im)
        @test abs(c - im * 0.06226111848298833) <= 1e-12     # the hard literal, cpp:91-93
        @test real(c) == 0
        # cross-check vs the orbit adapter's Ω = u^φ/u^t (independent route, sub-ulp)
        @test abs(c - 2im * CircularOrbit(0.5, 10.0).Omega) <= 4 * eps()
        # (2) hyperboloidal, order 2 with (α, β, γ_r) = (1, 0.5, 2): n·F = −αβ/γ_r = −0.25
        bc2 = SommerfeldBC(1.0, 0.5, 10.0, 2, true, 2)
        c2v = bt_robin_coefficient(bc2, 1.0, 0.5 + 0im, 2.0 + 0im)
        @test abs(c2v - (-0.25)) <= 1e-12                    # the hard literal, cpp:123-124
        @test c2v == -0.25
        # (3) hyperboloidal, order 1: plain outflow n·F = 0 (Sommerfeld.cpp:43)
        bch1 = SommerfeldBC(1.0, 0.5, 10.0, 2, true, 1)
        @test bt_robin_coefficient(bch1, 1.0, 0.5 + 0im, 2.0 + 0im) == 0
        # (4) no slicing, order 2: (k² − β)/(γ_r − 2ik) (Sommerfeld.cpp:61-64, as-built —
        # spectre leaves β un-multiplied by α in this unused-in-production branch)
        bc22 = SommerfeldBC(1.0, 0.5, 10.0, 2, false, 2)
        k = 2 / (0.5 + sqrt(1000.0))
        @test bt_robin_coefficient(bc22, 1.0, 0.5 + 0im, 2.0 + 0im) ==
              (k^2 - (0.5 + 0im)) / ((2.0 + 0im) - 2im * k)
        # (5) unimplemented orders refused (spectre ERROR, Sommerfeld.cpp:47-51, 65-68)
        @test_throws ArgumentError bt_robin_coefficient(
            SommerfeldBC(1.0, 0.5, 10.0, 2, true, 3), 1.0, 0.5 + 0im, 2.0 + 0im)
        @test_throws ArgumentError bt_robin_coefficient(
            SommerfeldBC(1.0, 0.5, 10.0, 2, false, 3), 1.0, 0.5 + 0im, 2.0 + 0im)
    end

    @testset "(e1) main.tex:616 through the operator's own u-region coefficients" begin
        bc = SommerfeldBC(1.0, 0.0, 10.0, 1, true, 2)   # spin/m overwritten per case below
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0)), m in (1, 2, 3)
            orb = CircularOrbit(a, r0)
            v = vtu_slicing(orb, 2.8, 24.4)
            fc = fluxform_callables(orb, m, v, :u)
            bcm = SommerfeldBC(1.0, a, r0, m, true, 2)
            R = 1.0e4
            for z in (0.0, 0.3, 0.8)
                al, be, ga = fc.A[1](R, z), fc.beta(R, z), fc.gamma[1](R, z)
                robin = bt_robin_coefficient(bcm, al, be, ga)
                # the printed choice (main.tex:616), identical expression: −Δ/(r²+a²)·β/γ^r
                @test robin == -al * be / ga
                # derivation ("solved for ∂Ψ"): the imposed slope ∂_rΨ = robin·Ψ/A^rr
                # zeroes β_mΨ + γ^r_m ∂_rΨ — i.e. the interior PDE with ∂_rF^r = 0, which
                # is EXACTLY ∂²_{r_*}Ψ = 0 (gate (c)) substituted into eq:first-order
                for psi in (1.0 + 0im, 0.7 - 0.4im)
                    @test abs(be * psi + ga * (robin * psi / al)) <= 1e-13 * abs(be * psi)
                end
            end
        end
    end

    @testset "(e2) ghost-policy exactness on a manufactured outgoing solution" begin
        # a = 0: the operator is exactly z-independent (β's a²(1−z²) term vanishes), so the
        # radial reduced ODE is the whole story on the outer face. u-sliced (H = +1) a = 0
        # closed forms (sealed tsliced/W6-E coefficients; re-pinned in calibration):
        #   A^rr = 1 − 2M/r,  β_m = (m(m+1) + 2M/r)/r²,  γ^r_m = −2imΩ.
        # Outgoing series Ψ = 1 + b₁/r + b₂/r² (u-sliced image of e^{imΩr_*}(1 + …)/r·r —
        # the mode function has 1/r factored, so the leading outgoing behavior is O(1)):
        # matching the ODE −∂_r(A^rr∂_rΨ) + βΨ + γ^r∂_rΨ = 0 at orders 1/r², 1/r³ gives
        #   b₁ = i m(m+1)/(2mΩ),   b₂ = −[(m(m+1) − 2) b₁ + 2M]/(4imΩ).
        # The ghost's boundary residual on this series is the PHYSICAL second-order-BT
        # reflection scale (m+1)/(2mΩ²R³) — calibrated 0.99995× at R = 1e4.
        orb = CircularOrbit(0.0, 10.0)
        Om, M = orb.Omega, orb.M
        v = vtu_slicing(orb, 2.8, 24.4)
        # face machinery: the production u-block outer face (Inverse 1/r map, eq:compaction)
        blk = Block(4, InverseRadialMap1D(24.4, 1.0e4), LinearMap1D(0.0, 1.0))
        em = element_matrices(Float64, blk, ElementId(4, SegmentId(0, 0), SegmentId(0, 0)),
                              (4, 4))
        fm = em.faces[1][2]                       # radial UPPER face (d = 1, s = +1)
        zf = [0.1, 0.3, 0.5, 0.7]                 # face z-points (z-independent at a = 0)
        for m in (1, 2)
            cfg = crp_paper_config(spin = 0.0, m_mode = m)
            fc = fluxform_callables(orb, m, v, :u)
            alv = [fc.A[1](1.0e4, z) for z in zf]
            bev = [fc.beta(1.0e4, z) for z in zf]
            gav = [fc.gamma[1](1.0e4, z) for z in zf]
            ghost = outer_boundary_ghost(cfg, alv, bev, gav)
            @test ghost isa BaylissTurkelGhost
            @test ghost isa GhostBoundary          # the EllipticCore seam subtype
            b1 = im * m * (m + 1) / (2 * m * Om)
            b2 = -((m * (m + 1) - 2) * b1 + 2M) / (4im * m * Om)
            deltas = Dict{Float64,Float64}()
            for R in (1.0e4, 1.0e5)
                psi  = 1 + b1 / R + b2 / R^2
                dpsi = -b1 / R^2 - 2 * b2 / R^3
                alR = fc.A[1](R, 0.3)
                phi = alR * dpsi                   # (n·F)_int, n = +ê_r (convention N1)
                robv = [bt_robin_coefficient(cfg.domain.outer_radial_bc, fc.A[1](R, z),
                                             fc.beta(R, z), fc.gamma[1](R, z)) for z in zf]
                gR = R == 1.0e4 ? ghost : BaylissTurkelGhost(robv)
                int = FaceState(fill(complex(psi), 4), fill(complex(phi), 4))
                ext = ghost_state(gR, int)
                # Ψ_ext = Ψ_int (main.tex:618) — exact, so the penalty term is INERT
                @test ext.psi == int.psi
                # (n·F)_ext = 2(n·F)_N − (n·F)_int in receiver terms = φ_int − 2·robin·ψ in
                # the N3 sender convention
                @test ext.nflux == int.nflux .- 2 .* (gR.robin .* int.psi)
                # linearized == full, bitwise: the BC is linear-homogeneous in the trial
                # data — it belongs to the GMRES matrix (spectre Sommerfeld.cpp:71-81)
                lin = ghost_state_linearized(gR, int)
                @test lin.psi == ext.psi && lin.nflux == ext.nflux
                # through the numerical flux: Ψ* == ψ_int exact; (nF)* == robin·ψ
                # (Neumann-on-average); penalty flux == central flux bitwise (zero jump)
                psis, nfs = numerical_flux(PenaltyFlux(), fm, complex(alv), int, ext)
                psic, nfc = numerical_flux(CentralFlux(), fm, complex(alv), int, ext)
                @test psis == int.psi && psic == int.psi
                @test nfs == nfc
                @test maximum(abs.(nfs .- gR.robin .* int.psi)) <=
                      1e-13 * maximum(abs.(gR.robin .* int.psi))
                # the boundary residual sits AT the physical BT reflection scale
                delta = maximum(abs.(nfs .- int.nflux))
                scale = (m + 1) / (2 * m * Om^2 * R^3)
                @test 0.9 * scale <= delta <= 1.1 * scale    # calibrated 0.99995×
                deltas[R] = delta
            end
            # 1/R³ decay of the reflection error (calibrated ratio 1.0000e-3 per decade)
            @test 0.9e-3 <= deltas[1.0e5] / deltas[1.0e4] <= 1.1e-3
            # FALSIFIER: the ingoing wave Ψ = e^{−2imΩ(r_*−r_*(R))} (u-sliced image of
            # e^{−imΩr_*}) is REJECTED: ψ(R) = 1, (n·F)_int = −2imΩ, residual ~2mΩ
            R = 1.0e4
            int_i = FaceState(fill(1.0 + 0im, 4), fill(-2im * m * Om, 4))
            ext_i = ghost_state(ghost, int_i)
            _, nfs_i = numerical_flux(PenaltyFlux(), fm, complex(alv), int_i, ext_i)
            delta_in = maximum(abs.(nfs_i .- int_i.nflux))
            @test delta_in >= 1e-2                            # measured 0.063 (m=1)
            @test delta_in / deltas[1.0e4] >= 1e6             # measured 6.3e7
        end
    end

    @testset "(f) m = 0 no-outer-BC flag + config plumbing + BigFloat(256)" begin
        # m = 0: the yaml Jinja2 branch (yaml:115-125) ⇒ no outer radial BC at all —
        # RunConfig carries nothing (sealed W4-E), and the ghost entry honors it
        cfg0 = crp_paper_config(m_mode = 0)
        @test cfg0.domain.outer_radial_bc === nothing
        @test outer_boundary_ghost(cfg0, [1.0], [0.5 + 0im], [2.0 + 0im]) === nothing
        # m ≠ 0: the config's order-2 hyperboloidal SommerfeldBC drives the ghost, robin
        # broadcast bitwise
        cfg = crp_paper_config(spin = 0.5, m_mode = 2)
        bc = cfg.domain.outer_radial_bc
        @test bc isa SommerfeldBC && bc.hyperboloidal_slicing && bc.order == 2
        alv, bev, gav = [1.0, 0.9], [0.5 + 0.1im, 0.4 - 0.2im], [2.0 + 1im, 1.5 - 0.5im]
        og = outer_boundary_ghost(cfg, alv, bev, gav)
        @test og isa BaylissTurkelGhost
        @test og.robin == bt_robin_coefficient.(Ref(bc), alv, bev, gav)
        # BigFloat(256): precision inherited end to end (package type policy)
        setprecision(BigFloat, 256) do
            cfgB = crp_paper_config(spin = big"0.5", m_mode = 2)
            orbB = CircularOrbit(big"0.5", big"10.0")
            vB = vtu_slicing(orbB, big"2.8", big"24.4")
            fcB = fluxform_callables(orbB, 2, vB, :u)
            rB = big"1.0e4"
            cB = bt_robin_coefficient(cfgB.domain.outer_radial_bc, fcB.A[1](rB, big"0.3"),
                                      fcB.beta(rB, big"0.3"), fcB.gamma[1](rB, big"0.3"))
            @test cB isa Complex{BigFloat} && precision(real(cB)) == 256
            gB = BaylissTurkelGhost([cB])
            sB = ghost_state(gB, FaceState([big"1.0" + 0im], [big"0.5" + 0im]))
            @test eltype(sB.nflux) == Complex{BigFloat}
            @test precision(real(sB.nflux[1])) == 256
            # the t-sliced branches inherit precision through k = mΩ too
            kB = bt_robin_coefficient(SommerfeldBC(big"1.0", big"0.5", big"10.0", 2,
                                                   false, 1),
                                      big"1.0", big"0.0" + 0im, big"0.0" + 0im)
            @test kB isa Complex{BigFloat} && precision(imag(kB)) == 256
        end
    end
end
