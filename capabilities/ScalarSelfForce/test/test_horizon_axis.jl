# test_horizon_axis.jl — ESSF node `fix-horizon-axis-layer` (GLOBAL_DAG row 58), written
# BEFORE the implementation (TDD). The node delivers `Solve/horizon_axis.jl` — the instrument
# that MEASURES regularity of the solved `Ψ_m` on the two DEGENERATE boundaries of the CRP
# domain (the horizon `r = r₊` and the polar axis `z² = 1`), the control parameter of the
# spurious layer that lives on their common corner, and the one refinement rule that removes
# it.
#
# ORACLE (the read-only mirrors, transcribed by hand):
#   * `main.tex:365-383` — `𝔄^rr = Δ/(r²+a²)`,
#     `γ^r_m = 2[a²Δ/(r(r²+a²)²) − ima/(r²+a²) − imΩH]`.
#   * `main.tex:385-386` — the normal flux vanishes at `r = r₊` and at `z² = 1`; the `1/Δ`
#     terms of `β_m` cancel at `H = −1`.
#   * `main.tex:396-401` — "The vanishing of `F^r_m` at the horizon means that we do not have
#     to explicitly impose boundary conditions there. … These regularity conditions are
#     particularly suitable for spectral numerical methods, because they expand fields in a
#     polynomial series that inherently selects regular solutions."
#   * `main.tex:614-615` — "on boundaries where the perpendicular flux vanishes,
#     `n_i𝔄^i_m = 0`, neither of the two boundary terms contributes."
#   * `main.tex:227-231` (eq:mode_sum) — the SOLVED variable already carries `sin^{|m|}θ`
#     divided out, so `Ψ_m` approaches a FINITE, generically non-zero limit at the axis. (
#     `arXiv:2206.07031` §BC imposes `Ψ_{m≠0}|_{θ=0,π} = 0` on ITS mode variable precisely
#     because that one does not carry the factoring; main.tex:399-401 calls the two
#     statements equivalent.)
#   * `arXiv:2206.07031` eq. horBC — the ingoing near-horizon condition
#     `(∂_{r_*}Ψ_m − i k^-Ψ_m)|_{r_*→−∞} = 0`, `k^- = −mΩ`, which under the CRP's `vtu`
#     v-domain boost `H = −1` is `∂_{r_*}Ψ_m = 0`, i.e. `n·F = 𝔄^rr∂_rΨ_m = 0` — a condition
#     EVERY polynomial satisfies identically at `r₊`, which is why no face condition can
#     select the branch and the trial space must.
#
# GATES. (a)-(b) pin the instrument's arithmetic against closed forms and against the sealed
# coefficients. (c) is the healthy signature where the layer is absent. (d) is the node's own
# gate: with the horizon order hold the layer is gone at the paper's refinement depth, and
# WITHOUT it (the pre-row-58 sweep, reproduced bitwise by `horizon_hold = false`) it is
# present — the falsifier, so a future change cannot make the gate pass vacuously.

using Test
using Printf
using ScalarSelfForce
using EllipticCore
using Spacetime: kerr_isco

@testset "horizon-axis regularity (row 58)" begin

    # ---- (a) the branch index against its closed form ------------------------------------
    @testset "horizon_branch_index closed form" begin
        for a in (0.0, 0.5, -0.5, 0.9, 0.998, -0.998), m in (0, 1, 2, 17)
            orb = CircularOrbit(a, kerr_isco(a))
            w = orb.rp^2 + orb.a^2
            OmH = orb.a / w                       # horizon angular velocity a/(2Mr₊)
            expect = 2 * m * abs(orb.Omega - OmH) * w / (orb.rp - orb.rm)
            @test horizon_branch_index(orb, m) ≈ expect rtol = 1e-13
        end
        @test horizon_branch_index(CircularOrbit(0.9, kerr_isco(0.9)), 0) == 0
        orb0 = CircularOrbit(0.0, 6.0)            # a = 0: Ω_H = 0, r₊ = 2, r₋ = 0 ⇒ ν = 4mΩ
        @test horizon_branch_index(orb0, 17) ≈ 4 * 17 * orb0.Omega rtol = 1e-13
        orbk = CircularOrbit(0.998, kerr_isco(0.998))
        @test horizon_branch_index(orbk, 34) ≈ 2 * horizon_branch_index(orbk, 17) rtol = 1e-13
        @test 27 < horizon_branch_index(orbk, 17) < 28      # the layer's own configuration
    end

    # ---- (b) it IS the indicial root of the sealed coefficients ---------------------------
    @testset "branch index == |γ^r(r₊)| / (dA^rr/dr)(r₊)" begin
        for a in (0.5, 0.998, -0.998)
            orb = CircularOrbit(a, kerr_isco(a))
            h = 1e-6
            c_fd = (ScalarSelfForce.flux_coeff_r(orb, orb.rp + h) -
                    ScalarSelfForce.flux_coeff_r(orb, orb.rp)) / h
            g = ScalarSelfForce.gamma_r_m(orb, 17, orb.rp, -1)
            @test isapprox(horizon_branch_index(orb, 17), abs(g) / c_fd; rtol = 1e-5)
            @test real(g) == 0                 # H = −1: the horizon value is purely imaginary
        end
    end

    # ---- (c) the healthy signature, where ν is small --------------------------------------
    @testset "healthy profile at a = 0, m = 17 (ν = 4.6)" begin
        a, m, r0 = 0.0, 17, 6.0
        orb = CircularOrbit(a, r0)
        cfg = crp_paper_config(; spin = a, orbital_radius = r0, m_mode = m,
                               inner_radius = orb.rp, outer_radius = 1.0e4,
                               worldtube_radius = 0.6 * (r0 - orb.rp),
                               worldtube_dz2 = 0.25, P = 4, L = 0)
        hist = solve_mode_amr(orb, cfg, m; iterations = 2, horizon_hold = true)
        p = horizon_axis_profile(hist.solution)
        @printf("  [row 58] a = 0, m = 17, nu = %.2f: r = %.8f, ratio = %.4e, q = %+.3f, layer/plane = %.4f\n",
                horizon_branch_index(orb, m), p.r, p.ratio, p.exponent, p.layer_ratio)
        @test p.r ≈ orb.rp atol = 0.05
        @test length(p.absu) == length(p.zeta) >= 8
        @test issorted(p.zeta)
        @test 1e-2 < p.ratio < 1e2             # O(1): the mode variable is regular at the axis
        @test abs(p.exponent) < 1              # no (1−z²)^{−|m|} branch
        @test p.layer_ratio <= 1.05            # the horizon strip never dominates the plane
    end

    # ---- (d) the node's gate, with its falsifier ------------------------------------------
    @testset "a = +0.998 at the ISCO, m = 17, 7 AMR grids" begin
        a, m = 0.998, 17
        r0 = kerr_isco(a)
        orb = CircularOrbit(a, r0)
        cfg = crp_paper_config(; spin = a, orbital_radius = r0, m_mode = m,
                               inner_radius = orb.rp, outer_radius = 1.0e4,
                               worldtube_radius = 0.6 * (r0 - orb.rp),
                               worldtube_dz2 = 0.25, P = 4, L = 0)
        dom = block_domain(orb, cfg)

        # the rule itself: exactly the elements on the inner edge of the inner radial column
        hold = horizon_order_hold(dom, dom.mesh)
        @test length(hold) == count(el -> block_lattice_index(dom, el.block)[1] == 0,
                                    dom.mesh.elements)
        @test all(==( (true, true) ), values(hold))
        @test all(el -> block_lattice_index(dom, el.block)[1] == 0, keys(hold))

        fixed = solve_mode_amr(orb, cfg, m; dom = dom, iterations = 6,
                               horizon_hold = true)
        pf = horizon_axis_profile(fixed.solution)
        @printf("  [row 58] a = 0.998, ISCO, m = 17, nu = %.2f, HELD: ratio = %.4e, q = %+.3f, layer/plane = %.4f, ndofs = %d\n",
                horizon_branch_index(orb, m), pf.ratio, pf.exponent, pf.layer_ratio,
                fixed.solution.mesh.ndofs)
        @test pf.ratio < 50                    # measured 7.6 (was 1.9e4)
        @test abs(pf.exponent) < 1             # measured −0.3 (was −1.63)
        @test pf.layer_ratio <= 1.05           # measured 1.00 (was 1.3e2)

        # FALSIFIER — the DEFAULT sweep (`horizon_hold = false`, bitwise the pre-row-58
        # composition) still has the layer, so the gate above is not vacuous
        old = solve_mode_amr(orb, cfg, m; dom = dom, iterations = 6)   # the default
        po = horizon_axis_profile(old.solution)
        @printf("  [row 58] a = 0.998, ISCO, m = 17, DEFAULT: ratio = %.4e, q = %+.3f, layer/plane = %.4f, ndofs = %d\n",
                po.ratio, po.exponent, po.layer_ratio, old.solution.mesh.ndofs)
        @test po.layer_ratio > 10
        @test po.ratio > 1e3
    end
    # ---- (e) the refinement rule: structure, neutrality, preconditions --------------------
    #
    # No solve here: everything below is a MESH transformation, so it is cheap and it isolates
    # the rule from the physics. The default branch of `amr_step!` is literally
    # `h_refine(p_refine(mesh; hold = closed, increment), closed)` with
    # `closed = balanced_flags(mesh, flags)`, which is the DEFINITION of `hp_refine`; gate (e1)
    # pins that identity on the objects the sweep carries.
    @testset "horizon_order_hold + the default path" begin
        a, r0 = 0.5, 10.0
        orb = CircularOrbit(a, r0)
        mkcfg(; kw...) = crp_paper_config(; spin = a, orbital_radius = r0, m_mode = 2,
                                          inner_radius = orb.rp, outer_radius = 1.0e4,
                                          worldtube_radius = 4.8, worldtube_dz2 = 0.25,
                                          P = 4, L = 0, kw...)
        dom = block_domain(orb, mkcfg())
        mesh = dom.mesh
        isinner(el) = block_lattice_index(dom, el.block)[1] == 0 && el.segments[1].index == 0

        # (e0) on the initial lattice: exactly the inner radial column, both angular rows
        hold0 = horizon_order_hold(dom, mesh)
        @test length(hold0) == 2                       # 5 x 2 lattice, L = 0
        @test all(isinner, keys(hold0))
        @test all(==((true, true)), values(hold0))

        # the FULL-z lattice has three angular rows, hence three horizon elements
        domz = block_domain(orb, mkcfg(equatorial_symmetry = false))
        @test domz.chart === :z
        holdz = horizon_order_hold(domz, domz.mesh)
        @test length(holdz) == 3
        @test all(el -> block_lattice_index(domz, el.block)[1] == 0, keys(holdz))

        # (e1) THE BITWISE-NEUTRALITY GATE — the default path IS `hp_refine`, for raw flags and
        # for flags the 2:1 balance expands, at two increments
        punct = amr_flags(mesh, dom, orb, mkcfg())
        for flags in (punct, Dict(mesh.elements[end] => (true, true))), inc in (1, 2)
            closed = balanced_flags(mesh, flags)
            got = h_refine(p_refine(mesh; hold = closed, increment = inc), closed)
            ref = hp_refine(mesh, flags; increment = inc)
            @test got.elements == ref.elements
            @test got.orders == ref.orders
            @test got.offsets == ref.offsets
            @test got.ndofs == ref.ndofs
        end

        # (e2) under a RADIAL split of a horizon element only the INNER child keeps the horizon
        # face — the discriminating clause, which the unrefined lattice cannot exercise
        hel = first(sort(collect(keys(hold0)); by = el -> el.block))
        for dims in ((true, false), (false, true), (true, true))
            closed = balanced_flags(mesh, Dict(hel => dims))
            m2 = h_refine(p_refine(mesh; hold = closed, increment = 1), closed)
            h2 = horizon_order_hold(dom, m2)
            @test !isempty(h2)
            @test all(isinner, keys(h2))               # never an outer child, never column > 0
            @test all(==((true, true)), values(h2))
            kids = [el for el in m2.elements if el.block == hel.block]
            heldhere = [el for el in keys(h2) if el.block == hel.block]
            # the held set of that block IS its inner-radial-index set — balance-proof, and
            # strictly stronger than a count
            @test Set(heldhere) == Set(el for el in kids if isinner(el))
            # a RADIAL split makes an OUTER child, and that child is NOT held: this is the
            # clause the unrefined lattice cannot exercise, because there the two conditions
            # (inner lattice column, inner radial segment) are satisfied by the same element
            dims[1] && @test any(el -> !isinner(el), kids)
            # order inheritance, spectre's rule (UpdateAmrDecision.cpp:82-87): a split
            # dimension HOLDS its order, an unsplit one takes the p step, and every child of
            # this parent carries that same pair
            po = mesh.orders[mesh.index[hel]]
            want = (dims[1] ? po[1] : po[1] + 1, dims[2] ? po[2] : po[2] + 1)
            @test all(el -> m2.orders[m2.index[el]] == want, heldhere)
        end

        # (e3) a BALANCE-INDUCED split: flag the inner column's radial neighbour twice, so the
        # 2:1 closure splits the horizon element itself
        nb = first(el for el in mesh.elements if block_lattice_index(dom, el.block)[1] == 1)
        c1 = balanced_flags(mesh, Dict(nb => (true, true)))
        m3 = h_refine(p_refine(mesh; hold = c1, increment = 1), c1)
        nb2 = first(el for el in m3.elements
                    if block_lattice_index(dom, el.block)[1] == 1 && el.segments[1].index == 0)
        c2 = balanced_flags(m3, Dict(nb2 => (true, true)))
        m4 = h_refine(p_refine(m3; hold = c2, increment = 1), c2)
        h4 = horizon_order_hold(dom, m4)
        @test all(isinner, keys(h4))
        @test all(==((true, true)), values(h4))
        @test length(h4) >= length(hold0)

        # (e4) PRECONDITIONS. The horizon layer is the horizon layer only because the domain
        # pins its inner bound; `crp_paper_config` does not.
        @test_throws ArgumentError block_domain(orb, mkcfg(inner_radius = orb.rp * 1.05))
        # and a same-SHAPED mesh from a different BlockDomain is refused, not silently held
        orb2 = CircularOrbit(a, 12.0)
        dom2 = block_domain(orb2, crp_paper_config(; spin = a, orbital_radius = 12.0,
                                                   m_mode = 2, inner_radius = orb2.rp,
                                                   outer_radius = 1.0e4, worldtube_radius = 4.8,
                                                   worldtube_dz2 = 0.25, P = 4, L = 0))
        @test length(dom2.mesh.elements) == length(mesh.elements)
        @test_throws ArgumentError horizon_order_hold(dom, dom2.mesh)

        # (e5) the branch index is derived with `β_m` dropped, which holds only at H = −1
        @test_throws ArgumentError horizon_branch_index(orb, 17; H = 0)
        @test_throws ArgumentError horizon_branch_index(orb, 17; H = 1)
        @test horizon_branch_index(orb, 17; H = -1) == horizon_branch_index(orb, 17)
        # the SIGNED root: purely imaginary, sign = sign(Ω − Ω_H), modulus = the index
        g = ScalarSelfForce.gamma_r_m(orb, 17, orb.rp, -1)
        c = (orb.rp - orb.rm) / (orb.rp^2 + orb.a^2)
        @test real(g) == 0
        @test sign(imag(g / c)) == sign(orb.Omega - orb.a / (orb.rp^2 + orb.a^2))
        @test abs(g / c) ≈ horizon_branch_index(orb, 17) rtol = 1e-14
    end
end
