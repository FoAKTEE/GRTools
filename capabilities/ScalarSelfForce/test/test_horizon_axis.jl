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
end
