# Grounding for Radiation.Einstein — the linearized (Lichnerowicz) Einstein tensor G⁽¹⁾_μν[h] of the
# reconstructed coordinate h, and the vacuum residual gate (G-M4).
#
# (1) ANALYTIC VACUUM IDENTITY (the physics gate, corpus-free): for a CCK radiation-gauge reconstruction of a
#     VACUUM perturbation about Kerr, G⁽¹⁾_μν[h] ≡ 0. With the EXACT analytic h-derivatives (RTJet projection) +
#     analytic Christoffels, ‖G⁽¹⁾‖ vanishes to machine precision RELATIVE to ‖h‖ at interior BL points, in both
#     gauges. (Measured: ORG r=8 ‖G‖≈9e-17 vs ‖h‖≈0.08; IRG r=8 ‖G‖≈3e-12 vs ‖h‖≈274.)
# (2) DELEGATION IDENTITY (live, corpus-free): christoffel_BL and _dg_BL reproduce Spacetime.christoffel /
#     Spacetime._dg_BL exactly (they ARE the delegation).
# (3) FD CROSS-CHECK: the finite-difference residual path (deriv=:fd) is small (~1e-8) — an independent confirmation
#     that the vanishing is not an artifact of the analytic differentiation.
# The legacy MetricReconstruction.lin_einstein_residual likewise vanishes to ~1e-16 (machine); both reconstructions
# are vacuum to their solver precision, so the gate is the absolute vanishing, not a cross-fixture magnitude.

using Test
import MetricReconstruction.Radiation as R
import Spacetime

@testset "Radiation.Einstein (linearized vacuum residual G-M4)" begin
    a = 0.5; ω = 0.2; l = 2; m = 2; θ = π/3; t = 0.0; φ = 0.0
    w = R.weyl_amplitudes(:unit, ω, l, m, a)

    @testset "(1) analytic vacuum identity ‖G⁽¹⁾‖/‖h‖ → 0 (machine)" begin
        for (g, r) in ((:ORG, 8.0), (:ORG, 12.0), (:IRG, 8.0))
            H = R.h_coordinate(g, :BL, a, ω, l, m, r, θ, t, φ, w; bc=:Up)
            G = R.lin_einstein_residual(g, a, ω, l, m, r, θ, t, φ, w; bc=:Up, deriv=:analytic)
            scale = maximum(abs, H)
            @test maximum(abs, G) < 1e-10 * scale     # vacuum to machine precision relative to ‖h‖
        end
    end

    @testset "(2) delegation identity vs Spacetime (exact)" begin
        r = 8.0
        @test maximum(abs, R.christoffel_BL(a, r, θ) .- Spacetime.christoffel(a, r, θ; coords=:BL)) == 0
        gr_d, gθ_d = R._dg_BL(a, r, θ)
        gr_s, gθ_s = Spacetime._dg_BL(a, r, θ)
        @test maximum(abs, gr_d .- gr_s) == 0 && maximum(abs, gθ_d .- gθ_s) == 0
    end

    @testset "(3) FD residual cross-check is small (independent)" begin
        for (g, r) in ((:ORG, 8.0), (:ORG, 12.0))
            H = R.h_coordinate(g, :BL, a, ω, l, m, r, θ, t, φ, w; bc=:Up)
            Gf = R.lin_einstein_residual(g, a, ω, l, m, r, θ, t, φ, w; bc=:Up, deriv=:fd)
            @test maximum(abs, Gf) < 1e-5 * maximum(abs, H)
        end
    end
end
