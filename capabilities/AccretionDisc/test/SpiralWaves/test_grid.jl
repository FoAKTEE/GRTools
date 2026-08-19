# SpiralWaves/test_grid.jl — EGD::sw-grid (arXiv:2608.17003 App. C, appendix.tex:53-128).
# Tests FIRST for `grid.jl`: the corotation radius (61), the three-block adaptive radial grid
# (left WKB | sinh-clustered corotation (C4)-(C5) | right WKB), the fiducial / criterion-selected
# domain (C1)-(C3), the WKB wavenumber and its trust floor (C2),(C6), the trial steps (C7)-(C8),
# smoothMin (C9)-(C10) and the cubic-spline node placement. GROUNDING: internal exactness (the
# central block reproduces (C4) to 1e-14; the WKB blocks land EXACTLY on the domain/central edges
# with the realised n_W ∈ [n_L, n_L+1] (C7)); the (C3) phase rule re-integrated on the finished grid;
# limits of smoothMin; the m = 1 vs m = 100 localisation described in A:69-71; BigFloat inherits.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using Spacetime: kerr_isco

@testset "SpiralWaves.grid (arXiv:2608.17003 App. C)" begin

    bg = DiscBackground(χ = 0.0, aEMRI = 15.0, h0 = 0.05)                 # Figs 1-3 style disc, 𝔭=𝔮=0
    bgp = DiscBackground(χ = 0.9, aEMRI = 15.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)

    @testset "smoothMin (C9)-(C10)" begin
        @test smoothmin(2.0, 2.0) == 2.0
        @test smoothmin(1.0, 10.0) == smoothmin(10.0, 1.0)                   # symmetric
        @test 1.0 ≤ smoothmin(1.0, 10.0) ≤ 1.01                               # soft minimum, sharpness 8
        # NOT a hard minimum for widely separated arguments: smoothMin(a, b) → a + b e^{−8} for b ≫ a
        @test smoothmin(1.0, 1000.0) ≈ 1 + 1000 * exp(-8) rtol = 1e-2
        @test smoothmin(1.0, 1e6) ≈ 1e6 * exp(-8) rtol = 1e-2
        @test min(3.0, 5.0) ≤ smoothmin(3.0, 5.0) ≤ max(3.0, 5.0)
        # closed form check at one point: f = −max/8; (a e^{a/f} + b e^{b/f})/(e^{a/f}+e^{b/f})
        a, b = 3.0, 5.0; f = -5.0 / 8
        @test smoothmin(a, b) ≈ (a * exp(a / f) + b * exp(b / f)) / (exp(a / f) + exp(b / f)) rtol = 1e-15
        @test smoothmin(big"3", big"5") isa BigFloat
    end

    @testset "corotation radius (61): ω_d(r_c) = 0" begin
        rc = corotation_radius(bg, 3)
        @test rc ≈ 15.0 rtol = 1e-13                                           # 𝔭 = 𝔮 = 0 ⇒ r_c = a• exactly
        @test abs(omega_d(bg, 3, rc)) ≤ 1e-15
        @test corotation_radius(bg, 1) == corotation_radius(bg, 100)           # m-independent
        rcp = corotation_radius(bgp, 5)
        @test abs(omega_d(bgp, 5, rcp)) ≤ 1e-14 * 5 * Omega_emri(bgp)
        @test Omega(bgp, rcp) ≈ Omega_emri(bgp) rtol = 1e-13
        @test rcp < bgp.aEMRI                                                  # sub-Keplerian ⇒ inward shift
        @test kerr_isco(bgp.χ) < rcp
        # BigFloat inherits
        rcb = corotation_radius(DiscBackground(χ = big"0.9", aEMRI = 15, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1), 5)
        @test rcb isa BigFloat && abs(rcb - rcp) ≤ 1e-12
    end

    @testset "WKB wavenumber (C2): k² = (ω_d² − κ²)/c_s², trust floor (C6)" begin
        m = 20
        # at corotation ω_d = 0 ⇒ k² = −κ²/c_s² < 0 ⇒ evanescent, k = 0
        @test wkb_wavenumber2(bg, m, 15.0) < 0
        @test wkb_wavenumber(bg, m, 15.0) == 0
        # far outside the OLR k is real and equals √((ω_d²−κ²))/c_s
        r = 25.0
        k2 = (omega_d(bg, m, r)^2 - kappa2(bg, r)) / csound(bg, r)^2
        @test wkb_wavenumber2(bg, m, r) ≈ k2 rtol = 1e-15
        @test wkb_wavenumber(bg, m, r) ≈ sqrt(k2) rtol = 1e-15
        # trust floor k̄ = max(k, k_trust/(h0 r)): equals k where trusted, the floor where not
        @test wkb_wavenumber_floored(bg, m, 15.0) ≈ 10 / (0.05 * 15.0) rtol = 1e-15
        kk = wkb_wavenumber(bg, 100, 25.0)
        @test kk * 25.0 > 10 / 0.05
        @test wkb_wavenumber_floored(bg, 100, 25.0) == kk
    end

    @testset "corotation block (C4)-(C5): nodes exact, symmetric, count n_C" begin
        rc = 15.0; nC = 1001
        blk = corotation_block(rc, 1e-3 * rc, 4.0, nC)
        @test length(blk) == nC
        for (j, x) in enumerate(range(-1.0, 1.0, length = nC))
            @test blk[j] ≈ rc + 1e-3 * rc * sinh(4x) / sinh(4.0) rtol = 1e-14
        end
        @test isapprox(blk[1], rc - 1e-3 * rc; rtol = 1e-14) && isapprox(blk[end], rc + 1e-3 * rc; rtol = 1e-14)
        @test blk[(nC + 1) ÷ 2] == rc
        @test all(diff(blk) .> 0)
        # sinh clustering: minimum spacing at r_c, maximum at the edges, ratio cosh(s_c) (up to the
        # O(s_c Δx) discreteness of the first/last cell)
        d = diff(blk)
        @test argmin(d) ∈ (nC ÷ 2, nC ÷ 2 + 1)
        @test d[end] / d[nC ÷ 2] ≈ cosh(4.0) rtol = 1e-2
        Δx = 2 / (nC - 1)
        @test d[end] ≈ 1e-3 * rc * 4 * coth(4.0) * Δx rtol = 1e-2                  # edge spacing σ s_c coth(s_c) Δx
        @test d[nC ÷ 2] ≈ 1e-3 * rc * 4 / sinh(4.0) * Δx rtol = 1e-4               # centre spacing σ s_c Δx/sinh(s_c)
    end

    @testset "fiducial domain (C1) and full grid structure" begin
        g = build_grid(bg, 2; nL = 30, nR = 30, nC = 3000, extent = :fiducial)
        @test g isa SpiralGrid{Float64}
        @test g.r[1] == kerr_isco(0.0) == g.r_in
        @test g.r[end] == 5 * 15.0 == g.r_out
        @test all(diff(g.r) .> 0)                                              # strictly increasing
        @test g.nC == 3000 && length(g.blocks.C) == 3000
        @test g.r[g.blocks.C] == corotation_block(g.r_c, 1e-3 * g.r_c, 4.0, 3000)
        @test g.r_Cin == g.r[first(g.blocks.C)] && g.r_Cout == g.r[last(g.blocks.C)]
        @test length(g.r) == g.NL + g.nC + g.NR                                # shared edges counted once
        @test g.blocks.L == 1:(g.NL + 1) && g.blocks.R == (g.NL + g.nC):(g.NL + g.nC + g.NR)
        @test last(g.blocks.L) == first(g.blocks.C) && last(g.blocks.C) == first(g.blocks.R)
        @test 30 ≤ g.nW_L ≤ 31 && 30 ≤ g.nW_R ≤ 31                              # (C7) n_W ∈ [n_L, n_L+1]
        # the realised cell counts are the integer the n_W optimisation targets: N = ⌈N(n_L)⌉ with
        # N(n_W) = ∫ dr/Δr_trial(r; n_W) the continuous cell count of the region, and N(n_W*) == N
        regL = wkb_region(bg, 2, g, :L); regR = wkb_region(bg, 2, g, :R)
        @test g.NL == ceil(Int, cell_count(regL, 30.0)) && g.NR == ceil(Int, cell_count(regR, 30.0))
        @test cell_count(regL, g.nW_L) ≈ g.NL atol = 1e-8
        @test cell_count(regR, g.nW_R) ≈ g.NR atol = 1e-8
        @test cell_count(regL, 31.0) ≥ g.NL && cell_count(regR, 31.0) ≥ g.NR
        # the first WKB step off the central block is continuous with the central spacing (C8):
        # smoothMin(Δr_WKB, Δr_in) ≤ Δr_in ≈ the central edge spacing
        dC = g.r_Cin - g.r[first(g.blocks.C) - 1]
        @test dC ≤ 1.05 * (g.r[first(g.blocks.C) + 1] - g.r_Cin)
        dCr = g.r[last(g.blocks.C) + 1] - g.r_Cout
        @test dCr ≤ 1.05 * (g.r_Cout - g.r[last(g.blocks.C) - 1])
        # determinism
        g2 = build_grid(bg, 2; nL = 30, nR = 30, nC = 3000, extent = :fiducial)
        @test g2.r == g.r
    end

    @testset "reference triple (240, 240, 24000): exact counts, monotone, landed edges" begin
        g = build_grid(bgp, 10; nL = 240, nR = 240, nC = 24000)
        @test g.nC == 24000 && length(g.blocks.C) == 24000
        @test length(g.r) == g.NL + 24000 + g.NR
        @test all(diff(g.r) .> 0)
        @test g.r[1] == g.r_in && g.r[end] == g.r_out
        @test 240 ≤ g.nW_L ≤ 241 && 240 ≤ g.nW_R ≤ 241
        regL = wkb_region(bgp, 10, g, :L); regR = wkb_region(bgp, 10, g, :R)
        @test g.NL == ceil(Int, cell_count(regL, 240.0)) && g.NR == ceil(Int, cell_count(regR, 240.0))
        @test cell_count(regL, g.nW_L) ≈ g.NL atol = 1e-6
        @test cell_count(regR, g.nW_R) ≈ g.NR atol = 1e-6
        @test g.NL ≥ 240 * 30 && g.NR ≥ 240 * 30                               # ≥ n_W × 30 radians (C3)
        println("  (240,240,24000) m=10 grid: N_L=$(g.NL) n_C=$(g.nC) N_R=$(g.NR) total=$(length(g.r)) ",
                "domain [$(g.r_in), $(g.r_out)] r_c=$(g.r_c)")
    end

    @testset "adaptive extent (C2)-(C3): m = 1 spans the disc, m = 100 is localised" begin
        g1 = build_grid(bg, 1; nL = 30, nR = 30, nC = 3000)
        g100 = build_grid(bg, 100; nL = 30, nR = 30, nC = 3000)
        # m = 1, h0 = 0.05, a• = 15: no point satisfies k r > 10/h0 ⇒ fiducial boundaries retained
        @test isnan(g1.rWKB_L) && isnan(g1.rWKB_R)
        @test g1.r_in == kerr_isco(0.0) && g1.r_out == 75.0
        # m = 100: WKB-trusted on both sides, 30 radians of phase ⇒ narrow domain around the perturber
        @test !isnan(g100.rWKB_L) && !isnan(g100.rWKB_R)
        @test g100.r_in > kerr_isco(0.0) && g100.r_out < 75.0
        @test (g100.r_out - g100.r_in) < (75.0 - 6.0) / 3
        @test g100.r_in < g100.rWKB_L < g100.r_c < g100.rWKB_R < g100.r_out
        @test wkb_wavenumber(bg, 100, g100.rWKB_L) * g100.rWKB_L > 10 / 0.05
        @test wkb_wavenumber(bg, 100, g100.rWKB_R) * g100.rWKB_R > 10 / 0.05
        # (C3): the accumulated WKB phase between the domain edge and r_WKB is 30 radians (trapezoid
        # on the finished grid, which is much finer than the ancillary grid)
        function phase(g, lo, hi)
            idx = findall(x -> lo ≤ x ≤ hi, g.r)
            rr = g.r[idx]; kk = [wkb_wavenumber(bg, 100, x) for x in rr]
            sum((kk[2:end] .+ kk[1:end-1]) .* diff(rr)) / 2
        end
        @test phase(g100, g100.r_in, g100.rWKB_L) ≈ 30 rtol = 2e-2
        @test phase(g100, g100.rWKB_R, g100.r_out) ≈ 30 rtol = 2e-2
        @test length(g100.r) < length(g1.r)
        println("  m=1: N=$(length(g1.r)) domain [$(g1.r_in), $(g1.r_out)];  m=100: N=$(length(g100.r)) ",
                "domain [$(g100.r_in), $(g100.r_out)], r_WKB=($(g100.rWKB_L), $(g100.rWKB_R))")
    end

    @testset "BigFloat inherits (small grid)" begin
        bgb = DiscBackground(χ = big"0.5", aEMRI = 12, h0 = 0.1, 𝔭 = 1)
        g = build_grid(bgb, 8; nL = 4, nR = 4, nC = 201, n_anc = 400, extent = :fiducial)
        @test g isa SpiralGrid{BigFloat} && eltype(g.r) == BigFloat
        @test all(diff(g.r) .> 0)
        @test g.r[1] == kerr_isco(big"0.5") && g.r[end] == 60
        @test length(g.r) == g.NL + 201 + g.NR
        @test big"4" ≤ g.nW_L ≤ big"5"
    end
end
