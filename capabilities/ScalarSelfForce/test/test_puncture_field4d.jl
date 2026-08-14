# test_puncture_field4d.jl — ESSF node `puncture-field-4d` (schedule W3-D), written BEFORE
# the implementation (TDD). Gates (GLOBAL_DAG #33 oracle column + node contract):
#   * characterization vs the BAKED effsource_PhiS 4D oracle (test/fixtures/effsource/
#     phis4d_elliptic.jl — captured from the SAME kerr-circular.c whose
#     effsource_set_particle block is transcribed; phis4d_comoving.jl cross-checks the
#     identical 4D field compiled into the comoving variant), rtol 1e-10, ALL 6 frozen
#     configs (incl. a = -0.5 and a = -0.998) × npts points × dphi ∈ {0, 0.5, 2.0}.
#     The W2-D near-worldline caveat (effsource_configs.jl header) concerns the calc_m
#     direct SOURCE only — the header states no reduced tolerance for PhiS values, so
#     the 1e-10 gate applies at every fixture point. Evaluation is at the capture
#     driver's EFFECTIVE offsets fl(r0+dr)−r0 etc. (see the in-testset note);
#   * A_ijk = 0 for odd j or k on circular equatorial orbits (main.tex:456-457):
#     structural — the coefficient name set is EXACTLY {A_ijk : j, k even,
#     6 ≤ i+j+k ≤ 9} (50 names, the C statics of kerr-circular.c:18) — and functional —
#     Φ^P is even under dtheta → -dtheta and under dphi → -dphi;
#   * ϱ²/z_c² scale literals vs main.tex:445-455 formulas, vs the solid CircularOrbit
#     fields, and vs the C-side alpha20/alpha02/beta baked in EFFSOURCE_CONFIGS
#     (kerr-circular.c:844-847), rtol 1e-13;
#   * singularity scaling Φ^P ~ 1/λ as λ → 0 along random-ish directions (main.tex:422):
#     log-log slope −1 within 2%;
#   * BigFloat(256) smoke at one point vs Float64 (rtol 1e-12) + type inheritance.
using Test
using ScalarSelfForce

# fixture constants may already be loaded by test_effsource_fixtures.jl under Pkg.test;
# guard the includes so this file also runs standalone.
if !@isdefined(EFFSOURCE_CONFIGS)
    include(joinpath(@__DIR__, "fixtures", "effsource", "effsource_configs.jl"))
end
if !@isdefined(EFFSOURCE_PHIS4D_ELLIPTIC)
    include(joinpath(@__DIR__, "fixtures", "effsource", "phis4d_elliptic.jl"))
end
if !@isdefined(EFFSOURCE_PHIS4D_COMOVING)
    include(joinpath(@__DIR__, "fixtures", "effsource", "phis4d_comoving.jl"))
end

@testset "Puncture: 4D Detweiler-Whiting Φ^P (puncture-field-4d)" begin

    punctures = Dict(c.idx => Puncture(CircularOrbit(c.a, c.r0)) for c in EFFSOURCE_CONFIGS)

    @testset "characterization vs baked effsource_PhiS (rtol 1e-10, 6 configs)" begin
        # The capture driver (generator/capture.c:290-296) evaluates at the ABSOLUTE point
        # x = (r0 + dr, π/2 + dθ, φ_p + dφ), so the C code sees the EFFECTIVE offsets
        # fl(r0+dr) − r0 etc. At |dr| = 1e-6 M the binade rounding of r0 + dr shifts the
        # offset by up to ~7.5e-10 RELATIVE, and Φ^P ~ 1/λ carries exactly that into the
        # value (diagnosed against BigFloat(256): at the shifted points the C value is off
        # truth by 7.5e-10 while ours at the recorded dr is ≤ 7.3e-15 — see the node
        # evidence). Characterize at the driver's effective offsets: worst rel dev then
        # collapses from 7.5e-10 to 1.6e-13. The 1e-10 gate itself is untouched.
        r0s = Dict(c.idx => c.r0 for c in EFFSOURCE_CONFIGS)
        for (label, fixture) in (("elliptic", EFFSOURCE_PHIS4D_ELLIPTIC),
                                 ("comoving", EFFSOURCE_PHIS4D_COMOVING))
            worst = 0.0
            for col in axes(fixture, 2)
                cfg = Int(fixture[1, col])
                dr, dtheta, dphi, oracle = fixture[3, col], fixture[4, col],
                                           fixture[5, col], fixture[6, col]
                r0 = r0s[cfg]
                dr_eff = (r0 + dr) - r0                                  # capture.c:290
                dth_eff = (pi / 2 + dtheta) - pi / 2                     # capture.c:290
                dph_eff = (EFFSOURCE_PHIP + dphi) - EFFSOURCE_PHIP       # capture.c:295
                ours = puncture_phi(punctures[cfg], dr_eff, dth_eff, dph_eff)
                @test isapprox(ours, oracle; rtol = 1e-10)
                worst = max(worst, abs(ours - oracle) / abs(oracle))
            end
            println("puncture-field-4d characterization [", label, "]: ",
                    size(fixture, 2), " points, worst rel dev = ", worst)
        end
    end

    @testset "A_ijk structure: even j,k only, 6 ≤ i+j+k ≤ 9 (main.tex:456-457)" begin
        anames = sort([String(f) for f in fieldnames(Puncture) if startswith(String(f), "A")])
        expected = sort(["A$i$j$k" for i in 0:9 for j in 0:2:8 for k in 0:2:8
                         if 6 <= i + j + k <= 9])
        @test length(expected) == 50            # the 50 C statics, kerr-circular.c:18
        @test anames == expected
    end

    @testset "Φ^P even under dtheta and dphi reflection (odd-j/odd-k vanishing)" begin
        for cfg in (2, 3, 6)                    # a = 0.5, 0.9, -0.998
            p = punctures[cfg]
            for (dr, dtheta, dphi) in ((0.01, 0.007, 0.3), (-0.05, 0.02, 1.1),
                                       (0.3, -0.15, 2.5), (1.0, 0.4, 0.05))
                v = puncture_phi(p, dr, dtheta, dphi)
                @test puncture_phi(p, dr, -dtheta, dphi) == v
                @test puncture_phi(p, dr, dtheta, -dphi) == v
                @test puncture_phi(p, dr, -dtheta, -dphi) == v
            end
            # dr enters with ODD powers too: no even symmetry in dr (guards against a
            # transcription that accidentally symmetrizes everything)
            @test puncture_phi(p, 0.1, 0.0, 0.5) != puncture_phi(p, -0.1, 0.0, 0.5)
        end
    end

    @testset "ϱ²/z_c² scale literals (main.tex:445-455) vs CircularOrbit vs C, rtol 1e-13" begin
        for c in EFFSOURCE_CONFIGS
            o = punctures[c.idx].orbit
            a, r0 = c.a, c.r0
            # main.tex:452  ϱ² = r_p²/(r_p² − 2M r_p + a²)·Δr² + r_p²·Δθ²   (M = 1)
            @test isapprox(o.rho2_dr2, r0^2 / (r0^2 - 2 * r0 + a^2); rtol = 1e-13)
            @test isapprox(o.rho2_dtheta2, r0^2; rtol = 1e-13)
            # main.tex:453  z_c² = 4(L² + r_p² + a² + 2M a²/r_p)
            @test isapprox(o.zc2, 4 * (o.Lz^2 + r0^2 + a^2 + 2 * a^2 / r0); rtol = 1e-13)
            # C-side R² pieces (kerr-circular.c:844-847), baked into the fixture configs:
            # alpha20 = r²/(a²+r(r−2M)), alpha02 = r², beta = 4(L²+r²+a²(r+2M)/r)
            @test isapprox(o.rho2_dr2, c.alpha20; rtol = 1e-13)
            @test isapprox(o.rho2_dtheta2, c.alpha02; rtol = 1e-13)
            @test isapprox(o.zc2, c.beta; rtol = 1e-13)
        end
    end

    @testset "Φ^P ~ 1/λ as λ → 0 (log-slope −1 within 2%)" begin
        # deterministic random-ish direction cosines (no Random dep in the test target)
        dirs = ((0.62, 0.53, 0.58), (-0.41, 0.27, -0.87), (0.11, -0.90, 0.42))
        for cfg in (1, 4, 5), n in dirs         # a = 0, 0.998@isco, -0.5
            p = punctures[cfg]
            lams = [1e-3 / 2.0^k for k in 0:10]
            vals = [abs(puncture_phi(p, l * n[1], l * n[2], l * n[3])) for l in lams]
            slope = (log(vals[end]) - log(vals[1])) / (log(lams[end]) - log(lams[1]))
            @test abs(slope + 1) < 0.02
        end
    end

    @testset "BigFloat(256) smoke vs Float64 (rtol 1e-12) + type inheritance" begin
        c = EFFSOURCE_CONFIGS[2]                # a = 0.5, r0 = 10
        p64 = punctures[c.idx]
        @test p64 isa Puncture{Float64}
        v64 = puncture_phi(p64, 0.01, 0.005, 0.5)
        @test v64 isa Float64
        setprecision(BigFloat, 256) do
            pb = Puncture(CircularOrbit(big(c.a), big(c.r0)))
            @test pb isa Puncture{BigFloat}
            vb = puncture_phi(pb, big(0.01), big(0.005), big(0.5))  # the SAME binary64 point
            @test vb isa BigFloat
            @test isapprox(v64, Float64(vb); rtol = 1e-12)
        end
        # promotion smoke: BigFloat point through a Float64 table stays BigFloat
        @test puncture_phi(p64, big(1) / 100, big(5) / 1000, big(1) / 2) isa BigFloat
    end
end
