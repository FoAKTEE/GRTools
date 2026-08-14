# Grounding for Perturbation.Radial renorm-angmom ν (RenormalizedAngularMomentum + MST machinery).
# Oracle: ref_teuk_nu.jl = the BHPToolkit (SpinWeightedSpheroidalHarmonics`) ν table frozen as data
# (378 rows). The legacy Teukolsky test asserts its RenormalizedAngularMomentum against exactly this
# table at these tolerances; we reproduce that characterization for the rewritten Perturbation form
# (which reuses NumericsCore.cf_finite for the CF tails) + analytic identities + a live-Teukolsky
# cross-check that the new value == the legacy value to full precision.

include("ref_teuk_nu.jl")  # const REF_NU
include("fixtures/leg_renorm_angmom.jl")  # const REF_NU_LEG (baked Teukolsky.RenormalizedAngularMomentum, Phase 7-B)
import NumericsCore: cf_finite

@testset "Radial: RenormalizedAngularMomentum (FindRoot) real branch vs BHPToolkit oracle" begin
    nreal = 0
    for r in REF_NU
        (r.ok == 1 && abs(r.nu_im) < 1e-10) || continue
        nu = Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w, r.lam_re; method=:FindRoot)
        @test isapprox(nu, r.nu_re; rtol=1e-9, atol=1e-10)
        nreal += 1
    end
    @test nreal == 330   # guard against silently testing nothing
end

@testset "Radial: renorm-angmom identities + branches" begin
    # 5-arg auto-lambda (via Harmonics SWSH spectral eigenvalue) close to oracle (λ differs ~1e-11)
    for r in REF_NU[1:60]
        (r.ok == 1 && abs(r.nu_im) < 1e-10) || continue
        @test isapprox(Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w), r.nu_re; rtol=1e-6, atol=1e-7)
    end
    @test Perturbation.RenormalizedAngularMomentum(-2, 2, 2, 0.5, 0.0, 4.0) == 2.0  # static ω=0 → l
    @test Perturbation.RenormalizedAngularMomentum(-2, 1, 0, 0.5, 0.1, 0.0) == 0.0  # l<|s| → 0
    # :Series — closed-form O(ω⁴) seed; ≈ FindRoot at small ω (~1e-3, the series order)
    for r in REF_NU[1:60]
        (r.ok == 1 && abs(r.nu_im) < 1e-10 && r.w <= 0.1) || continue
        @test isapprox(Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w, r.lam_re; method=:Series),
                       r.nu_re; rtol=2e-3, atol=1e-4)
    end
    # :Monodromy (confluent-Heun monodromy) returns the same ν as :FindRoot
    @test isapprox(Perturbation.RenormalizedAngularMomentum(-2, 4, 4, 0.9, 0.1, 17.141170617488463; method=:Monodromy),
                   Perturbation.RenormalizedAngularMomentum(-2, 4, 4, 0.9, 0.1, 17.141170617488463; method=:FindRoot); rtol=1e-9)
    @test_throws ArgumentError Perturbation.RenormalizedAngularMomentum(-2, 3, 2, 0.5, 0.1, 1.0; method=:Bogus)
end

@testset "Radial: RenormalizedAngularMomentum complex-ν branch (½+iν_i / iν_i) vs oracle" begin
    ncx = 0
    for r in REF_NU
        (r.ok == 1 && abs(r.nu_im) >= 1e-10) || continue
        nu = Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w, r.lam_re; method=:FindRoot)
        @test isapprox(nu, complex(r.nu_re, r.nu_im); rtol=1e-8, atol=1e-9)
        ncx += 1
    end
    @test ncx == 6
end

@testset "Radial: RenormalizedAngularMomentum :Monodromy vs oracle (329 real + 6 cx)" begin
    # Matches FindRoot on every ok:1 row EXCEPT the degenerate boundary mode (-1,3,-3,0.9,0.5)
    # where Re[Cos2πν]<-1 and BHPToolkit's own two methods disagree.
    nreal = 0; ncx = 0
    for r in REF_NU
        r.ok == 1 || continue
        (r.s == -1 && r.l == 3 && r.m == -3 && r.a == 0.9 && r.w == 0.5) && continue
        nu = Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w, r.lam_re; method=:Monodromy)
        if abs(r.nu_im) < 1e-10
            @test isapprox(nu, r.nu_re; rtol=1e-9, atol=1e-10); nreal += 1
        else
            @test isapprox(nu, complex(r.nu_re, r.nu_im); rtol=1e-9, atol=1e-10); ncx += 1
        end
    end
    @test nreal == 329
    @test ncx == 6
end

@testset "Radial: cf_finite golden-ratio identity grounds the ν characteristic CF tails" begin
    # The MST characteristic CF tails are NumericsCore.cf_finite; the all-ones K-form → 1/golden.
    @test isapprox(cf_finite(k -> 1, k -> 1, 1, 60), (sqrt(5) - 1) / 2; rtol=1e-10)
end

@testset "Radial: new ν == legacy Teukolsky ν (full precision cross-check, baked, captured once)" begin
    # Direct characterization vs the (baked) LIGHT oracle on a representative slice of each branch.
    for r in REF_NU
        r.ok == 1 || continue
        (r.s == -1 && r.l == 3 && r.m == -3 && r.a == 0.9 && r.w == 0.5) && continue
        legpair = REF_NU_LEG[(r.s, r.l, r.m, r.a, r.w)]
        for (meth, leg) in ((:FindRoot, legpair.fr), (:Monodromy, legpair.mo))
            new = Perturbation.RenormalizedAngularMomentum(r.s, r.l, r.m, r.a, r.w, r.lam_re; method=meth)
            @test isapprox(new, leg; rtol=1e-12, atol=1e-13)
        end
    end
end
