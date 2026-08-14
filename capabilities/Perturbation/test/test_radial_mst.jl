# Grounding for Perturbation.Radial MST minimal-solution coefficients f_n (mst.jl).
# Oracle: ref_teuk_fn.jl = BHPToolkit f_n table (144 modes, n=-10..10) frozen as data. The legacy
# Teukolsky test asserts its _mst_fn_array against exactly this table at rtol=1e-8/atol=1e-10; we
# reproduce that for the rewritten form (which reuses NumericsCore.cf_convergent), plus the
# minimal-solution normalization f[0]==1, the cf_convergent golden-ratio identity, the recurrence
# residual self-consistency, and a live-Teukolsky cross-check.

include("ref_teuk_fn.jl")  # const REF_FN
include("fixtures/leg_radial_mst.jl")  # const REF_FN_LEG (baked Teukolsky._mst_fn_array, Phase 7-B)
import NumericsCore: cf_convergent

@testset "Radial: MST f_n coefficients vs BHPToolkit oracle (144 modes, n=-10..10)" begin
    for r in REF_FN
        r.ok == 1 || continue
        f, _, _, _ = Perturbation._mst_fn_array(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, 10)
        @test f[0] == 1                                  # minimal-solution normalization f_0 = 1
        for n in -10:10
            j = n + 11
            @test isapprox(real(f[n]), r.fre[j]; rtol=1e-8, atol=1e-10)
            @test isapprox(imag(f[n]), r.fim[j]; rtol=1e-8, atol=1e-10)
        end
    end
end

@testset "Radial: f_n recurrence residual (self-consistency)" begin
    # α_n f_{n+1} + β_n f_n + γ_n f_{n-1} = 0 at the renorm ν
    r = REF_FN[1]
    f, alpha, beta, gamma = Perturbation._mst_fn_array(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, 10)
    for n in -8:8
        res = alpha(n) * f[n+1] + beta(n) * f[n] + gamma(n) * f[n-1]
        @test abs(res) < 1e-9
    end
end

@testset "Radial: cf_convergent golden-ratio identity grounds the f_n ratios" begin
    # The f_n upward/downward ratios are NumericsCore.cf_convergent; all-ones K-form → 1/golden.
    @test isapprox(cf_convergent(k -> 1.0, k -> 1.0, 1), (sqrt(5) - 1) / 2; rtol=1e-12)
end

@testset "Radial: new f_n == legacy Teukolsky f_n (cross-check, baked, captured once)" begin
    for r in REF_FN
        r.ok == 1 || continue
        fnew, = Perturbation._mst_fn_array(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, 10)
        fleg = REF_FN_LEG[(r.s, r.l, r.m, r.a, r.w)]   # baked vector, n=-10..10 → index n+11
        for n in -10:10
            @test isapprox(fnew[n], fleg[n+11]; rtol=1e-12, atol=1e-14)
        end
    end
end

@testset "Radial: single-coefficient _mst_fn matches the array entry" begin
    r = REF_FN[1]
    f, = Perturbation._mst_fn_array(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, 5)
    for n in -5:5
        @test Perturbation._mst_fn(r.s, r.m, r.a, r.w, r.nu_re, r.lam_re, n) == f[n]
    end
end
