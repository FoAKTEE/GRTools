# Grounding for Perturbation.Radial hypergeometric MST layer (hypergeom.jl).
# Oracle: ref_teuk_hyp.jl = BHPToolkit exact 2F1 (In) / U (Up) per-n term table (6 rows, n=-3..3)
# frozen as data. The legacy Teukolsky test asserts _h2f1_exact / _hu_exact against exactly this
# table at rtol=1e-8/atol=1e-9; reproduced here for the rewritten form. AHA: the MST param wrappers
# (_aF/_bF/_cF/_aU) and the _safe_2f1 |z|>6 / near-integer branch heuristic stay LOCAL; the generic
# 2F1 z→∞ / Pfaff-Maclaurin fallbacks delegate DOWN to BigSpecial.f21_zinf/f21_maclaurin (already
# characterization-gated to the Teukolsky 2F1 in BigSpecial's own core suite). Also: BigFloat-vs-
# Float64 self-consistency, the f21_maclaurin closed-form identity, and the _hu_array/_dhu_array
# recurrence-seed identity + live-Teukolsky cross-check.

include("ref_teuk_hyp.jl")  # const REF_HYP
include("fixtures/leg_radial_hypergeom.jl")  # const REF_HYP_LEG (baked Teukolsky._hu_array/_h2f1_exact/_hu_exact, Phase 7-B)
import BigSpecial

@testset "Radial: exact 2F1 / U MST terms vs BHPToolkit oracle (6 modes, n=-3..3)" begin
    for r in REF_HYP
        eps = 2 * r.w; kappa = sqrt(1 - r.a^2); tau = (eps - r.m * r.a) / kappa
        x = (1 + kappa - r.r) / (2 * kappa)                       # In coordinate
        zhat = eps * r.r / 2 - eps * (1 - kappa) / 2              # Up coordinate
        for n in -3:3
            j = n + 4
            h2 = Perturbation._h2f1_exact(n, r.s, r.nu_re, tau, eps, x)
            hu = Perturbation._hu_exact(n, r.s, r.nu_re, tau, eps, zhat)
            @test isapprox(real(h2), r.h2re[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(h2), r.h2im[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(real(hu), r.hure[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(hu), r.huim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: BigFloat ₂F₁ at |z|>1 (Pfaff+Maclaurin fallback) vs oracle + Float64" begin
    # The MST In-coordinate x=(1+κ−r)/(2κ) is < −1 for every REF_HYP row, so ₂F₁ takes the
    # analytic-continuation branch that throws on Complex{BigFloat} upstream → the local _safe_2f1
    # routes through BigSpecial.f21_maclaurin (Pfaff). Assert it (a) matches the oracle, (b)
    # reproduces the Float64 path.
    setprecision(BigFloat, 512) do
        nbig = 0
        for r in REF_HYP
            eps = 2 * r.w; kappa = sqrt(1 - r.a^2); tau = (eps - r.m * r.a) / kappa
            x = (1 + kappa - r.r) / (2 * kappa)
            @test x < -1
            xb = BigFloat(x); taub = BigFloat(tau); epsb = BigFloat(eps); nub = Complex{BigFloat}(r.nu_re)
            for n in -3:3
                j = n + 4
                h64  = Perturbation._h2f1_exact(n, r.s, r.nu_re, tau, eps, x)
                hbig = Perturbation._h2f1_exact(n, r.s, nub, taub, epsb, xb)
                @test hbig isa Complex{BigFloat}
                @test isapprox(real(hbig), r.h2re[j]; rtol=1e-8, atol=1e-9)
                @test isapprox(imag(hbig), r.h2im[j]; rtol=1e-8, atol=1e-9)
                @test isapprox(ComplexF64(hbig), h64; rtol=1e-10, atol=1e-12)
                nbig += 1
            end
        end
        @test nbig == 6 * 7
    end
    # the delegated Maclaurin primitive vs a closed form: ₂F₁(1,1,2;z) = −log(1−z)/z, full precision
    setprecision(BigFloat, 256) do
        for z in (BigFloat(-3)/10, BigFloat(-85)/100, BigFloat(2)/10)
            v = BigSpecial.f21_maclaurin(Complex{BigFloat}(1), Complex{BigFloat}(1), Complex{BigFloat}(2), complex(z))
            @test isapprox(real(v), -log(1 - z) / z; rtol=1e-50)
            @test abs(imag(v)) < 1e-60
        end
    end
    # _safe_2f1 is byte-identical to upstream for Float64 (never enters the fallback there)
    @test Perturbation._safe_2f1(1.0+0im, 0.5+0im, 2.0+0im, -3.5+0im) ===
          Perturbation.HypergeometricFunctions._₂F₁(1.0+0im, 0.5+0im, 2.0+0im, -3.5+0im)
end

@testset "Radial: _hu_array / _dhu_array recurrence-seed identity + live-Teukolsky cross-check" begin
    for r in REF_HYP
        eps = 2 * r.w; kappa = sqrt(1 - r.a^2); tau = (eps - r.m * r.a) / kappa
        zhat = eps * r.r / 2 - eps * (1 - kappa) / 2
        Nmax = 12
        HU  = Perturbation._hu_array(r.s, r.nu_re, tau, eps, zhat, Nmax)
        dHU = Perturbation._dhu_array(r.s, r.nu_re, tau, eps, zhat, HU, Nmax)
        # seeds exact (the recurrence is seeded by the exact HU[0],HU[1] / dHU[0],dHU[1])
        @test HU[0] == Perturbation._hu_exact(0, r.s, r.nu_re, tau, eps, zhat)
        @test HU[1] == Perturbation._hu_exact(1, r.s, r.nu_re, tau, eps, zhat)
        @test dHU[0] == Perturbation._dhu_exact(0, r.s, r.nu_re, tau, eps, zhat)
        @test dHU[1] == Perturbation._dhu_exact(1, r.s, r.nu_re, tau, eps, zhat)
        # the full recurrence array matches the (baked) legacy Teukolsky array bit-for-bit
        leg = REF_HYP_LEG[(r.s, r.l, r.m, r.a, r.w)]   # baked HU/dHU, index n+Nmax+1 for n=-Nmax:Nmax
        for n in -Nmax:Nmax
            @test HU[n] == leg.hu[n+Nmax+1]
            @test dHU[n] == leg.dhu[n+Nmax+1]
        end
    end
end

@testset "Radial: _safe_2f1 large-|z| z→∞ connection branch (|z|>6) vs BigSpecial + analytic" begin
    # The |z|≤6 oracle rows above exercise ONLY the Pfaff/Maclaurin route. This grounds the
    # OTHER branch of _safe_2f1 — the z→∞ connection (DLMF 15.8.2), reached for |z|>6 with a−b
    # away from an integer, delegated to BigSpecial.f21_zinf.
    a = 1.3 + 0.0im; b = 0.4 + 0.0im; c = 2.1 + 0.0im   # a−b = 0.9 ∉ ℤ (non-degenerate)

    # (1) f21_zinf vs HypergeometricFunctions._₂F₁ at Float64 |z|>6 (upstream succeeds there).
    for z in (-8.0 + 0im, 7.0 + 3im, -6.0 - 4im, 10.0 + 0.5im, 9.0 + 5im)
        @test abs(z) > 6
        ref = Perturbation.HypergeometricFunctions._₂F₁(a, b, c, z)
        @test isapprox(BigSpecial.f21_zinf(a, b, c, z), ref; rtol=1e-10)
    end

    # (2) _safe_2f1 at |z|>6 Complex{BigFloat} agrees with the delegated f21_zinf primitive to
    # full precision. (With HypergeometricFunctions ≥0.3.28 upstream now covers Complex{BigFloat}
    # for these args, so _safe_2f1 may return the upstream value rather than entering the catch;
    # either way the public result must equal the z→∞ primitive's value — that is what is grounded
    # here, route-agnostically — and match the independent Float64 upstream value.)
    setprecision(BigFloat, 256) do
        ab = Complex{BigFloat}(a); bb = Complex{BigFloat}(b); cb = Complex{BigFloat}(c)
        for z in (Complex{BigFloat}(9, 5), Complex{BigFloat}(-8, 0))
            v = Perturbation._safe_2f1(ab, bb, cb, z)
            @test v isa Complex{BigFloat}
            @test isapprox(v, BigSpecial.f21_zinf(ab, bb, cb, z); rtol=1e-50)
            @test isapprox(ComplexF64(v),
                Perturbation.HypergeometricFunctions._₂F₁(a, b, c, ComplexF64(z)); rtol=1e-10)
        end
        # (3) analytic cross-check: DLMF 15.8.2 (f21_zinf) vs the INDEPENDENT Pfaff continuation
        # (1−z)^(−a)·₂F₁(a,c−b,c; z/(z−1)) at Complex{BigFloat} — distinct analytic routes, agree
        # to full precision where Re z<1/2 (so the Pfaff |w|<1).
        for z in (Complex{BigFloat}(-8, 0), Complex{BigFloat}(-7, -3))
            @test real(z) < 0.5
            w = z / (z - 1)
            pfaff = (one(z) - z)^(-ab) * BigSpecial.f21_maclaurin(ab, cb - bb, cb, w)
            @test isapprox(BigSpecial.f21_zinf(ab, bb, cb, z), pfaff; rtol=1e-50)
        end
    end
end

@testset "Radial: new _h2f1_exact/_hu_exact == legacy Teukolsky (cross-check, baked exact, captured once)" begin
    for r in REF_HYP
        eps = 2 * r.w; kappa = sqrt(1 - r.a^2); tau = (eps - r.m * r.a) / kappa
        x = (1 + kappa - r.r) / (2 * kappa)
        zhat = eps * r.r / 2 - eps * (1 - kappa) / 2
        leg = REF_HYP_LEG[(r.s, r.l, r.m, r.a, r.w)]   # baked h2e/hue, index n+4 for n=-3:3
        for n in -3:3
            @test Perturbation._h2f1_exact(n, r.s, r.nu_re, tau, eps, x) == leg.h2e[n+4]
            @test Perturbation._hu_exact(n, r.s, r.nu_re, tau, eps, zhat) == leg.hue[n+4]
        end
    end
end
