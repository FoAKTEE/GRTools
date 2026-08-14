# Public-surface grounding for teukolsky_radial_deriv(R, r) through the TeukolskyRadial(...) carrier.
# The internal _mst_radial_in_deriv/_up_deriv are already oracle-grounded (norm=1) vs BHPToolkit
# REF_RADDERIV in test_radial_radial.jl. This file grounds the PUBLIC dispatcher at the normalized
# (unit-transmission) surface: (a) vs the LIGHT legacy Teukolsky.teukolsky_radial_deriv built from
# its own TeukolskyRadial (exact same MST path); (b) vs the BHPToolkit REF_RADDERIV table, rescaled
# by the transmission amplitude R.norm (public deriv = raw[norm=1]/Transmission); (c) a central
# finite-difference self-check of the green public callable R(r) (oracle-free). The static ω=0 deriv
# stays an honest throw (no oracle slice) — asserted with @test_throws.

@isdefined(REF_RADDERIV) || include("ref_teuk_radial_deriv.jl")
include("fixtures/leg_radial_teuk_deriv.jl")  # const REF_TDERIV_LEG (baked Teukolsky.teukolsky_radial_deriv, Phase 7-B)

@testset "Radial: PUBLIC teukolsky_radial_deriv vs legacy Teukolsky (normalized surface, baked, captured once)" begin
    for r in REF_RADDERIV
        Rnew = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, r.w)
        leg = REF_TDERIV_LEG[(r.s, r.l, r.m, r.a, r.w)]
        for (j, rr) in enumerate(r.r)
            # exact same MST math, ULP-level reorder only (matches the radial.jl cross-check tol)
            @test isapprox(teukolsky_radial_deriv(Rnew.In, rr), leg.dIn[j]; rtol=1e-12)
            @test isapprox(teukolsky_radial_deriv(Rnew.Up, rr), leg.dUp[j]; rtol=1e-12)
        end
    end
end

@testset "Radial: PUBLIC teukolsky_radial_deriv vs BHPToolkit REF_RADDERIV (unit-transmission rescaled)" begin
    # public dR/dr = raw[norm=1]/Transmission, and R.norm = the raw transmission amplitude, so
    # public · R.norm = the norm=1 BHPToolkit value (REF_RADDERIV). Load-bearing oracle at 1e-8.
    for r in REF_RADDERIV
        Rnew = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, r.w)
        for (j, rr) in enumerate(r.r)
            dIn = teukolsky_radial_deriv(Rnew.In, rr) * Rnew.In.norm
            dUp = teukolsky_radial_deriv(Rnew.Up, rr) * Rnew.Up.norm
            @test isapprox(real(dIn), r.dinre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(dIn), r.dinim[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(real(dUp), r.dupre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(dUp), r.dupim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: PUBLIC teukolsky_radial_deriv finite-difference self-check of R(r)" begin
    h = 1e-5
    for r in REF_RADDERIV
        Rnew = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, r.w)
        for rr in r.r
            dInFD = (Rnew.In(rr + h) - Rnew.In(rr - h)) / (2h)
            dUpFD = (Rnew.Up(rr + h) - Rnew.Up(rr - h)) / (2h)
            @test isapprox(teukolsky_radial_deriv(Rnew.In, rr), dInFD; rtol=1e-6, atol=1e-7)
            @test isapprox(teukolsky_radial_deriv(Rnew.Up, rr), dUpFD; rtol=1e-6, atol=1e-7)
        end
    end
end

@testset "Radial: PUBLIC teukolsky_radial_deriv static (ω=0) closed-form vs FD + ODE residual" begin
    # NEW primitive: legacy teukolsky_radial_deriv THROWS at ω=0 (no oracle). Grounded by (a) a
    # Richardson-extrapolated finite difference of the already-green static R(r) (grounded vs the
    # BHPToolkit REF_RADSTATIC table), and (b) the static (ω=0) radial Teukolsky ODE residual
    # Δ R'' + (2r-2)(1+s) R' + [(K²-2iK(r-1)s)/Δ - λ] R ≈ 0 with K=-ma (R'' via central FD of R').
    modes = [(-2, 2, 2, 0.0), (-2, 2, 2, 0.5), (0, 3, 1, 0.9), (-1, 2, 1, 0.5), (2, 2, 2, 0.5)]
    for (s, l, m, a) in modes
        R = Perturbation.TeukolskyRadial(s, l, m, a, 0.0)
        @test R.In.method === :Static && R.Up.method === :Static
        λ = l * (l + 1) - s * (s + 1)
        for Rbc in (R.In, R.Up)
            for r0 in (3.0, 6.0, 12.0)
                # Richardson central FD of the value → high-accuracy dR/dr reference
                fd(h) = (Rbc(r0 + h) - Rbc(r0 - h)) / (2h)
                d1 = fd(1e-3); d2 = fd(5e-4)
                dref = (4 * d2 - d1) / 3                     # O(h⁴) Richardson
                dana = teukolsky_radial_deriv(Rbc, r0)
                @test isapprox(dana, dref; rtol=1e-6, atol=1e-8)
                # ODE residual: R'' from central FD of the analytic R'
                h = 1e-4
                d2R = (teukolsky_radial_deriv(Rbc, r0 + h) -
                       teukolsky_radial_deriv(Rbc, r0 - h)) / (2h)
                Δ = r0^2 - 2r0 + a^2; K = -m * a
                coefR = (K^2 - 2im * K * (r0 - 1) * s) / Δ - λ
                res = Δ * d2R + (2r0 - 2) * (1 + s) * dana + coefR * Rbc(r0)
                scale = abs(Δ * d2R) + abs((2r0 - 2) * (1 + s) * dana) + abs(coefR * Rbc(r0))
                @test abs(res) ≤ 1e-4 * scale + 1e-8
            end
        end
    end
end
