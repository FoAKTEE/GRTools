# Characterization of the first-law binding-energy closed forms (Hybrid.BindingEnergy, §4.9).
#
# GROUNDING (no WaSABI live dep — pure ANALYTIC):
#   * E_dm / E_dm_deriv: literal closed form -(y/3)(1-6y)/(1-3y)^{3/2}, with the analytic
#     derivative cross-checked against a central finite difference of E_dm (~1e-7).
#   * E_FLx / E_FLx_deriv: feed a KNOWN low-degree polynomial as the redshift Interp1 z. The
#     order-8 interpolant reproduces z, z', z'' exactly (deg ≤ 8), so the whole expression is
#     hand-computable from the docstring closed form; E_FLx_deriv is independently checked by a
#     central finite difference of E_FLx.
#   * dEdΩ: geodesic + secondary-spin terms vs the literal formula (isolated by ν = 0); the SF
#     channel is cross-checked via (dEdΩ - geo - spin)/(ν dx/dΩ) == E_FLx_deriv + δm1 E_dm_deriv.

using Test
using Hybrid.BindingEnergy
using Hybrid.Models: interp1_order8

@testset "BindingEnergy" begin

    # --- E_dm / E_dm_deriv: literal algebra + FD cross-check -------------------------------
    @testset "E_dm secondary-mass correction" begin
        Edm_ref(y) = -(y / 3) * (1 - 6y) / (1 - 3y)^(3 / 2)
        for y in (0.03, 0.08, 0.12, 0.18)
            @test E_dm(y) ≈ Edm_ref(y) rtol=0
            h = 1e-6
            fd = (E_dm(y + h) - E_dm(y - h)) / (2h)
            @test E_dm_deriv(y) ≈ fd rtol=1e-7
        end
    end

    # --- E_FLx / E_FLx_deriv with a known polynomial redshift -------------------------------
    @testset "E_FLx flux function (polynomial redshift)" begin
        # z(x) degree 3 ≤ 8 ⇒ order-8 Interp1 reproduces z, z', z'' exactly.
        zc = (1.0, 2.0, -0.5, 0.1)
        z(x)  = zc[1] + zc[2] * x + zc[3] * x^2 + zc[4] * x^3
        zp(x) = zc[2] + 2zc[3] * x + 3zc[4] * x^2
        xs = collect(range(0.02, 0.22; length = 12))
        zitp = interp1_order8(xs, z.(xs))

        # EFLx[x] closed form (1PAT1.m L14), recomputed independently here.
        EFLx_ref(x) = 0.5 * z(x) - (x / 3) * zp(x) +
                      (-1 + sqrt(1 - 3x) + (x / 6) * (7 - 24x) / (1 - 3x)^(3 / 2))
        for x in (0.05, 0.10, 0.15)
            @test E_FLx(zitp, x) ≈ EFLx_ref(x) rtol=1e-10
            # E_FLx_deriv independently grounded by a central FD of E_FLx.
            h = 1e-6
            fd = (E_FLx(zitp, x + h) - E_FLx(zitp, x - h)) / (2h)
            @test E_FLx_deriv(zitp, x) ≈ fd rtol=1e-6
        end
    end

    # --- dEdΩ: geodesic+spin literal, SF channel cross-check --------------------------------
    @testset "dEdΩ flux per unit frequency" begin
        zc = (1.0, 2.0, -0.5, 0.1)
        z(x) = zc[1] + zc[2] * x + zc[3] * x^2 + zc[4] * x^3
        xs = collect(range(0.02, 0.22; length = 12))
        zitp = interp1_order8(xs, z.(xs))

        # geodesic + secondary-spin literal closed form (1PAT1R.m L84-87, M=1).
        function geo_spin(r0, at1, at2)
            x = 1 / r0
            geo = (-1 + 6x) / (3 * (1 - 3x)^(3 / 2) * sqrt(x))
            spin = (at2 * x * (-5 + 12x)) / (3 * (1 - 3x)^(3 / 2)) -
                   (2 * at1 * x * (10 - 33x + 36x^2)) / (9 * (1 - 3x)^(5 / 2))
            return geo + spin
        end

        r0, at1, at2 = 12.0, 0.3, -0.2
        # ν = 0 isolates the geodesic + spin terms.
        @test dEdΩ(zitp, r0, 0.0, at1, at2, 0.0) ≈ geo_spin(r0, at1, at2) rtol=1e-12

        # SF channel: (dEdΩ - geo - spin)/(ν dx/dΩ) == E_FLx_deriv + δm1 E_dm_deriv.
        ν, δm1 = 1.0e-3, 0.5
        x = 1 / r0
        # Ω = x^{3/2} (M=1) ⇒ dx/dΩ = (2/3) x^{-1/2} = 2/(3√x) (analytic Jacobian).
        dxdΩ = 2 / (3 * sqrt(x))
        full = dEdΩ(zitp, r0, ν, at1, at2, δm1)
        sf = (full - geo_spin(r0, at1, at2)) / (ν * dxdΩ)
        @test sf ≈ E_FLx_deriv(zitp, x) + δm1 * E_dm_deriv(x) rtol=1e-9
    end

end
