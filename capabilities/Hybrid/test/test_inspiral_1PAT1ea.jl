# Characterization of the 1PAT1ea / 1PAT1R 7-var spin+δm forcing models (Hybrid.Models, §4.9).
#
# GROUNDING (no WaSABI live dep — ANALYTIC closed forms + synthetic fixtures, reusing the
# already-grounded BindingEnergy.E_dm/dEdΩ and the F0/F1 channels grounded in test_inspiral_1PAT1):
#   * flux_LH = r0^{3/2} F_H, flux_sigma = F2σI+F2σH, flux_a = FaI+FaH: additive identities.
#   * forcing_F1s1/F1s2/F1extra/F1dm: the literal 1PAT1ea coupling rationals.
#   * omega_crit == min(1/6.26^{3/2}, Ωcrit(ν,χt1,χt2)) over a synthetic Interp3Lagrange grid.
#   * _spin_eom! (1PAT1ea additive): du == the literal channel combination; ν,M,χt2 frozen.
#   * _spin_eom! (1PAT1R ratio): du[1] == _domega_times_dEdΩ / dEdΩ / M, num rebuilt from the
#     literal 1PAT1R.m formula, den from the grounded dEdΩ.

using Test
using Hybrid.Models
using Hybrid.BindingEnergy: E_dm, dEdΩ

@testset "Models.inspiral_1PAT1ea" begin
    # --- base 1PAT1 synthetic data (same scheme as test_inspiral_1PAT1) -----------------
    fI(r) = 0.02 + 1.0e-3 * r - 5.0e-5 * r^2 + 1.0e-6 * r^3
    fH(r) = 0.005 - 2.0e-4 * r + 3.0e-6 * r^2
    dfI(r) = 1.0e-3 - 1.0e-4 * r + 3.0e-6 * r^2
    dfH(r) = -2.0e-4 + 6.0e-6 * r
    rs = collect(range(6.2, 32.0; length = 14))
    FEH = interp1_order8(rs, fH.(rs))
    FEI = interp1_order8(rs, fI.(rs))
    a2, b2 = -3.0, -3.0
    g2(r) = exp(a2) * r^b2
    F2 = interp1_order8(log.(rs), [complex(a2 + b2 * log(r), 0.0) for r in rs])
    zfun(x) = 0.5 + 0.3 * x - 0.1 * x^2
    xs = collect(range(0.02, 0.18; length = 12))
    z = interp1_order8(xs, zfun.(xs))
    base = Forcing1PAT1(FEH, FEI, F2, z)

    # --- spin flux tables (synthetic polynomials) ---------------------------------------
    aI(r) = 0.01 + 5.0e-4 * r; aH(r) = 0.002 - 1.0e-4 * r
    sI(r) = 0.03 - 8.0e-4 * r; sH(r) = 0.001 + 2.0e-5 * r
    FaI = interp1_order8(rs, aI.(rs)); FaH = interp1_order8(rs, aH.(rs))
    F2σI = interp1_order8(rs, sI.(rs)); F2σH = interp1_order8(rs, sH.(rs))
    f = Forcing1PAT1ea(base, F2σI, F2σH, FaI, FaH)

    pref(r0) = 3 * (1 - 3 / r0)^(3 / 2) * sqrt(1 / r0) / (1 - 6 / r0)

    @testset "spin flux identities" begin
        for r0 in (8.0, 12.5, 20.0)
            @test flux_LH(f, r0) ≈ r0^(3 / 2) * fH(r0) rtol = 1e-10
            @test flux_a(f, r0) ≈ aI(r0) + aH(r0) rtol = 1e-10
            @test flux_sigma(f, r0) ≈ sI(r0) + sH(r0) rtol = 1e-10
        end
    end

    @testset "spin/δm forcing channels (literal closed forms)" begin
        for r0 in (8.0, 12.5, 20.0)
            F1 = fI(r0) + fH(r0)
            back1 = (2 * sqrt(1 - 3 / r0) * (36 - 33r0 + 10r0^2)) / (r0^2 * (r0 - 6)^2)
            @test forcing_F1s1(f, r0) ≈ pref(r0) * (aI(r0) + aH(r0)) - back1 * F1 rtol = 1e-10
            back2 = (3 / r0^2) * ((1 - 3 / r0)^(3 / 2) * (5 - 12 / r0)) / (1 - 6 / r0)^2
            @test forcing_F1s2(f, r0) ≈ pref(r0) * (sI(r0) + sH(r0)) - back2 * F1 rtol = 1e-10
            @test forcing_F1extra(f, r0) ≈
                  -(1 / r0)^(3 / 2) * ((5 - 12 / r0) / (1 - 6 / r0)) * fH(r0) rtol = 1e-10
            t1 = pref(r0) * (sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * (dfI(r0) + dfH(r0)))
            t2 = (sqrt(1 - 3 / r0) * sqrt(1 / r0) * (3 * (7 - 6 / r0) / r0 - 2)) /
                 (1 - 6 / r0)^2 * F1
            @test forcing_F1dm(f, r0) ≈ t1 + t2 rtol = 1e-9
        end
    end

    # --- 3D stop interpolant + 7-var EOMs -----------------------------------------------
    νg = [0.05, 0.1, 0.15]; c1g = [-0.5, 0.0, 0.5]; c2g = [-0.5, 0.0, 0.5]
    Ωc(ν, χ1, χ2) = 0.05 + 0.1ν + 0.03χ1 + 0.02χ2 + 0.05ν * χ1
    grid = [((ν, χ1, χ2), Ωc(ν, χ1, χ2)) for ν in νg for χ1 in c1g for χ2 in c2g]
    itp3 = interp3_lagrange(grid)
    mea = Model1PAT1ea(f, itp3)
    mR = Model1PAT1R(f, itp3)
    floor1pa = 1 / 6.26^(3 / 2)

    @testset "omega_crit min(floor, Ωcrit)" begin
        for (ν, χ1, χ2) in ((0.1, 0.0, 0.0), (0.05, 0.3, -0.2), (0.15, 0.5, 0.4))
            @test omega_crit(mea, ν, χ1, χ2) ≈ min(floor1pa, Ωc(ν, χ1, χ2)) rtol = 1e-10
            @test omega_crit(mR, ν, χ1, χ2) ≈ min(floor1pa, Ωc(ν, χ1, χ2)) rtol = 1e-10
        end
    end

    @testset "_spin_eom! 1PAT1ea additive" begin
        ν, M, χt1, χt2, δm = 0.1, 1.0, 0.3, -0.2, 0.4
        Ω = 0.05; r0 = Ω^(-2 / 3)
        u = [Ω, 1.2, ν, M, χt1, χt2, δm]; du = zeros(7)
        _spin_eom!(du, u, mea)
        exp1 = (1 / M) * (ν * forcing_F0(f, r0) +
                          ν^2 * (forcing_F1(f, r0) + forcing_F1extra(f, r0) + δm * forcing_F1dm(f, r0)) +
                          ν * χt1 * forcing_F1s1(f, r0) + ν * χt2 * forcing_F1s2(f, r0))
        @test du[1] ≈ exp1 rtol = 0
        @test du[2] == Ω / M
        @test du[3] == 0.0 && du[4] == 0.0 && du[6] == 0.0
        @test du[5] ≈ (ν^2 / M) * flux_LH(f, r0) rtol = 0
        @test du[7] ≈ (ν / M) * fH(r0) rtol = 1e-12
    end

    @testset "_spin_eom! 1PAT1R ratio" begin
        ν, M, χt1, χt2, δm = 0.1, 1.0, 0.3, -0.2, 0.4
        Ω = 0.05; r0 = Ω^(-2 / 3); x = 1 / r0
        u = [Ω, 1.2, ν, M, χt1, χt2, δm]; du = zeros(7)
        _spin_eom!(du, u, mR)
        dm_term = δm * (sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * (dfI(r0) + dfH(r0)))
        num = -ν * flux_1SF(f, r0) +
              ν^2 * (-flux_2SF(f, r0) - dm_term - χt1 * flux_a(f, r0) / ν -
                     χt2 * flux_sigma(f, r0) / ν -
                     (2 * x^(5 / 2) * (-2 + 3x) * flux_LH(f, r0)) / (3 * (1 - 3x)^(3 / 2)) -
                     E_dm(x) * fH(r0))
        den = dEdΩ(z, r0, ν, χt1, χt2, δm)
        @test du[1] ≈ (1 / M) * num / den rtol = 1e-12
        @test du[5] ≈ (1 / M) * ν^2 * flux_LH(f, r0) rtol = 0
        @test du[7] ≈ (ν / M) * fH(r0) rtol = 1e-12
    end
end
