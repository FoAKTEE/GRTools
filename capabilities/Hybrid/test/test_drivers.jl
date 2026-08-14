# Characterization of the N-var inspiral ODE driver (Hybrid.Drivers, §4.9).
#
# GROUNDING (no live WaSABI / no OrdinaryDiffEq monolith — Vern9 vs analytic ODE solutions +
# rhs collocation + serial==threaded bit-identical parity):
#   * 2-var (Ω,φ): a synthetic CONSTANT-forcing model has the analytic solution Ω(t)=Ω0+Ct,
#     φ(t)=φ0+(Ω0 t + C t²/2)/M, terminating at the literal stop_r0 crossing; Vern9 must match
#     at the upstream reltol/abstol=1e-10.
#   * Kerr (ω,φ): termination at the (grounded) omega_crit; dω/dt == Fω/M by central-difference
#     collocation against the grounded forcing_Fω.
#   * 7-var spin/δm: termination at min(floor,Ωcrit); frozen ν,M,χt2; χt1' and δm' collocate to
#     (ν²/M)F_LH and (ν/M)F_H; serial==threaded bit-identical.

using Test
using Hybrid.Drivers
using Hybrid.Models
import Hybrid.Models: omega_dot, stop_r0

# Synthetic constant-forcing 2-var model: dΩ/dt = C, terminate at r0 = rstop.
struct ConstForcing
    C::Float64
    rstop::Float64
end
omega_dot(m::ConstForcing, r0::Real, ν::Real, M::Real) = m.C
stop_r0(m::ConstForcing, ν::Real) = m.rstop

@testset "Drivers" begin
    @testset "2-var Vern9 vs analytic constant-forcing" begin
        C, rstop, r0_init, ν, M, φ0 = 1.0e-3, 6.5, 10.0, 1.0e-3, 1.0, 0.4
        m = ConstForcing(C, rstop)
        s = int_inspiral(m, r0_init, ν, M; φ0 = φ0)
        Ω0 = r0_init^(-3 / 2)
        tstar = (rstop^(-3 / 2) - Ω0) / C
        @test duration(s) ≈ tstar rtol = 1e-7
        @test r0_at(s, duration(s)) ≈ rstop rtol = 1e-7
        for t in (0.2 * tstar, 0.5 * tstar, 0.9 * tstar)
            @test omega_at(s, t) ≈ Ω0 + C * t rtol = 1e-8
            @test phi_at(s, t) ≈ φ0 + (Ω0 * t + C * t^2 / 2) / M rtol = 1e-8
        end
    end

    @testset "2-var batch serial == threaded (bit-identical)" begin
        m = ConstForcing(1.0e-3, 6.5)
        ics = [(9.0, 1.0e-3, 1.0), (10.0, 2.0e-3, 1.0), (12.0, 1.0e-3, 1.0), (8.0, 3.0e-3, 1.0)]
        ser = int_inspiral_batch(m, ics; threaded = false)
        thr = int_inspiral_batch(m, ics; threaded = true)
        for i in eachindex(ics)
            @test duration(ser[i]) == duration(thr[i])
            @test omega_at(ser[i], 0.5 * duration(ser[i])) == omega_at(thr[i], 0.5 * duration(thr[i]))
        end
    end

    # --- Kerr (ω,φ) driver ---------------------------------------------------------------
    @testset "Kerr driver: termination + rhs collocation" begin
        Xn = log.(1 .- [0.1, 0.3, 0.5, 0.7]); Yn = [0.2, 0.4, 0.6, 0.8, 1.0]
        PI(X, Y) = 2.0 + 0.5X + 0.2Y; PH(X, Y) = 2.0 - 0.3X + 0.4Y
        gI = [((X, Y), PI(X, Y)) for X in Xn for Y in Yn]
        gH = [((X, Y), PH(X, Y)) for X in Xn for Y in Yn]
        model = Forcing0PAKerr(interp2_hermite(gI), interp2_hermite(gH))
        ν, M, χ1, ω0 = 1.0e-3, 1.0, 0.1, 0.05
        @test forcing_Fω(model, ω0, ν, χ1) > 0          # chirps up ⇒ terminates
        s = int_inspiral(model, ω0, ν, M, χ1)
        @test omega_at(s, duration(s)) ≈ omega_crit(model, χ1) rtol = 1e-7
        td = duration(s)
        for t in (0.3 * td, 0.6 * td)
            h = 1.0e-5 * td
            dωdt = (omega_at(s, t + h) - omega_at(s, t - h)) / (2h)
            @test dωdt ≈ forcing_Fω(model, omega_at(s, t), ν, χ1) / M rtol = 1e-5
        end
    end

    # --- 7-var spin/δm driver ------------------------------------------------------------
    @testset "7-var spin driver: termination + frozen/evolving collocation + parity" begin
        rs = collect(range(6.2, 32.0; length = 14))
        fI(r) = 0.02 + 1.0e-3 * r - 5.0e-5 * r^2 + 1.0e-6 * r^3
        fH(r) = 0.005 - 2.0e-4 * r + 3.0e-6 * r^2
        a2, b2 = -3.0, -3.0
        FEH = interp1_order8(rs, fH.(rs)); FEI = interp1_order8(rs, fI.(rs))
        F2 = interp1_order8(log.(rs), [complex(a2 + b2 * log(r), 0.0) for r in rs])
        z = interp1_order8(collect(range(0.02, 0.18; length = 12)),
                           [0.5 + 0.3x - 0.1x^2 for x in range(0.02, 0.18; length = 12)])
        base = Forcing1PAT1(FEH, FEI, F2, z)
        FaI = interp1_order8(rs, [0.01 + 5e-4r for r in rs]); FaH = interp1_order8(rs, [0.002 - 1e-4r for r in rs])
        F2σI = interp1_order8(rs, [0.03 - 8e-4r for r in rs]); F2σH = interp1_order8(rs, [0.001 + 2e-5r for r in rs])
        f = Forcing1PAT1ea(base, F2σI, F2σH, FaI, FaH)
        νg = [0.05, 0.1, 0.15]; c1g = [-0.5, 0.0, 0.5]; c2g = [-0.5, 0.0, 0.5]
        Ωc(ν, χ1, χ2) = 0.05 + 0.1ν + 0.03χ1 + 0.02χ2 + 0.05ν * χ1
        grid = [((ν, χ1, χ2), Ωc(ν, χ1, χ2)) for ν in νg for χ1 in c1g for χ2 in c2g]
        model = Model1PAT1ea(f, interp3_lagrange(grid))

        # Weak spins keep the additive forcing positive (monotonic chirp) up to Ωstop: the
        # primary-spin back-reaction channel ∝ χt1/(r0-6)² would otherwise flip dΩ/dt<0 and the
        # synthetic inspiral would never reach the stop frequency.
        ics = (; Ω = 0.05, φ = 0.0, ν = 0.1, M = 1.0, χt1 = 0.02, χt2 = -0.02, δm = 0.4)
        s = int_inspiral(model, ics)
        Ωstop = omega_crit(model, ics.ν, ics.χt1, ics.χt2)
        @test omega_at(s, duration(s)) ≈ Ωstop rtol = 1e-7
        td = duration(s)
        for t in (0.3 * td, 0.7 * td)
            @test s.sol(t)[3] == ics.ν           # ν frozen
            @test s.sol(t)[4] == ics.M           # M frozen
            @test s.sol(t)[6] == ics.χt2         # χt2 frozen
            r0 = r0_at(s, t)
            h = 1.0e-5 * td
            dχ1 = (chit1_at(s, t + h) - chit1_at(s, t - h)) / (2h)
            @test dχ1 ≈ (ics.ν^2 / ics.M) * flux_LH(f, r0) rtol = 1e-4
            ddm = (dm_at(s, t + h) - dm_at(s, t - h)) / (2h)
            @test ddm ≈ (ics.ν / ics.M) * FEH(r0) rtol = 1e-4
        end
        @test chit1_at(s, td) > ics.χt1 && dm_at(s, td) > ics.δm   # both grow

        batch_ics = [ics, merge(ics, (; δm = 0.2)), merge(ics, (; δm = 0.6))]
        ser = int_inspiral_batch(model, batch_ics; threaded = false)
        thr = int_inspiral_batch(model, batch_ics; threaded = true)
        for i in eachindex(batch_ics)
            @test duration(ser[i]) == duration(thr[i])
            @test omega_at(ser[i], 0.5duration(ser[i])) == omega_at(thr[i], 0.5duration(thr[i]))
        end
    end
end
