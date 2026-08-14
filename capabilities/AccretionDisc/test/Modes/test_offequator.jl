# Modes.OffEquator — off-equatorial (z=cosθ) disc background + RT-A2 Q^{-1}_z/Q_z + Gaussian
# height-average. GROUNDING (Phase-7B: no live packages/* dep):
#   • closed-form identities: z=0 reductions to the equatorial Background/Velocity pieces, the
#     metric-normalization identity α = 1/√(−(g_tt+2Ωg_tφ+Ω²g_φφ)) at z=0, the exact Gaussian
#     vertical profile ratio e^{−z²/2h²}, U·U=−1 at z=0, and the Gaussian moments ⟨z²⟩=h²,
#     ⟨z⁴⟩=3h⁴ of the height-average.
#   • Q_z·Q_z^{-1} = I to <1e-10 off the equator — verbatim the legacy DiscModes q_resid gate.
#   • CAPTURED one-shot legacy DiscModes oracle literals (packages/DiscModes OmegaK_z/alpha_z/
#     Uvec_z/pressure_z/eden_z/build_Qinv_z/build_Q_z/heightavg at a=0.6,h0=0.1,r=12,z=0.05,
#     r_p=15,m=2 and a=0, captured 2026-07-02 by scratchpad/capture_s55_accretiondisc_force.jl).
#     Tolerances: exact-formula quantities rtol 1e-12 (identical algebra); Q=inv(Q^{-1}) entries
#     rtol 1e-6 — the ONLY implementation delta is the analytic Spacetime.christoffel replacing the
#     legacy FD kerr_christoffel fork (the Velocity.jl dedup decision); measured discrepancy ≤2.5e-8,
#     1e-6 keeps a ~40× margin while pinning 6 digits.
using Test
using AccretionDisc
using AccretionDisc.Modes: RelativisticDisc, disc_OmegaK, csound, eden0,
                           disc_OmegaK_z, alpha_z, Uvec_z, pressure_z, eden_z,
                           qinv_operator, qinv_operator_z, q_operator_z, heightavg
import Spacetime
using LinearAlgebra

@testset "Modes.OffEquator (z=cosθ background + Q_z + height-average)" begin
    d = RelativisticDisc(a = 0.6, h0 = 0.1, r_I = 1.0, eden_I = 1.0)
    r = 12.0; z = 0.05; rp = 15.0; m = 2

    @testset "z=0 reductions + closed-form identities" begin
        @test disc_OmegaK_z(d, r, 0.0) == disc_OmegaK(d, r)          # Ω^k(r,0) = Ω^k_eq exactly
        @test pressure_z(d, r, 0.0) == csound(d, r)^2 * eden0(d, r)  # p(r,0) = c_s²ē
        @test eden_z(d, r, 0.0) == eden0(d, r)
        # exact Gaussian vertical profile: p(r,z)/p(r,0) = e^{−z²/2h₀²}
        @test pressure_z(d, r, z)/pressure_z(d, r, 0.0) ≈ exp(-z^2/(2*0.1^2)) rtol = 1e-14
        @test eden_z(d, r, z)/eden0(d, r) ≈ exp(-z^2/(2*0.1^2)) rtol = 1e-14
        # α at z=0 == the metric-normalization of the circular 4-velocity (Spacetime.kerr_metric)
        Ω = disc_OmegaK(d, r); g = Spacetime.kerr_metric(0.6, r, π/2)
        αmet = 1/sqrt(-(g[1,1] + 2Ω*g[1,4] + Ω^2*g[4,4]))
        @test alpha_z(d, r, 0.0) ≈ αmet rtol = 1e-12
        # U·U = −1 at z=0 (exact normalization of the equatorial circular 4-velocity)
        U = Uvec_z(d, r, 0.0)
        @test abs(sum(g[μ,ν]*U[μ]*U[ν] for μ in 1:4, ν in 1:4) + 1) < 1e-12
    end

    @testset "captured legacy DiscModes literals (a=0.6, r=12, z=0.05, r_p=15, m=2)" begin
        @test disc_OmegaK_z(d, r, z) ≈ 0.023743824489938635 rtol = 1e-13
        @test alpha_z(d, r, z) ≈ 1.1496047289588343 rtol = 1e-13
        U = Uvec_z(d, r, z)
        @test U[1] ≈ 1.1496047289588343 rtol = 1e-13
        @test U[2] == 0.0 && U[3] == 0.0
        @test U[4] ≈ 0.027296012917202035 rtol = 1e-13
        @test pressure_z(d, r, z) ≈ 6.2339534814522242e-6 rtol = 1e-13
        @test eden_z(d, r, z) ≈ 0.0061284507123930236 rtol = 1e-13
        Qi = qinv_operator_z(d, r, z, rp, m)
        @test Qi[1,1] ≈ (0.0 - 0.015419707179787552im) rtol = 1e-12
        @test Qi[2,2] ≈ (0.0 - 0.016872317229274791im) rtol = 1e-12
        @test Qi[4,2] ≈ (0.00089260769075949085 + 0.0im) rtol = 1e-6   # Γ-carrying entry (christoffel swap)
        @test Qi[2,1] ≈ (0.013155549281713854 + 0.0im) rtol = 1e-6
        Q = q_operator_z(d, r, z, rp, m)
        @test Q[1,1] ≈ (0.0 + 78.08866241808488im) rtol = 1e-6
        @test Q[2,2] ≈ (0.0 - 73.638557081059844im) rtol = 1e-6
        @test Q[1,4] ≈ (0.0 - 565.5273515730986im) rtol = 1e-6
        @test Q[4,2] ≈ (-4.4155603468010778 + 0.0im) rtol = 1e-6
        @test Q[3,3] ≈ (0.0 + 60.201107905973032im) rtol = 1e-6
        @test Q[2,1] ≈ (-62.848522243232452 + 0.0im) rtol = 1e-6
        # a=0 probe in the Force-pipeline region (r=11, z=0.02, r_p=8, m=2)
        d0 = RelativisticDisc(a = 0.0, h0 = 0.1, r_I = 1.0, eden_I = 1.0)
        Q0 = q_operator_z(d0, 11.0, 0.02, 8.0, 2)
        @test Q0[1,1] ≈ (0.0 - 24.241214760727651im) rtol = 1e-6
        @test Q0[2,2] ≈ (0.0 - 35.222532757073516im) rtol = 1e-6
    end

    @testset "Q_z·Q_z^{-1} = I off the equator (legacy q_resid gate) + z→0 reduction" begin
        mx = 0.0
        for (rr, zz) in [(12.0, 0.05), (20.0, 0.08), (8.0, 0.03)]
            Qi = qinv_operator_z(d, rr, zz, rp, m)
            Q  = q_operator_z(d, rr, zz, rp, m)
            @test all(isfinite, Qi)
            mx = max(mx, maximum(abs.(Q * Qi - I)))
        end
        @test mx < 1e-10                                   # the legacy q_resid<1e-10 identity
        # z=0 reduces to the equatorial qinv_operator (measured 3.5e-18; same FD steps, same Γ)
        @test maximum(abs.(qinv_operator_z(d, r, 0.0, rp, m) .- qinv_operator(d, r, rp, m))) < 1e-12
    end

    @testset "heightavg — exact Gaussian moments (+ captured legacy samples)" begin
        @test heightavg(zz -> zz^2, 0.1) ≈ 0.01 rtol = 1e-10           # ⟨z²⟩ = h²
        @test heightavg(zz -> one(zz), 0.1) == 1.0                      # ⟨1⟩ = 1
        @test heightavg(zz -> zz^4, 0.1; order = 4) ≈ 3e-4 rtol = 1e-4  # ⟨z⁴⟩ = 3h⁴
        # captured legacy DiscModes.heightavg samples (identical formula ⇒ tight)
        @test heightavg(zz -> cos(zz), 0.1) ≈ 0.995 rtol = 1e-10
        @test heightavg(zz -> cos(zz), 0.1; hz = 0.02) ≈ 0.99500016666444391 rtol = 1e-12
    end

    @testset "BigFloat inheritance" begin
        db = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
        @test disc_OmegaK_z(db, big"12.0", big"0.05") isa BigFloat
        @test Uvec_z(db, big"12.0", big"0.05") isa Vector{BigFloat}
        Qib = qinv_operator_z(db, big"12.0", big"0.05", big"15.0", 2)
        Qb  = q_operator_z(db, big"12.0", big"0.05", big"15.0", 2)
        @test eltype(Qib) == Complex{BigFloat}
        @test maximum(abs.(Qb * Qib - I)) < 1e-30           # BigFloat identity far tighter than Float64
        @test heightavg(zz -> zz^2, big"0.1") isa BigFloat
    end
end
