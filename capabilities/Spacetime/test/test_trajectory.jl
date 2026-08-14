# Geodesics/trajectory.jl — BAKED characterization vs the KerrGeodesics oracle (values captured
# once from KerrGeoTrajectory{Equatorial,Radial,Polar}/KerrGeoTrajectory/KerrGeoProperTime, Phase7-B)
# over representative + edge (a,p,e,x) and several radial/polar Mino phases. The tolerance reflects
# the NumericsCore-vs-KerrGeodesics elliptic-layer agreement already established by the frequencies
# test (~1e-9). Plus the closed-form Schwarzschild circular check (e=0 ⇒ r≡p, deviations ≡ 0).
using Test
using Spacetime
include(joinpath(@__DIR__, "fixtures", "spacetime_trajectory_ref.jl"))  # TRAJ_{EQ,RAD,POL,FULL,PT}_REF

@testset "Geodesics: closed-form trajectories" begin
    qrs = (0.0, 0.7, 1.6, 3.0, 4.5, 6.0)
    qzs = (0.0, 0.9, 2.1, 3.7, 5.5)

    @testset "geodesic_trajectory_equatorial == KerrGeoTrajectoryEquatorial (baked): (a,p,e)=($a,$p,$e)" for (ci, (a, p, e)) in
        enumerate(((0.5, 10.0, 0.2), (0.9, 8.0, 0.4), (0.0, 12.0, 0.3), (-0.7, 11.0, 0.25)))
        f  = geodesic_trajectory_equatorial(a, p, e)
        for (qi, qr) in enumerate(qrs)
            g = f(qr); go = TRAJ_EQ_REF[ci][qi]
            @test isapprox(g.r,   go.r;   rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.Ur,  go.Ur;  rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.Δtr, go.Δtr; rtol = 1e-9, atol = 1e-10)
            @test isapprox(g.Δφr, go.Δφr; rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.Δτr, go.Δτr; rtol = 1e-9, atol = 1e-10)
        end
    end

    @testset "geodesic_trajectory_radial == KerrGeoTrajectoryRadial (baked): (a,p,e,x)=$c" for (ci, c) in
        enumerate(((0.5, 10.0, 0.2, 0.7), (0.9, 12.0, 0.1, 0.5), (0.0, 10.0, 0.3, 1.0), (0.998, 7.0, 0.15, 0.8)))
        a, p, e, x = c
        f  = geodesic_trajectory_radial(a, p, e, x)
        for (qi, qr) in enumerate(qrs)
            g = f(qr); go = TRAJ_RAD_REF[ci][qi]
            @test isapprox(g.r,   go.r;   rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.Ur,  go.Ur;  rtol = 1e-9, atol = 1e-10)
            @test isapprox(g.Δtr, go.Δtr; rtol = 1e-9, atol = 1e-10)
            @test isapprox(g.Δφr, go.Δφr; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "geodesic_trajectory_polar == KerrGeoTrajectoryPolar (baked): (a,p,e,x)=$c" for (ci, c) in
        enumerate(((0.5, 10.0, 0.2, 0.7), (0.9, 12.0, 0.1, 0.5), (0.998, 7.0, 0.15, 0.8)))
        a, p, e, x = c
        f  = geodesic_trajectory_polar(a, p, e, x)
        for (qi, qz) in enumerate(qzs)
            g = f(qz); go = TRAJ_POL_REF[ci][qi]
            @test isapprox(g.z,   go.z;   rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.θ,   go.θ;   rtol = 1e-9, atol = 1e-11)
            @test isapprox(g.Uz,  go.Uz;  rtol = 1e-9, atol = 1e-10)
            @test isapprox(g.Δtz, go.Δtz; rtol = 1e-9, atol = 1e-10)
            @test isapprox(g.Δφz, go.Δφz; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "geodesic_trajectory (full t/r/θ/φ) == KerrGeoTrajectory (baked): (a,p,e,x)=$c" for (ci, c) in
        enumerate(((0.5, 10.0, 0.2, 0.7), (0.9, 12.0, 0.1, 0.5)))
        a, p, e, x = c
        f  = geodesic_trajectory(a, p, e, x)
        ref = TRAJ_FULL_REF[ci]                            # t/φ flattened row-major: (qr outer, qz inner)
        for (qi, qr) in enumerate(qrs)
            @test isapprox(f.r(qr), ref.r[qi]; rtol = 1e-9, atol = 1e-11)
        end
        for (zi, qz) in enumerate(qzs)
            @test isapprox(f.θ(qz), ref.θ[zi]; rtol = 1e-9, atol = 1e-11)
        end
        for (qi, qr) in enumerate(qrs), (zi, qz) in enumerate(qzs)
            k = (qi - 1) * length(qzs) + zi
            @test isapprox(f.t(0.3, qr, qz), ref.t[k]; rtol = 1e-9, atol = 1e-10)
            @test isapprox(f.φ(0.2, qr, qz), ref.φ[k]; rtol = 1e-9, atol = 1e-10)
        end
    end

    @testset "geodesic_proper_time == KerrGeoProperTime (inclined, baked): (a,p,e,x)=$c" for (ci, c) in
        enumerate(((0.5, 10.0, 0.2, 0.7), (0.9, 12.0, 0.1, 0.5), (0.998, 7.0, 0.15, 0.8)))
        a, p, e, x = c
        pt  = geodesic_proper_time(a, p, e, x)
        pto = TRAJ_PT_REF[ci]                              # ProperTime flattened row-major: (qr outer, qz inner)
        @test isapprox(pt.ProperTimeFrequency, pto.ProperTimeFrequency; rtol = 1e-9, atol = 1e-11)
        @test isapprox(pt.Υτz, pto.Υτz; rtol = 1e-9, atol = 1e-11)
        for (qi, qr) in enumerate(qrs)
            @test isapprox(pt.Δτr(qr), pto.Δτr[qi]; rtol = 1e-9, atol = 1e-10)
        end
        for (zi, qz) in enumerate(qzs)
            @test isapprox(pt.Δτz(qz), pto.Δτz[zi]; rtol = 1e-9, atol = 1e-10)
        end
        for (qi, qr) in enumerate(qrs), (zi, qz) in enumerate(qzs)
            k = (qi - 1) * length(qzs) + zi
            @test isapprox(pt.ProperTime(0.4, qr, qz), pto.ProperTime[k]; rtol = 1e-9, atol = 1e-10)
        end
    end

    @testset "analytic circular (e=0): r≡p, Ur≡0 (all qr); deviations ≡ 0 on the principal half-period" begin
        # The equatorial convention parametrises r via asin(sn), valid on qr ∈ [0,π]; the oscillating
        # deviations of t,φ,τ then vanish identically for a circular orbit there (physical truth).
        # (Beyond π the asin branch reflects — a known convention limitation the live e>0 tests cover.)
        for (a, p) in ((0.0, 10.0), (0.7, 8.0))
            f = geodesic_trajectory_equatorial(a, p, 0.0)
            for qr in qrs
                @test isapprox(f(qr).r, p; rtol = 1e-12)   # circular radius is constant for all qr
                @test isapprox(f(qr).Ur, 0.0; atol = 1e-12)
            end
            for qr in (0.0, 0.7, 1.6, 3.0)                 # qr ≤ π: deviations identically zero
                g = f(qr)
                @test isapprox(g.Δtr, 0.0; atol = 1e-10)
                @test isapprox(g.Δφr, 0.0; atol = 1e-10)
                @test isapprox(g.Δτr, 0.0; atol = 1e-10)
            end
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, p, e = big"0.5", big"10.0", big"0.2"
        g = geodesic_trajectory_equatorial(a, p, e)(big"1.6")
        @test g.r isa BigFloat
        @test g.Δtr isa BigFloat
    end
end
