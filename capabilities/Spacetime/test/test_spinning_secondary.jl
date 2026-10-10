# SpinningSecondary — BAKED characterization vs the SpinningOrbit oracle (values captured once from
# KerrSpinLinearParts/KerrSpinOrbit/KerrSpinOrbitEquatorial, Phase7-B) at its frozen tolerances:
# CTilde linear parts (rtol 1e-9/1e-8 per runtests L47-52) and KerrSpinOrbit equatorial
# scalars+freqs+trajectory (rtol 1e-9/1e-8 per runtests). Plus the spin→0 reduction identity:
# spar=0 ⇒ spin_orbit constants/frequencies equal the in-package Geodesics.
using Test
using Spacetime
using .EquatorialQuadratureOracle: quadrature_orbit
using Spacetime.SpinningSecondary: spin_linear_constants, spin_linear_frequencies,
                                   spin_mino_time_correction, spin_orbit, spin_orbit_equatorial
include(joinpath(@__DIR__, "fixtures", "spacetime_spinning_ref.jl"))  # SPIN_{LIN,ORBIT,ORBIT_EQ}_REF

@testset "SpinningSecondary: linear-in-spin secondary" begin
    # representative + edge: equatorial (x=1), inclined, near-extremal, retrograde
    lincases = [
        (0.5,  10.0, 0.2, 1.0),
        (0.9,  10.0, 0.0, 1.0),
        (0.9,  12.0, 0.1, 0.5),
        (0.998, 7.0, 0.15, 0.8),
        (-0.7, 11.0, 0.2, 0.6),
    ]

    @testset "spin_linear_constants/frequencies == KerrSpinLinearParts CTilde (baked): (a,p,e,x)=$c" for (ci, c) in enumerate(lincases)
        a, p, e, x = c
        lp = SPIN_LIN_REF[ci]
        dc = spin_linear_constants(a, p, e, x)
        for i in 1:3
            @test isapprox(dc[i], lp.C1CTilde[i]; rtol = 1e-9, atol = 1e-12)
        end
        df = spin_linear_frequencies(a, p, e, x)
        for i in 1:3
            @test isapprox(df[i], lp.Omega1CTilde[i]; rtol = 1e-8, atol = 1e-12)
        end
        @test isapprox(spin_mino_time_correction(a, p, e, x), lp.Upsilont1CTilde; rtol = 1e-9, atol = 1e-12)
    end

    @testset "spin_orbit (equatorial): baked scalars + independent signed-spin trajectory: (a,p,e,spar,gauge)=$c" for (ci, c) in enumerate([
            (0.5, 10.0, 0.2, 0.1, "pex"), (0.9, 8.0, 0.3, 0.05, "pex"),
            (0.5, 10.0, 0.2, 0.1, "CTilde"), (-0.6, 11.0, 0.15, 0.08, "pex")])
        a, p, e, spar, gauge = c
        o  = spin_orbit(a, p, e, 1.0, spar; gauge = gauge)
        oo = SPIN_ORBIT_REF[ci]
        independent = nothing
        if a < 0
            # Replace the wrong-chart scalar fixture using the radial-root spin
            # displacement and a separately integrated shifted geodesic.
            ab,pb,eb,sb=BigFloat.((a,p,e,spar))
            base=quadrature_orbit(ab,pb,eb)
            K=(base.L-ab*base.E)^2
            ra=pb/(1-eb); rp=pb/(1+eb)
            dra=-((ra^2+ab^2)*base.E-ab*base.L)/(ra*sqrt(K))
            drp=-((rp^2+ab^2)*base.E-ab*base.L)/(rp*sqrt(K))
            pT=pb-sb*2(rp^2*dra+ra^2*drp)/(ra+rp)^2
            eT=eb-sb*2(rp*dra-ra*drp)/(ra+rp)^2
            independent=quadrature_orbit(ab,pT,eT)
            Et,Lt=independent.E,independent.L
            Kt=(Lt-ab*Et)^2
            oo=(EnTilde=Et,JzTilde=Lt,KTilde=Kt,pTilde=pT,eTilde=eT,
                Energy=Et-sb*(1-Et^2)/(2sqrt(Kt)),
                AngularMomentum=Lt-sb*(ab-Et*Lt/2)/sqrt(Kt),
                CarterConstantK=Kt-sb*(3ab*(Lt-ab*Et)-Et*Kt)/sqrt(Kt),
                Frequencies=(Ur=independent.Υr,Uph=independent.Υφ,
                    Ut=independent.Γ-3sb*independent.Υτ/(2sqrt(Kt))),
                ProperTimeFrequency=independent.Υτ)
        end
        for f in (:EnTilde, :JzTilde, :KTilde, :pTilde, :eTilde, :Energy, :AngularMomentum, :CarterConstantK)
            @test isapprox(getfield(o, f), getfield(oo, f); rtol = 1e-9, atol = 1e-12)
        end
        @test isapprox(o.Frequencies.Ur,  oo.Frequencies.Ur;  rtol = 1e-8, atol = 1e-11)
        @test isapprox(o.Frequencies.Uph, oo.Frequencies.Uph; rtol = 1e-8, atol = 1e-11)
        @test isapprox(o.Frequencies.Ut,  oo.Frequencies.Ut;  rtol = 1e-8, atol = 1e-11)
        @test isapprox(o.ProperTimeFrequency, oo.ProperTimeFrequency; rtol = 1e-9, atol = 1e-11)
        # Independently integrate the corrected shifted geodesic and apply
        # Δt_spin = Δt_geo - 3s/(2√K̃) Δτ_geo.
        for (qi, qr) in enumerate((0.0, 0.7, 1.6, 3.0))
            reference = if independent === nothing
                oo.trajectory[qi]
            else
                g=independent.at_phase(BigFloat(qr))
                (r=g.r, Δtr=g.Δtr-3BigFloat(spar)/(2sqrt(BigFloat(oo.KTilde)))*g.Δτr, Δφr=g.Δφr)
            end
            @test isapprox(o.trajectory(qr).r, reference.r; rtol = 1e-9, atol = 1e-11)
            @test isapprox(o.trajectory(qr).Δtr, reference.Δtr; rtol = 1e-8, atol = 1e-10)
            @test isapprox(o.trajectory(qr).Δφr, reference.Δφr; rtol = 1e-8, atol = 1e-11)
        end
    end

    @testset "spin_orbit_equatorial == KerrSpinOrbitEquatorial (baked): (a,p,e)=$c" for (ci, c) in
        enumerate(((0.5, 10.0, 0.2), (0.9, 8.0, 0.3), (0.0, 12.0, 0.25)))
        a, p, e = c
        o  = spin_orbit_equatorial(a, p, e)
        oo = SPIN_ORBIT_EQ_REF[ci]
        @test o.Inclination == 1
        for f in (:Energy, :AngularMomentum, :CarterConstant, :ProperTimeFrequency)
            @test isapprox(getfield(o, f), getfield(oo, f); rtol = 1e-9, atol = 1e-12)
        end
        for f in (:Ur, :Uth, :Uph, :Ut)
            @test isapprox(getfield(o.Frequencies, f), getfield(oo.Frequencies, f); rtol = 1e-9, atol = 1e-11)
        end
        for (qi, qr) in enumerate((0.0, 0.7, 1.6, 3.0))
            @test isapprox(o.trajectory(qr).r, oo.trajectory_r[qi]; rtol = 1e-9, atol = 1e-11)
        end
    end

    @testset "spin→0 reduction: spin_orbit constants/frequencies == in-package Geodesics" begin
        for (a, p, e) in ((0.0, 10.0, 0.2), (0.9, 10.0, 0.1), (0.5, 8.0, 0.0))
            o = spin_orbit(a, p, e, 1.0, 0.0)             # spar = 0 ⇒ pure geodesic
            En, Lz, Q = constants_of_motion(a, p, e, 1.0)
            @test isapprox(o.Energy, En; rtol = 1e-13)
            @test o.EnTilde == En                          # tilde == geodesic at spar=0
            @test isapprox(o.AngularMomentum, Lz; rtol = 1e-13)
            @test isapprox(o.CarterConstantK, Q + (Lz - a * En)^2; rtol = 1e-13)
            fm = kerr_frequencies(a, p, e, 1.0; Time = "Mino")
            @test isapprox(o.Frequencies.Ur,  fm.Ur;  rtol = 1e-13)
            @test isapprox(o.Frequencies.Uph, fm.Uph; rtol = 1e-13)
            @test isapprox(o.Frequencies.Ut,  fm.Ut;  rtol = 1e-13)   # δΥt · 0 vanishes
        end
    end

    @testset "out-of-slice + bad gauge throw" begin
        @test_throws ArgumentError spin_orbit(0.9, 10.0, 0.1, 0.5, 0.1)        # inclined x²≠1
        @test_throws ArgumentError spin_orbit(0.5, 10.0, 0.2, 1.0, 0.1; gauge = "bogus")
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, p, e = big"0.5", big"10.0", big"0.2"
        @test spin_linear_constants(a, p, e, big"1.0")[1] isa BigFloat
        @test spin_orbit(a, p, e, big"1.0", big"0.1").Energy isa BigFloat
    end
end
