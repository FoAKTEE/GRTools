# Grounding for hertz_potential_mode (RadHzpotmode, Sourceincoordinates.nb).
#
# ORACLE: frozen `wolframscript` values from the upstream `RadHzpotmode` (BHPToolkit `Teukolsky`
# paclet) at fixed (l,m,n,k,r) on circular and spherical Schwarzschild geodesics (NOT legacy-code
# output — the paclet captured InputForm values):
#   circ a=0 p=10 l=2 m=2: ω=0.06324555320336758,
#     RadHZ(7)=-0.014036995034207932 - 0.008501341967457264i,
#     dRadHZ(7)=0.012436577381918407 + 0.006341609877163205i
#   sph  a=0 p=12 x=0.5 l=2 m=2 k=0: ω=0.04811252243246881,
#     RadHZ(8)=-0.0056928017370048815 - 0.002819766511723147i
# Agreement is MST-truncation limited (Julia Nmax=80 vs paclet), ~1e-7 relative.
#
# ORACLE-INDEPENDENT regression: the minus-mode self-symmetry used by `Sourcesamplitudes`,
#   RadHZ(l,-m) = (-1)^(l+kloc) conj(RadHZ(l,+m)),  ω(l,-m) = -ω(l,+m),
# holds to machine precision (Teukolsky m→-m conjugation), independent of any oracle.

@testset "hertz_potential_mode" begin
    @testset "circular oracle (a=0,p=10,l=2,m=2)" begin
        RadHZ, dRadHZ, ω = hertz_potential_mode(2, 2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        @test ω ≈ 0.06324555320336758 rtol=1e-10
        @test RadHZ(7.0) ≈ (-0.014036995034207932 - 0.008501341967457264im) rtol=1e-7
        @test dRadHZ(7.0) ≈ (0.012436577381918407 + 0.006341609877163205im) rtol=1e-7
    end

    @testset "spherical oracle (a=0,p=12,x=0.5,l=2,m=2,k=0)" begin
        RadHZ, _, ω = hertz_potential_mode(2, 2, 0, 0, 0.0, 12.0, 0.0, 0.5)
        @test ω ≈ 0.04811252243246881 rtol=1e-10
        @test RadHZ(8.0) ≈ (-0.0056928017370048815 - 0.002819766511723147im) rtol=1e-7
    end

    @testset "minus-mode self-symmetry (oracle-independent)" begin
        l, kloc, r0 = 2, 0, 7.0
        RadHZp, _, ωp = hertz_potential_mode(l,  2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        RadHZm, _, ωm = hertz_potential_mode(l, -2, 0, 0, 0.0, 10.0, 0.0, 1.0)
        @test ωm ≈ -ωp atol=1e-14
        @test RadHZm(r0) ≈ (-1)^(l + kloc) * conj(RadHZp(r0)) rtol=1e-10
    end

    @testset "eccentric orbit is not ported (honest throw)" begin
        @test_throws ArgumentError hertz_potential_mode(2, 2, 0, 0, 0.0, 10.0, 0.1, 1.0)
    end

    @testset "static mode is rejected instead of returning NaN" begin
        @test_throws DomainError hertz_potential_mode(2, 0, 0, 0, 0.0, 10.0, 0.0, 1.0)
    end
end
