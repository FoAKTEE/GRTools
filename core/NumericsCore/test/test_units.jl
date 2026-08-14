# Characterization of split-named physical constants + convention helpers.
#
# Heavy KerrBackwards/FastEMRIWaveforms test-deps DROPPED (see Project.toml); the constants are
# re-grounded WITHOUT a live dep:
# Oracle 1 (KB canonical literals): exact decimal literals READ from
#   packages/KerrBackwards/src/constants.jl (the LAL/global.h #defines) — pasted here as expected.
# Oracle 2 (FEW re-derivation): the lisaconstants 2.0.2 derivation EXPRESSIONS read from
#   packages/FastEMRIWaveforms/src/constants.jl. FEW's GM_SUN differs slightly, so MTSUN/MRSUN agree
#   with the canonical LAL values only to ~1e-7 (asserted approximately); YRSID/Gpc agree exactly.
# Oracle 3 (cross-identities): MRSUN_SI = MTSUN_SI*c, Gpc_in_s = Gpc_in_m/c (c=299792458).
# Oracle 4 (p_to_y): closed form log(p - 2e - 21/10) — this IS the KB source formula read from
#   packages/KerrBackwards/src/utility.jl:15 (`p_to_y(p,e) = log(p - 2*e - 21//10)`).
# Oracle 5 (eps<->omega / q<->a / imagunit): inverse round-trip + analytic convention definitions
#   (eps=2*M*omega, a=q*M, imagunit^2=-1).
# Every expected value is a source literal or an analytic relation — never the code-under-test output.

using NumericsCore

@testset "units" begin
    c = 299792458.0

    @testset "constants exact vs KerrBackwards source literals (canonical LAL)" begin
        # Literals copied verbatim from packages/KerrBackwards/src/constants.jl (LAL/global.h #defines).
        @test NumericsCore.MTSUN_SI == 4.925491025873693e-06
        @test NumericsCore.YRSID_SI == 31558149.763545603
        @test NumericsCore.MRSUN_SI == 1476.6250615036158
        @test NumericsCore.AUsec    == 499.004783836156412
        @test NumericsCore.Gpc_in_m == 3.0856775814913674e25     # KB `Gpc` (meters axis)
        @test NumericsCore.Gpc_in_s == 1.029271251e17            # KB `GPCINSEC` (seconds axis)
        @test NumericsCore.GPCINSEC == 1.029271251e17            # legacy alias retained
        @test NumericsCore.PI   == 3.141592653589793238462643383279502884
        @test NumericsCore.PI_2 == 1.570796326794896619231321691639751442
        @test NumericsCore.PI_4 == 0.785398163397448309615660845819875721
        @test NumericsCore.PI == Float64(pi)
    end

    @testset "constants approx vs FastEMRIWaveforms source re-derivation" begin
        # Derivation expressions read verbatim from packages/FastEMRIWaveforms/src/constants.jl
        # (lisaconstants 2.0.2 base values: GM_SUN=1.3271244e20, PARSEC=3.0856775814913674e16, ...).
        few_MTSUN = 1.3271244e20 / (299792458.0^3)
        few_MRSUN = 0.5 * (2 * 1.3271244e20 / (299792458.0^2))
        few_YRSID = 365.256363004 * 60 * 60 * 24
        few_Gpc   = 1e9 * 3.0856775814913674e16
        @test NumericsCore.MTSUN_SI ≈ few_MTSUN rtol=1e-7   # FEW GM_SUN differs from LAL at ~1.6e-8
        @test NumericsCore.MRSUN_SI ≈ few_MRSUN rtol=1e-7
        @test NumericsCore.YRSID_SI ≈ few_YRSID rtol=1e-7
        @test NumericsCore.Gpc_in_m == few_Gpc              # meters axis agrees exactly (PARSEC*1e9)
    end

    @testset "analytic cross-identities" begin
        @test NumericsCore.MRSUN_SI ≈ NumericsCore.MTSUN_SI*c rtol=1e-15   # GM/c^2 = (GM/c^3)*c
        @test NumericsCore.Gpc_in_s ≈ NumericsCore.Gpc_in_m/c rtol=1e-9    # light-time = length/c
        @test NumericsCore.PI_2 ≈ NumericsCore.PI/2 rtol=1e-15
        @test NumericsCore.PI_4 ≈ NumericsCore.PI/4 rtol=1e-15
    end
end

@testset "conventions" begin
    @testset "p_to_y closed form (== KerrBackwards utility.jl:15)" begin
        # KB source: `p_to_y(p, e) = log(p - 2*e - 21//10)` (packages/KerrBackwards/src/utility.jl:15).
        for (p, e) in ((10.0, 0.1), (6.5, 0.3), (8.0, 0.0))
            @test NumericsCore.p_to_y(p, e) ≈ log(p - 2*e - 21//10) rtol=1e-15
        end
    end

    @testset "eps<->omega (eps=2*M*omega) round-trip" begin
        for eps in (0.05, 0.3, 1.2)
            @test NumericsCore.eps_to_omega(eps) ≈ eps/2 rtol=1e-15        # M=1
            @test NumericsCore.omega_to_eps(NumericsCore.eps_to_omega(eps)) ≈ eps rtol=1e-15
        end
        @test NumericsCore.omega_to_eps(0.4, 2.0) ≈ 2*2.0*0.4 rtol=1e-15   # general M
        @test NumericsCore.eps_to_omega(NumericsCore.omega_to_eps(0.7, 3.0), 3.0) ≈ 0.7 rtol=1e-15
    end

    @testset "q<->a (a=q*M) round-trip" begin
        for q in (0.0, 0.5, 0.998)
            @test NumericsCore.q_to_a(q) == q                              # M=1 identity
            @test NumericsCore.a_to_q(NumericsCore.q_to_a(q)) ≈ q rtol=1e-15
        end
        @test NumericsCore.q_to_a(0.5, 4.0) ≈ 0.5*4.0 rtol=1e-15
        @test NumericsCore.a_to_q(NumericsCore.q_to_a(0.9, 5.0), 5.0) ≈ 0.9 rtol=1e-15
    end

    @testset "imagunit" begin
        @test NumericsCore.imagunit() == im
        @test NumericsCore.imagunit()^2 == -1
        @test NumericsCore.imagunit(BigFloat) isa Complex{BigFloat}
        @test NumericsCore.imagunit(BigFloat)^2 == -1
    end
end
