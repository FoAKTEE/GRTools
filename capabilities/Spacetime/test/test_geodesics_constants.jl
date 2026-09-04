# Geodesics/constants.jl — grounded by BAKED characterization vs the KerrGeodesics oracle
# (values captured once from KerrGeoEnergy/AngularMomentum/CarterConstant, Phase7-B) at the legacy
# runtests tolerance, plus the closed-form Schwarzschild circular-equatorial limit. Spacetime's
# geodesic_energy/angular_momentum/carter_constant must equal the frozen oracle references.
using Test
using Spacetime
include(joinpath(@__DIR__, "fixtures", "spacetime_constants_ref.jl"))  # const CONSTANTS_REF

@testset "Geodesics: constants of motion" begin
    # Representative (a,p,e,x): a=0 Schwarzschild, circular e=0, equatorial x=1, generic inclined,
    # near-extremal, retrograde a<0. x = cos(inclination).
    cases = [
        (0.0,   10.0, 0.0, 1.0),  # Schwarzschild circular equatorial
        (0.0,    8.0, 0.3, 1.0),  # Schwarzschild eccentric equatorial
        (0.9,   10.0, 0.0, 1.0),  # Kerr circular equatorial prograde
        (0.9,   10.0, 0.1, 0.5),  # generic inclined eccentric
        (0.998,  6.0, 0.2, 0.8),  # near-extremal generic
        (-0.7,  12.0, 0.4, 0.6),  # retrograde generic
        (0.5,    7.0, 0.0, 0.3),  # circular inclined
    ]
    @testset "baked vs KerrGeodesics oracle (captured once): (a,p,e,x)=$c" for (i, c) in enumerate(cases)
        a, p, e, x = c
        ref = CONSTANTS_REF[i]
        @test isapprox(geodesic_energy(a, p, e, x),           ref.E;   rtol = 1e-12)
        @test isapprox(geodesic_angular_momentum(a, p, e, x), ref.Lz;  rtol = 1e-12, atol = 1e-13)
        @test isapprox(carter_constant(a, p, e, x),           ref.Q;   rtol = 1e-12, atol = 1e-13)
        com = constants_of_motion(a, p, e, x)
        @test com.E  == geodesic_energy(a, p, e, x)
        @test com.Lz == geodesic_angular_momentum(a, p, e, x)
        @test com.Q  == carter_constant(a, p, e, x)
    end

    @testset "analytic Schwarzschild circular-equatorial limit (a=0,e=0,x=1)" begin
        # E = (1-2M/r)/sqrt(1-3M/r) = (p-2)/sqrt(p(p-3)); Lz = sqrt(Mr)/sqrt(1-3M/r) = p/sqrt(p-3); Q = 0.
        for p in (6.5, 8.0, 10.0, 20.0)
            @test isapprox(geodesic_energy(0.0, p, 0.0, 1.0),           (p - 2) / sqrt(p * (p - 3)); rtol = 1e-13)
            @test isapprox(geodesic_angular_momentum(0.0, p, 0.0, 1.0), p / sqrt(p - 3);             rtol = 1e-13)
            @test isapprox(carter_constant(0.0, p, 0.0, 1.0),           0.0;                         atol = 1e-12)
        end
    end

    @testset "optional-arg aliases ignore passed E/Lz (KerrBackwards-wrapper signature)" begin
        a, p, e, x = 0.9, 10.0, 0.1, 0.5
        E  = geodesic_energy(a, p, e, x)
        Lz = geodesic_angular_momentum(a, p, e, x)
        @test geodesic_angular_momentum(a, p, e, x, E) == Lz
        @test carter_constant(a, p, e, x, E)     == carter_constant(a, p, e, x)
        @test carter_constant(a, p, e, x, E, Lz) == carter_constant(a, p, e, x)
    end

    @testset "exact polar endpoint x=0 (BHPToolkit formulas)" begin
        # Captured from the authoritative local BHPToolkit ConstantsOfMotion.m polar branches.
        circular = constants_of_motion(0.9, 10.0, 0.0, 0.0)
        @test circular.E  ≈ 0.9559116358210666 rtol = 2e-15
        @test circular.Lz === 0.0
        @test circular.Q  ≈ 14.175086162300348 rtol = 2e-15

        eccentric = constants_of_motion(0.9, 10.0, 0.2, 0.0)
        @test eccentric.E  ≈ 0.9574725410737306 rtol = 2e-15
        @test eccentric.Lz === 0.0
        @test eccentric.Q  ≈ 14.25031969387193 rtol = 2e-15

        # Independent turning-point check of the Kerr radial potential at peri/apoapsis.
        radial_potential(r, c) =
            (c.E * (r^2 + 0.9^2) - 0.9 * c.Lz)^2 -
            (r^2 - 2r + 0.9^2) * (r^2 + (c.Lz - 0.9 * c.E)^2 + c.Q)
        for r in (10.0 / 1.2, 10.0 / 0.8)
            @test abs(radial_potential(r, eccentric)) < 1e-10
        end

        # The exact branch joins the generic Schmidt expression continuously from both signs of x.
        for e in (0.0, 0.2)
            polar = constants_of_motion(0.9, 10.0, e, 0.0)
            for x in (-1e-6, 1e-6)
                nearby = constants_of_motion(0.9, 10.0, e, x)
                @test nearby.E ≈ polar.E rtol = 1e-8
                @test nearby.Q ≈ polar.Q rtol = 5e-7
                @test abs(nearby.Lz) < 4e-6
            end
        end
    end

    @testset "BigFloat inheritance (type-generic {T})" begin
        a, p, e, x = big"0.9", big"10.0", big"0.1", big"0.5"
        @test geodesic_energy(a, p, e, x) isa BigFloat
        @test constants_of_motion(a, p, e, x).Q isa BigFloat
        polar = constants_of_motion(a, p, e, zero(BigFloat))
        @test polar.E isa BigFloat
        @test polar.Lz isa BigFloat && iszero(polar.Lz)
        @test polar.Q isa BigFloat
    end
end

# =================================================================================================
# The REMOVABLE 0/0 of the circular closed forms on the equator at p = 2 (node
# `fix-circular-constants-p2`). This testset was written BEFORE the fix.
#
# Defect: `_angular_momentum_circular` divides by ((p−2)p − a²(x²−1)), which on the equator
# (x² = 1) is exactly p(p−2); at p = 2 the numerator vanishes too, so Float64 returned NaN, and
# Q = (1−x²)·(…)² inherited it as 0·NaN. `constants_of_motion(a, 2.0, 0, 1).Lz` was therefore NaN
# at every spin — and (a, r0) = (0.998, 2M) is arXiv:2606.04998's own fig:modes prograde
# configuration. `geodesic_energy` is NOT affected: its closed form carries no such factor
# (verified below against all three oracles, and to 1e-76 at 256 bits).
#
# Three INDEPENDENT oracles for the removed value, none of them the implementation:
#   1. the Bardeen–Press–Teukolsky (1972, ApJ 178 347, Eqs. 2.12–2.13) equatorial circular closed
#      forms, transcribed here with SIGNED a (package convention: a > 0 prograde, a < 0
#      retrograde, x = +1 always);
#   2. the four-velocity identities E = (1 − 2M/r)u^t + (2Ma/r)u^φ and
#      Lz = (r² + a² + 2Ma²/r)u^φ − (2Ma/r)u^t over the SEALED `kerr_circular_ut_uphi`;
#   3. the two-sided p → 2 limit of the UNTOUCHED Schmidt form itself, evaluated at 256 bits.
#
# r = 2M carries a circular geodesic only for a > 1/√2 = 0.70710678…: the radicand common to both
# closed forms is r^{3/2} − 3M r^{1/2} + 2a√M = 2a − √2 there, i.e. r = 2M IS the prograde photon
# orbit at a = 1/√2. Below that spin both forms must raise the ordinary `DomainError` of a
# negative radicand — never a silent NaN (pre-fix, a = 0 returned E = 0.0 with Lz = Q = NaN at
# what is its own horizon). Above it, p = 2 is a stable circular orbit for a ≳ 0.9525 and an
# unstable one (inside the ISCO) below — both are legitimate (E, Lz, Q) points.
@testset "equatorial circular constants at the removable p = 2 singularity" begin
    # oracle 1 — Bardeen–Press–Teukolsky, signed a, M = 1 (u = √r)
    function bpt_circeq(a, p)
        u = sqrt(p)
        den = u * sqrt(u) * sqrt(u^3 - 3u + 2a)
        (E = (u^3 - 2u + a) / den, Lz = (p^2 - 2a * u + a^2) / den)
    end
    # oracle 2 — the u^t/u^φ route, which never touches constants.jl
    function utuphi_circeq(a, p)
        ut, uphi = kerr_circular_ut_uphi(a, p)
        (E = (1 - 2 / p) * ut + (2a / p) * uphi,
         Lz = (p^2 + a^2 + 2a^2 / p) * uphi - (2a / p) * ut)
    end
    circular_exists(a, p) = (u = sqrt(p); u^3 - 3u + 2a > 0)

    p2_spins = (0.72, 0.8, 0.9, 0.95, 0.998, 0.9999)      # a > 1/√2: r = 2M is a circular orbit
    p2_no_orbit = (0.0, 0.5, 0.7, -0.5, -0.9, -0.998)     # r = 2M at/inside the photon orbit

    # (a) finite, and equal to both independent closed-form oracles.
    @testset "finite and Bardeen-exact at p = 2: a = $a" for a in p2_spins
        c = constants_of_motion(a, 2.0, 0.0, 1.0)
        @test isfinite(c.E) && isfinite(c.Lz) && isfinite(c.Q)
        b = bpt_circeq(a, 2.0)
        @test isapprox(c.E,  b.E;  rtol = 1e-14)          # measured max 2.3e-15 over the six spins
        @test isapprox(c.Lz, b.Lz; rtol = 1e-14)          # measured max 2.5e-15
        @test c.Q === 0.0                                 # exactly zero on the equator
        u = utuphi_circeq(a, 2.0)
        @test isapprox(c.E,  u.E;  rtol = 1e-14)
        @test isapprox(c.Lz, u.Lz; rtol = 1e-14)
        @test geodesic_energy(a, 2.0, 0.0, 1.0)           === c.E
        @test geodesic_angular_momentum(a, 2.0, 0.0, 1.0) === c.Lz
        @test carter_constant(a, 2.0, 0.0, 1.0)           === c.Q
    end

    # (a′) it is the LIMIT of the untouched generic form, approached from both sides.
    @testset "equals the two-sided p → 2 limit of the Schmidt form (256 bits): a = $a" for a in p2_spins
        setprecision(BigFloat, 256) do
            ab, one_b, zero_b = BigFloat(a), one(BigFloat), zero(BigFloat)
            lim = bpt_circeq(ab, BigFloat(2))
            for d in (big"1e-20", big"1e-25"), p in (BigFloat(2) - d, BigFloat(2) + d)
                # measured max deviation 2.0e-19 (the limit is approached linearly in |p−2|)
                @test isapprox(geodesic_energy(ab, p, zero_b, one_b),           lim.E;  rtol = 1e-14)
                @test isapprox(geodesic_angular_momentum(ab, p, zero_b, one_b), lim.Lz; rtol = 1e-14)
            end
            c = constants_of_motion(ab, BigFloat(2), zero_b, one_b)   # BigFloat inherited AT p = 2
            @test c.E  isa BigFloat && isapprox(c.E,  lim.E;  rtol = 1e-60)
            @test c.Lz isa BigFloat && isapprox(c.Lz, lim.Lz; rtol = 1e-60)
            @test c.Q  isa BigFloat && iszero(c.Q)
        end
    end

    @testset "x = −1 mirror at p = 2: a = $a" for a in p2_spins
        cm = constants_of_motion(-a, 2.0, 0.0, -1.0)      # the a → −a, Lz → −Lz image
        cp = constants_of_motion(a, 2.0, 0.0, 1.0)
        @test cm.E  === cp.E                              # E's branch is untouched and exactly even
        @test isapprox(cm.Lz, -cp.Lz; rtol = 1e-14)
        @test cm.Q === 0.0
    end

    @testset "no circular orbit at r = 2M below a = 1/√2: a = $a throws, never NaN" for a in p2_no_orbit
        @test !circular_exists(a, 2.0)
        @test_throws DomainError constants_of_motion(a, 2.0, 0.0, 1.0)
        @test_throws DomainError geodesic_angular_momentum(a, 2.0, 0.0, 1.0)
    end

    # (c) THE ORACLE THAT WOULD HAVE CAUGHT THIS: a fine scan over every radius that carries a
    # circular equatorial geodesic. The p = 2 defect sits on the a = ±0.998 curves (and on every
    # a > 1/√2), so this scan is red before the fix.
    @testset "finiteness scan over 1.1 ≤ p ≤ 12 in steps of 0.005: a = $a" for a in
            (0.0, 0.5, -0.5, 0.9, -0.9, 0.998, -0.998)
        grid = collect(1.1:0.005:12.0)
        @test 2.0 in grid
        nonfinite = 0
        scanned = 0
        for p in grid, x in (1.0, -1.0)
            circular_exists(x > 0 ? a : -a, p) || continue
            scanned += 1
            c = constants_of_motion(a, p, 0.0, x)
            (isfinite(c.E) && isfinite(c.Lz) && isfinite(c.Q)) || (nonfinite += 1)
        end
        @test scanned > 1000
        @test nonfinite == 0
    end

    # (b) BITWISE NEUTRALITY. Every other (a, p, e, x) must reproduce the PRE-FIX value under
    # `===`. Captured from commit 5243f42 in a detached worktree before this node edited
    # anything. The last three rows are p = 2 on branches the fix must NOT touch (inclined
    # x = 0.5 and x = 0.9 circular, eccentric e = 0.2 equatorial).
    @testset "bitwise-identical to the pre-fix Schmidt evaluation away from the singularity" begin
        pre_fix_pin = [
        # --- circular equatorial (e = 0, x = ±1), away from p = 2 ---
        (a = 0.0, p = 6.0, e = 0.0, x = 1.0, E = 0.9428090415820634, Lz = 3.464101615137755, Q = 0.0),
        (a = 0.0, p = 6.0, e = 0.0, x = -1.0, E = 0.9428090415820634, Lz = -3.464101615137755, Q = 0.0),
        (a = 0.0, p = 10.0, e = 0.0, x = 1.0, E = 0.9561828874675149, Lz = 3.7796447300922735, Q = 0.0),
        (a = 0.0, p = 10.0, e = 0.0, x = -1.0, E = 0.9561828874675149, Lz = -3.7796447300922735, Q = 0.0),
        (a = 0.5, p = 2.5, e = 0.0, x = 1.0, E = 1.4184247781443615, Lz = 5.406167605043682, Q = 0.0),
        (a = 0.5, p = 3.0, e = 0.0, x = 1.0, E = 0.9791813544770079, Lz = 3.2980591876185925, Q = 0.0),
        (a = 0.5, p = 6.0, e = 0.0, x = 1.0, E = 0.9296817250679497, Lz = 3.051450949862335, Q = 0.0),
        (a = 0.5, p = 6.0, e = 0.0, x = -1.0, E = 0.9625863042965734, Lz = -4.006427531603818, Q = 0.0),
        (a = 0.5, p = 10.0, e = 0.0, x = 1.0, E = 0.9537754836255025, Lz = 3.5893904590904584, Q = 0.0),
        (a = 0.5, p = 10.0, e = 0.0, x = -1.0, E = 0.9592011926920925, Lz = -4.000011893806728, Q = 0.0),
        (a = -0.5, p = 2.5, e = 0.0, x = -1.0, E = 1.4184247781443615, Lz = -5.406167605043682, Q = 0.0),
        (a = -0.5, p = 3.0, e = 0.0, x = -1.0, E = 0.9791813544770079, Lz = -3.2980591876185925, Q = 0.0),
        (a = -0.5, p = 6.0, e = 0.0, x = 1.0, E = 0.9625863042965734, Lz = 4.006427531603818, Q = 0.0),
        (a = -0.5, p = 6.0, e = 0.0, x = -1.0, E = 0.9296817250679497, Lz = -3.051450949862335, Q = 0.0),
        (a = -0.5, p = 10.0, e = 0.0, x = 1.0, E = 0.9592011926920925, Lz = 4.000011893806728, Q = 0.0),
        (a = -0.5, p = 10.0, e = 0.0, x = -1.0, E = 0.9537754836255025, Lz = -3.5893904590904584, Q = 0.0),
        (a = 0.9, p = 1.9, e = 0.0, x = 1.0, E = 0.8841158772165063, Lz = 2.249118203605822, Q = 0.0),
        (a = 0.9, p = 2.5, e = 0.0, x = 1.0, E = 0.8463300697960878, Lz = 2.1095807456750633, Q = 0.0),
        (a = 0.9, p = 3.0, e = 0.0, x = 1.0, E = 0.8606311772539785, Lz = 2.188259195548883, Q = 0.0),
        (a = 0.9, p = 6.0, e = 0.0, x = 1.0, E = 0.9225996262796217, Lz = 2.7942783614832116, Q = 0.0),
        (a = 0.9, p = 6.0, e = 0.0, x = -1.0, E = 0.9853496579901126, Lz = -4.564553289908999, Q = 0.0),
        (a = 0.9, p = 10.0, e = 0.0, x = 1.0, E = 0.9522402386495982, Lz = 3.457299296190152, Q = 0.0),
        (a = 0.9, p = 10.0, e = 0.0, x = -1.0, E = 0.962112819266394, Lz = -4.199774823890682, Q = 0.0),
        (a = -0.9, p = 1.9, e = 0.0, x = -1.0, E = 0.8841158772165063, Lz = -2.249118203605822, Q = 0.0),
        (a = -0.9, p = 2.5, e = 0.0, x = -1.0, E = 0.8463300697960878, Lz = -2.1095807456750633, Q = 0.0),
        (a = -0.9, p = 3.0, e = 0.0, x = -1.0, E = 0.8606311772539785, Lz = -2.188259195548883, Q = 0.0),
        (a = -0.9, p = 6.0, e = 0.0, x = 1.0, E = 0.9853496579901126, Lz = 4.564553289908999, Q = 0.0),
        (a = -0.9, p = 6.0, e = 0.0, x = -1.0, E = 0.9225996262796217, Lz = -2.7942783614832116, Q = 0.0),
        (a = -0.9, p = 10.0, e = 0.0, x = 1.0, E = 0.962112819266394, Lz = 4.199774823890682, Q = 0.0),
        (a = -0.9, p = 10.0, e = 0.0, x = -1.0, E = 0.9522402386495982, Lz = -3.457299296190152, Q = 0.0),
        (a = 0.998, p = 1.9, e = 0.0, x = 1.0, E = 0.7673689196223052, Lz = 1.6546293484694026, Q = 0.0),
        (a = 0.998, p = 2.5, e = 0.0, x = 1.0, E = 0.8193693044494106, Lz = 1.8737109831683405, Q = 0.0),
        (a = 0.998, p = 3.0, e = 0.0, x = 1.0, E = 0.8477143209755932, Lz = 2.0303872428276937, Q = 0.0),
        (a = 0.998, p = 6.0, e = 0.0, x = 1.0, E = 0.9212350826689453, Lz = 2.739722458617393, Q = 0.0),
        (a = 0.998, p = 6.0, e = 0.0, x = -1.0, E = 0.9921792037808781, Lz = -4.722477702716197, Q = 0.0),
        (a = 0.998, p = 10.0, e = 0.0, x = 1.0, E = 0.9519128613708321, Lz = 3.4275270765924644, Q = 0.0),
        (a = 0.998, p = 10.0, e = 0.0, x = -1.0, E = 0.9628998146930459, Lz = -4.252091602776, Q = 0.0),
        (a = -0.998, p = 1.9, e = 0.0, x = -1.0, E = 0.7673689196223052, Lz = -1.6546293484694026, Q = 0.0),
        (a = -0.998, p = 2.5, e = 0.0, x = -1.0, E = 0.8193693044494106, Lz = -1.8737109831683405, Q = 0.0),
        (a = -0.998, p = 3.0, e = 0.0, x = -1.0, E = 0.8477143209755932, Lz = -2.0303872428276937, Q = 0.0),
        (a = -0.998, p = 6.0, e = 0.0, x = 1.0, E = 0.9921792037808781, Lz = 4.722477702716197, Q = 0.0),
        (a = -0.998, p = 6.0, e = 0.0, x = -1.0, E = 0.9212350826689453, Lz = -2.739722458617393, Q = 0.0),
        (a = -0.998, p = 10.0, e = 0.0, x = 1.0, E = 0.9628998146930459, Lz = 4.252091602776, Q = 0.0),
        (a = -0.998, p = 10.0, e = 0.0, x = -1.0, E = 0.9519128613708321, Lz = -3.4275270765924644, Q = 0.0),
        # --- other branches (inclined circular, eccentric equatorial, polar) ---
        (a = 0.9, p = 10.0, e = 0.0, x = 0.5, E = 0.9538433478686694, Lz = 1.7955565292763689, Q = 9.726855841624962),
        (a = 0.5, p = 7.0, e = 0.0, x = 0.3, E = 0.9419868559907482, Lz = 1.0146170980765363, Q = 10.434491975785827),
        (a = -0.7, p = 12.0, e = 0.4, x = 0.6, E = 0.9669180949945181, Lz = 2.389883252673195, Q = 10.174258138914741),
        (a = 0.9, p = 10.0, e = 0.1, x = 0.5, E = 0.9542699747531124, Lz = 1.7964716981097308, Q = 9.736223241613233),
        (a = 0.998, p = 6.0, e = 0.2, x = 0.8, E = 0.9263968995765239, Lz = 2.2696367957704515, Q = 2.9484187821323125),
        (a = 0.0, p = 8.0, e = 0.3, x = 1.0, E = 0.95253975012074, Lz = 3.6103494592521925, Q = 0.0),
        (a = 0.9, p = 10.0, e = 0.0, x = 0.0, E = 0.9559116358210666, Lz = 0.0, Q = 14.175086162300348),
        (a = 0.9, p = 10.0, e = 0.2, x = 0.0, E = 0.9574725410737306, Lz = 0.0, Q = 14.25031969387193),
        (a = 0.998, p = 2.0, e = 0.0, x = 0.5, E = 1.0124596764364187, Lz = 1.601177715351133, Q = 7.672579429669359),
        (a = 0.998, p = 2.0, e = 0.2, x = 1.0, E = 0.7859903523518139, Lz = 1.7057174336893255, Q = 0.0),
        (a = 0.9, p = 2.0, e = 0.0, x = 0.9, E = 0.9069632892440683, Lz = 2.1923468967595356, Q = 1.1547281896117398),        ]
        @test length(pre_fix_pin) == 55
        for r in pre_fix_pin
            @test geodesic_energy(r.a, r.p, r.e, r.x)           === r.E
            @test geodesic_angular_momentum(r.a, r.p, r.e, r.x) === r.Lz
            @test carter_constant(r.a, r.p, r.e, r.x)           === r.Q
        end
    end

    # p = 0 also zeroes the equatorial denominator p(p−2), but it is NOT a removable root — the
    # numerator does not vanish there, and the ring singularity is no orbit's radius. It must keep
    # the generic form's NaN rather than gain a spurious finite value from the Bardeen branch
    # (which would return Lz = Inf, Q = 0 at a ≠ 0).
    @testset "p = 0 is excluded: the other root of p(p−2) is not removable" begin
        for a in (0.0, 0.9, -0.9), p in (0.0, -0.0), x in (1.0, -1.0)
            c = constants_of_motion(a, p, 0.0, x)
            @test isnan(c.E) && isnan(c.Lz) && isnan(c.Q)
        end
    end

    # Float64 continuity through p = 2 from both sides. NOTE the generic Schmidt form loses about
    # log10|p−2|⁻¹ digits in a NEIGHBOURHOOD of the singularity (measured relative error in Lz vs
    # a 256-bit evaluation, a = 0.998: 1.1e−14 at |p−2| = 1e−2, 6.2e−11 at 1e−6, 1.5e−6 at 1e−10)
    # — an accuracy-only defect of the untouched generic form, not a NaN, and out of this node's
    # scope; closing it would mean re-evaluating every equatorial radius. Hence the tolerance
    # below is the physical variation 100·|p−2| (measured budget use 0.196), not machine epsilon.
    @testset "continuous through p = 2 in Float64: a = $a" for a in p2_spins
        c2 = constants_of_motion(a, 2.0, 0.0, 1.0)
        for d in (1e-4, 1e-6), p in (2.0 - d, 2.0 + d)
            c = constants_of_motion(a, p, 0.0, 1.0)
            @test isapprox(c.E,  c2.E;  rtol = 100d)
            @test isapprox(c.Lz, c2.Lz; rtol = 100d)
        end
    end

    # Everything downstream of the constants at p = 2 recovers with them: pre-fix, `polar_roots`,
    # `radial_roots` (r3, r4), `proper_frequency_factor` and `kerr_frequencies(…; Time="Proper")`
    # all returned NaN at (a, p) = (0.998, 2) although the Mino/BL circ-eq branch — which never
    # touches the constants — was finite there.
    @testset "downstream consumers finite at p = 2 (a = 0.998)" begin
        a, p = 0.998, 2.0
        @test all(isfinite, abs.(radial_roots(a, p, 0.0, 1.0)))
        zp, zm = polar_roots(a, p, 0.0, 1.0)
        @test isfinite(zp) && zm == 0.0
        # on the equator zp IS Υ_θ of the circ-eq branch: an independent cross-check that the
        # restored Lz is the right one (measured agreement 1.2e-16).
        @test isapprox(zp, kerr_frequencies(a, p, 0.0, 1.0; Time = "Mino").Uth; rtol = 1e-14)
        # ⟨dτ/dλ⟩ = Σ = r² + a²cos²θ, which on a circular equatorial orbit is exactly r₀²
        @test isapprox(proper_frequency_factor(a, p, 0.0, 1.0), p^2; rtol = 1e-12)
        @test all(isfinite, values(kerr_frequencies(a, p, 0.0, 1.0; Time = "Proper")))
    end
end
