# test_puncture_derivatives.jl — ESSF node `puncture-derivatives` (schedule W6-C re-run of
# W5-E; GLOBAL_DAG row 35), written BEFORE the implementation (TDD). Target:
# Puncture/derivatives.jl — the analytic FIRST derivatives of the m-mode puncture
# (arXiv:2606.04998 main.tex:481-487): ∂_r and ∂_θ by the chain rule through
# (Δr, Δθ, ϱ²), with the trivial t/φ̄ derivatives −imΩ / +im exposed on the same API; the
# ∂J_n/∂ϱ² primitive over the exact power reduction ∂R^n/∂ϱ² = (n/2)R^{n−2}
# (main.tex:483 "∇_α ϱ^n = (n/2)ϱ^{n−2}∇_α ϱ²"), with the d₃P/dz recurrence form
# (main.tex:486) as the mandated cross-oracle.
#
# Gates (GLOBAL_DAG row 35 + worker contract; calibrated by the PRE-test prototype
# re-measurement in the node evidence, proto_rerun.out — the original W5-E worker's
# calibration outputs were lost with its stream, so every number below is re-measured):
#   * ∂J_n/∂ϱ² power-reduction form vs the d₃P/dz-recurrence form (legendre_P3_dz over
#     Bourg Eq. 115): Float64 window-gated — far C1 ≥ 1e-2: 1e-10 (measured 9.0e-12; the
#     RECURRENCE form itself loses digits far out, z̃ → 1⁺ so (z̃²−1)⁻¹ amplifies —
#     the production form is separately quadrature-gated), mid C1 ∈ [1e-6, 1e-2): 5e-9
#     (measured 6.3e-10); NEAR the worldline the recurrence form cancels ~z̃²
#     (measured up to 4.1e-3 at C1 = 1e-13) so there the two forms are compared in
#     BigFloat(256) at 1e-40 (measured 2.9e-65) and the production (power-reduction)
#     form is gated vs BigFloat truth at 1e-12 (measured 6.4e-16);
#   * derivative-recurrence self-consistency vs the BAKED BigSpecial legendre_P3 fixture
#     values (mpmath dps=60, 40 digits): d₃P_n/dz = [(n−m+1)₃P_{n+1} − (n+1)z₃P_n]/(z²−1)
#     assembled from fixture value pairs (n, n+1) vs legendre_P3_dz, 1e-12 (mandated);
#   * n = −9 closure: the derivative range extends mode_integral_Rn to n = −9 which the
#     W4-D quadrature identity did not cover ([−7, 9]); trapezoid identity with the W4-D
#     noise floor (measured 1.5e-3 in gate units);
#   * bare ∂I_m/∂Δr, ∂I_m/∂Δθ vs BOTH baked effsource_calc_m first-derivative fixture
#     slots (READ the fixture headers: NaN sentinels mark exactly the 5 mixed-second-
#     derivative slots {tr, tθ, rθ, rφ, θφ} = rows 18-21, 26-29, 32-33, which are NEVER
#     read — asserted below), at the capture driver's effective offsets, window-gated per
#     measurement (comoving near/mid/far 3e-9/1e-11/1e-11, measured 5.1e-10/3.7e-13/
#     5.8e-13; elliptic 1e-8/1e-9/3e-8-m≤1, measured 1.5e-9/7.7e-11/4.9e-9; elliptic
#     far-field m ≥ 2 excluded, measured deviation O(1) = the header-documented elliptic
#     E/K far-field digit loss, adjudicated in W4-D);
#   * t/φ̄ bridge: fixture dPhiS_dt / dPhiS_dph columns equal −imΩ / +im times the value
#     columns — the SAME trivial factors the Ψ-level API exposes;
#   * Ψ-level ∂_r, ∂_θ vs Richardson-extrapolated central FD of the LANDED (W4-D)
#     puncture_psi_m at 30 points × 4 configs incl. a < 0, gate 1e-9 (mandated; measured
#     worst r 3.8e-12, θ 3.5e-10);
#   * conjugate symmetry dΨ^P_{−m} = conj(dΨ^P_m) (all five slots), φ̄_p rotation, and
#     BigFloat(256) near-worldline loss ≤ 1e-12 (measured ≤ 1.6e-14) with full type
#     inheritance.
using Test
using ScalarSelfForce
using ScalarSelfForce.BigSpecial: legendre_P3, legendre_P3_dz, pochhammer

# fixture constants may already be loaded by earlier test files under Pkg.test; guard the
# includes so this file also runs standalone.
if !@isdefined(EFFSOURCE_CONFIGS)
    include(joinpath(@__DIR__, "fixtures", "effsource", "effsource_configs.jl"))
end
if !@isdefined(EFFSOURCE_CALC_M_COMOVING)
    include(joinpath(@__DIR__, "fixtures", "effsource", "calc_m_comoving.jl"))
end
if !@isdefined(EFFSOURCE_CALC_M_ELLIPTIC)
    include(joinpath(@__DIR__, "fixtures", "effsource", "calc_m_elliptic.jl"))
end
if !@isdefined(LEGENDRE_P3_FIXTURE)   # BigSpecial's baked mpmath fixture (W1-D)
    include(joinpath(pkgdir(ScalarSelfForce.BigSpecial), "test", "fixtures_legendre_p3.jl"))
end

# ---- NaN-sentinel read discipline (W2-D header contract) ---------------------------------
# calc_m rows this node is allowed to touch: meta 1-5, PhiS 6-7, dPhiS(t,r,θ,φ) 8-15.
const _PD_READ_COLS = 1:15
# the 5 mixed second-derivative slots the C code never provides (kerr-circular.c:709-713):
# (tr, tθ, rθ, rφ, θφ) × (re, im) in the 16-35 d2 block = rows 18-21, 26-29, 32-33.
const _PD_NAN_SLOTS = Set([18, 19, 20, 21, 26, 27, 28, 29, 32, 33])
# ALL fixture reads in this file go through this accessor.
function _pd_read(fix, col, k)
    @assert k in _PD_READ_COLS "read of non-whitelisted calc_m slot $k"
    return fix[k, col]
end

# ---- in-test oracles ---------------------------------------------------------------------

# the d₃P/dz-recurrence form of ∂J_n/∂ϱ² (main.tex:481-487 differentiated verbatim):
# J_n = 2π(−1)^m/(1+ν)_m g^ν ₃P_ν^m(z̃), ν = n/2, g = ϱ√(ϱ²+z_c²), z̃ = (ϱ²+z_c²/2)/g ⇒
# ∂J_n/∂ϱ² = K g^ν [ν(2ϱ²+z_c²)/(2g²) ₃P − z_c⁴/(4g³) d₃P/dz̃]  (∂z̃/∂ϱ² = −z_c⁴/(4g³)).
function _pd_dRn_rec(n::Integer, m::Integer, rho2::Real, zc2::Real)
    T = float(promote_type(typeof(rho2), typeof(zc2)))
    u, be = T(rho2), T(zc2)
    g = sqrt(u) * sqrt(u + be)
    zt = (u + be / 2) / g
    nu = T(n) / 2
    K = 2 * T(pi) * (iseven(m) ? one(T) : -one(T)) / pochhammer(1 + nu, m)
    return K * g^nu * (nu * (2 * u + be) / (2 * g^2) * legendre_P3(nu, T(m), zt) -
                       be^2 / (4 * g^3) * legendre_P3_dz(nu, T(m), zt))
end

# trapezoid of ∫ R^n e^{-imΔφ} dΔφ (spectrally convergent for the periodic integrand)
function _pd_Rn_quad(n::Integer, m::Integer, rho2::Float64, zc2::Float64; k::Int = 18)
    N = 2^k
    s = complex(0.0)
    for j in 0:(N - 1)
        phi = -pi + 2 * pi * j / N
        s += (rho2 + zc2 * sin(phi / 2)^2)^(n / 2) * cis(-m * phi)
    end
    return s * (2 * pi / N)
end

# Richardson-extrapolated central difference: three dyadic steps, O(h^6) truncation.
function _pd_richardson(f, x, h)
    D(hh) = (f(x + hh) - f(x - hh)) / (2 * hh)
    d1 = D(h); d2 = D(h / 2); d3 = D(h / 4)
    r1 = (4 * d2 - d1) / 3
    r2 = (4 * d3 - d2) / 3
    return (16 * r2 - r1) / 15
end

# capture-driver effective offsets (generator/capture.c:290; same as test_mmode_puncture.jl)
_pd_eff(r0, dr, dth) = ((r0 + dr) - r0, (pi / 2 + dth) - pi / 2)

@testset "Puncture: analytic first derivatives (puncture-derivatives)" begin

    orbits = Dict(c.idx => CircularOrbit(c.a, c.r0) for c in EFFSOURCE_CONFIGS)
    puncs = Dict(c.idx => Puncture(orbits[c.idx]) for c in EFFSOURCE_CONFIGS)
    r0s = Dict(c.idx => c.r0 for c in EFFSOURCE_CONFIGS)
    phi_p = 0.31                       # the frozen-grid particle azimuth (fixture headers)

    @testset "NaN sentinel discipline: mixed-d2 slots never read" begin
        for fix in (EFFSOURCE_CALC_M_COMOVING, EFFSOURCE_CALC_M_ELLIPTIC)
            nanset = Set(k for k in 1:size(fix, 1) for c in 1:size(fix, 2) if isnan(fix[k, c]))
            @test nanset == _PD_NAN_SLOTS               # exactly the 5 mixed d2 slots
            @test isempty(intersect(Set(_PD_READ_COLS), nanset))   # our reads skip exactly those
        end
    end

    @testset "mode_integral_Rn_drho2: power reduction vs d₃P/dz recurrence + quadrature" begin
        # Float64 cross-form identity, window-gated per the measured recurrence-form
        # conditioning (proto_rerun.out): far 1e-10 (measured 9.0e-12 — the recurrence
        # oracle's own z̃ → 1⁺ far-field loss), mid 5e-9 (measured 6.3e-10).
        for be in (457.14285714285717, 24.294546275416153), lgr in (-6, -4, -2, 0, 2),
            n in -7:2:3, m in (0, 1, 2, 5, 10, 20)

            u = be * exp10(lgr)
            gate = lgr >= -2 ? 1e-10 : 5e-9
            a = mode_integral_Rn_drho2(n, m, u, be)
            b = _pd_dRn_rec(n, m, u, be)
            @test abs(a - b) <= gate * abs(a)
        end
        # near-worldline: the two analytic forms are the SAME function — BigFloat(256)
        # agreement 1e-40 (measured 2.9e-65); the production form additionally keeps
        # Float64 machine accuracy there (vs BigFloat truth 1e-12, measured 6.4e-16),
        # while the recurrence form cancels ~z̃² (documented, NOT gated tight here).
        setprecision(BigFloat, 256) do
            for be in (457.14285714285717, 24.294546275416153), lgr in (-13, -8, -3, 2),
                n in (-7, -3, 1), m in (0, 5, 20)

                u = be * exp10(lgr)
                aB = mode_integral_Rn_drho2(n, m, big(u), big(be))
                bB = _pd_dRn_rec(n, m, big(u), big(be))
                @test abs(aB - bB) <= big"1e-40" * abs(aB)
                if lgr < -6
                    a64 = mode_integral_Rn_drho2(n, m, u, be)
                    @test abs(a64 - aB) <= big"1e-12" * abs(aB)
                end
            end
        end
        # n = −9 quadrature closure (W4-D tested the Rn identity on [−7, 9] only; the
        # derivative range needs −9). W4-D noise floor: 1e-12·|I_0| for the Float64 sum.
        for be in (457.14285714285717, 24.294546275416153), ratio in (0.1, 1.0, 10.0),
            m in (0, 2, 5)

            rho2 = ratio^2 * be
            q = _pd_Rn_quad(-9, m, rho2, be)
            floor0 = 1e-12 * abs(_pd_Rn_quad(-9, 0, rho2, be))
            v = mode_integral_Rn(-9, m, rho2, be)
            @test abs(v - real(q)) <= max(1e-9 * abs(q), floor0)
        end
    end

    @testset "d₃P/dz recurrence vs BAKED BigSpecial fixture values (1e-12)" begin
        # fixture degree pairs (n, n+1) available at n ∈ {−0.5, +0.5}: assemble the
        # recurrence RHS from the 40-digit mpmath values and compare legendre_P3_dz.
        rows = LEGENDRE_P3_FIXTURE
        byidx = Dict((r.n, r.m, r.z) => r.val for r in rows)
        npairs = 0
        setprecision(BigFloat, 256) do
            for r in rows
                r.n in (-0.5, 0.5) || continue
                up = get(byidx, (r.n + 1, r.m, r.z), nothing)
                up === nothing && continue
                npairs += 1
                nB, mB, zB = big(r.n), big(r.m), big(r.z)
                Pn = parse(BigFloat, r.val)
                Pn1 = parse(BigFloat, up)
                rhs = ((nB - mB + 1) * Pn1 - (nB + 1) * zB * Pn) / (zB^2 - 1)
                lhs = legendre_P3_dz(nB, mB, zB)
                @test abs(lhs - rhs) <= big"1e-12" * max(abs(lhs), abs(rhs))
                # Float64 leg away from the large-z̃ cancellation (z ≤ 10): same identity
                # from Float64-rounded fixture values, still 1e-12.
                if r.z <= 10.0
                    rhs64 = ((r.n - r.m + 1) * Float64(Pn1) - (r.n + 1) * r.z * Float64(Pn)) /
                            (r.z^2 - 1)
                    lhs64 = legendre_P3_dz(r.n, Float64(r.m), r.z)
                    # term-magnitude floor tolerates a (rare) small-derivative cancellation
                    tmag = (abs((r.n - r.m + 1) * Float64(Pn1)) +
                            abs((r.n + 1) * r.z * Float64(Pn))) / (r.z^2 - 1)
                    @test abs(lhs64 - rhs64) <= 1e-12 * max(abs(lhs64), abs(rhs64), 1e-2 * tmag)
                end
            end
        end
        @test npairs == 90   # 2 degrees × 5 orders × 9 z values, all present
    end

    @testset "bare ∂I_m/∂Δr, ∂I_m/∂Δθ vs calc_m first-derivative fixtures (both variants)" begin
        # gates per window C1 = ϱ²/z_c² (measured worst in proto_rerun.out, quoted in the
        # header): the comoving near-worldline number is dominated by the C fixture's OWN
        # coordinate rounding (W3-D capture-driver diagnosis, error ledger 2026-07-28);
        # elliptic far-field m ≥ 2 is the E/K route's documented digit loss (excluded;
        # adjudicated in W4-D — the comoving/3P route tracks BigFloat truth there).
        for (variant, fix, gates) in
            (("comoving", EFFSOURCE_CALC_M_COMOVING, (3e-9, 1e-11, 1e-11)),
             ("elliptic", EFFSOURCE_CALC_M_ELLIPTIC, (1e-8, 1e-9, 3e-8)))
            worst = 0.0
            nskip = 0
            for col in 1:size(fix, 2)
                cfg = Int(_pd_read(fix, col, 1))
                m = Int(_pd_read(fix, col, 2))
                dr = _pd_read(fix, col, 4)
                dth = _pd_read(fix, col, 5)
                dre, dte = _pd_eff(r0s[cfg], dr, dth)
                p = puncs[cfg]
                o = p.orbit
                C1 = (o.rho2_dr2 * dre^2 + o.rho2_dtheta2 * dte^2) / o.zc2
                win = C1 < 1e-8 ? 1 : (C1 < 1e-2 ? 2 : 3)
                if variant == "elliptic" && win == 3 && m >= 2
                    nskip += 1
                    continue
                end
                rot = cis(m * phi_p)   # unwind the fixture's e^{−imφ_p} (real-type pairs)
                vr = real(complex(_pd_read(fix, col, 10), _pd_read(fix, col, 11)) * rot)
                vt = real(complex(_pd_read(fix, col, 12), _pd_read(fix, col, 13)) * rot)
                d = puncture_dphi_m(p, m, dre, dte)
                scale = max(abs(vr), abs(vt))
                gate = gates[win]
                @test abs(d.dr - vr) <= gate * scale
                if dte == 0
                    # θ-evenness of I_m: the derivative is an EXACT zero on the equator
                    @test d.dtheta == 0
                    @test vt == 0
                else
                    @test abs(d.dtheta - vt) <= gate * scale
                end
                worst = max(worst, max(abs(d.dr - vr), abs(d.dtheta - vt)) / scale)
            end
            println("puncture-derivatives vs ", variant, " calc_m: worst gated rel dev = ",
                    worst, (nskip > 0 ? " (skipped $nskip far-field m ≥ 2 elliptic rows)" : ""))
            variant == "elliptic" && @test nskip > 0
        end
    end

    @testset "t/φ̄ bridge: fixture dPhiS_dt, dPhiS_dph = −imΩ, +im × value columns" begin
        fix = EFFSOURCE_CALC_M_COMOVING
        for col in 1:size(fix, 2)
            cfg = Int(_pd_read(fix, col, 1))
            m = Int(_pd_read(fix, col, 2))
            Om = orbits[cfg].Omega
            v = complex(_pd_read(fix, col, 6), _pd_read(fix, col, 7))
            dt = complex(_pd_read(fix, col, 8), _pd_read(fix, col, 9))
            dp = complex(_pd_read(fix, col, 14), _pd_read(fix, col, 15))
            if m == 0
                @test dt == 0
                @test dp == 0
            else
                @test abs(dt - (-im * m * Om * v)) <= 1e-13 * abs(dt)
                @test abs(dp - im * m * v) <= 1e-13 * abs(dp)
            end
        end
    end

    @testset "Ψ-level: analytic ∂_r, ∂_θ vs Richardson FD of the LANDED puncture_psi_m" begin
        # 10 offsets × m ∈ {0, 1, 10} = 30 points per config; configs incl. a < 0 and the
        # a = 0.998 ISCO extreme. FD steps: r-step 5e-3 of the worldline distance; θ-step
        # additionally capped by the sin^{|m|}θ ansatz-factor scale 1/|m| (calibrated).
        pairs = ((-0.07, 0.03), (0.06, -0.12), (0.15, 0.30), (0.35, -0.20), (0.80, 0.10),
                 (-0.05, 0.40), (0.10, 0.02), (0.50, -0.35), (0.25, 0.18), (1.20, -0.05))
        worst = 0.0
        npts = 0
        for cfg in (1, 4, 5, 6), m in (0, 1, 10), (frac, z) in pairs
            p = puncs[cfg]
            o = p.orbit
            r = o.r0 * (1 + frac)
            @assert r > o.rp + 0.05    # FD stencil stays exterior
            npts += 1
            d = puncture_dpsi_m(p, m, r, z; varphi_p = phi_p)
            # the value slot IS the landed mode function (same assembly, same knowledge home)
            psi_l = puncture_psi_m(p, m, r, z; varphi_p = phi_p)
            @test d.psi == psi_l
            @test d.dt == -im * m * o.Omega * d.psi
            @test d.dphi == im * m * d.psi
            rho = sqrt(o.rho2_dr2 * (r - o.r0)^2 + o.rho2_dtheta2 * asin(z)^2)
            hr = 5e-3 * rho / sqrt(o.rho2_dr2)
            ht = min(5e-3 * rho / sqrt(o.rho2_dtheta2), 1e-2 / max(1, abs(m)))
            fdr = _pd_richardson(rr -> puncture_psi_m(p, m, rr, z; varphi_p = phi_p), r, hr)
            th0 = acos(z)
            fdt = _pd_richardson(th -> puncture_psi_m(p, m, r, cos(th); varphi_p = phi_p),
                                 th0, ht)
            er = abs(d.dr - fdr) / max(abs(d.dr), abs(d.psi))
            et = abs(d.dtheta - fdt) / max(abs(d.dtheta), abs(d.psi))
            @test er <= 1e-9
            @test et <= 1e-9
            worst = max(worst, max(er, et))
        end
        @test npts == 120
        println("puncture-derivatives Ψ-level vs Richardson FD: 120 points, worst rel dev = ",
                worst)
    end

    @testset "conjugate symmetry, φ̄_p rotation, m < 0" begin
        for cfg in (1, 4, 6), m in (1, 2, 10), (frac, z) in ((0.15, 0.30), (-0.05, -0.12))
            p = puncs[cfg]
            r = p.orbit.r0 * (1 + frac)
            dp_ = puncture_dpsi_m(p, m, r, z; varphi_p = phi_p)
            dm = puncture_dpsi_m(p, -m, r, z; varphi_p = phi_p)
            for slot in (:psi, :dt, :dr, :dtheta, :dphi)
                a = getfield(dp_, slot)
                b = getfield(dm, slot)
                @test abs(b - conj(a)) <= 1e-12 * max(abs(a), 1e-300)
            end
            # φ̄_p enters ALL slots as the common rotation e^{−imφ̄_p}
            d0 = puncture_dpsi_m(p, m, r, z)
            for slot in (:psi, :dt, :dr, :dtheta, :dphi)
                @test abs(getfield(dp_, slot) - getfield(d0, slot) * cis(-m * phi_p)) <=
                      1e-14 * max(abs(getfield(d0, slot)), 1e-300)
            end
            # bare integral derivatives are REAL and EVEN in m
            db = puncture_dphi_m(p, m, 0.15 * p.orbit.r0, -asin(z))
            dbm = puncture_dphi_m(p, -m, 0.15 * p.orbit.r0, -asin(z))
            @test db.dr isa Real && db.dtheta isa Real
            @test abs(dbm.dr - db.dr) <= 1e-12 * abs(db.dr)
            @test abs(dbm.dtheta - db.dtheta) <= 1e-12 * abs(db.dtheta)
        end
    end

    @testset "BigFloat(256): near-worldline accuracy + type inheritance" begin
        setprecision(BigFloat, 256) do
            cB = EFFSOURCE_CONFIGS[4]                       # a = 0.998 AT the ISCO
            oB = CircularOrbit(big(cB.a), big(cB.r0))
            pB = Puncture(oB)
            o64 = orbits[4]
            p64 = puncs[4]
            dr, z = 1e-6, 5.0e-7                            # deep near-worldline point
            m = 10
            r64 = o64.r0 + dr
            # SAME real-number field point in both precisions: the Float64 leg sits at the
            # EFFECTIVE offset fl(r0 + dr) − r0 (W3-D capture-driver diagnosis), so the
            # BigFloat leg must too — big(r0) + big(dr) would be a slightly different point
            # (relative dr mismatch ~3e-10) and would masquerade as precision loss.
            dr_eff = (o64.r0 + dr) - o64.r0
            rB = oB.r0 + big(dr_eff)
            dB = puncture_dpsi_m(pB, m, rB, big(z); varphi_p = big(phi_p))
            d64 = puncture_dpsi_m(p64, m, r64, z; varphi_p = phi_p)
            @test dB.psi isa Complex{BigFloat}
            @test dB.dr isa Complex{BigFloat}
            for slot in (:psi, :dt, :dr, :dtheta, :dphi)
                a = getfield(d64, slot)
                b = getfield(dB, slot)
                @test abs(a - b) <= big"1e-12" * abs(b)     # F64 loss (measured ≤ 1.6e-14)
            end
            bB = puncture_dphi_m(pB, m, big(dr), big(-5.0e-7))
            @test bB.dr isa BigFloat && bB.dtheta isa BigFloat
            b64 = puncture_dphi_m(p64, m, dr, -5.0e-7)
            @test abs(b64.dr - bB.dr) <= big"1e-12" * abs(bB.dr)
            @test abs(b64.dtheta - bB.dtheta) <= big"1e-12" * abs(bB.dtheta)
            # ∂J/∂ϱ² BigFloat inheritance
            @test mode_integral_Rn_drho2(-7, 3, big"0.04", big"24.5") isa BigFloat
        end
    end
end
