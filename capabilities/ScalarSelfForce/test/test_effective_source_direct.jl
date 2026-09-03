# test_effective_source_direct.jl — ESSF node `effective-source-direct` (schedule W7-C,
# GLOBAL_DAG row 36), written BEFORE the implementation (TDD). Target:
# `EffectiveSource/direct.jl` + `Puncture/second_derivatives.jl` — the POINTWISE direct
# evaluation of the m-mode effective source
#
#     S^eff_m(r, z) = the m-mode of −□Φ^P  =  −△_m Ψ^P_m        (main.tex:414-460)
#
# in exactly the conventions the SEALED flux-form operator consumes, plus the complete m-mode
# convention transformation chain that maps the baked C `effsource_calc_m` oracle onto it.
#
# THE CHAIN (derived here, gated below, stated once in the direct.jl docstring):
#   □Φ = −4πρ (eq:KG) and △̃_mΨ̃_m = S_m (eq:elliptic_eqn, main.tex:247-262) together fix the
#   mode-projection scale: for ANY field F, △̃_m[mode of F] = −(Σ/(r²+a²))·P_m[□F] with
#   P_m[f] = r e^{imΩt}/(2π sin^{|m|}θ)·∫f e^{−imφ̄}dφ̄ (eq:modeIntegral). With
#   □Φ^R = S^eff = −4πρ − □Φ^P (main.tex:428-432) this gives, off the worldline,
#       S^eff_m := −(Σ/(r²+a²))·P_m[S^eff] = +(Σ/(r²+a²))·P_m[□Φ^P] = −△_m Ψ^P_m ,
#   the spectre identity (Test_CircularOrbit.cpp:99-112, `-div F + sources == -effective_source`).
#   In the horizon-regular chart φ̄ = φ + χ(r) the mode function carries the FIELD-POINT
#   rotation e^{−imχ(r)} (χ = kerr_varphi_offset, χ′ = a/Δ), so
#       Ψ^P_m(r, z) = r e^{−im(φ_p + χ(r))} / (2π sin^{|m|}θ) · I_m(Δr, Δθ)
#   and, with `src_C` the C oracle's `effsource_calc_m` src (which carries its own e^{−imφ_p}
#   and equals +□_m I_m — verified below),
#       S^eff_m = e^{−imχ(r)} · r Σ / (2π (r²+a²) sin^{|m|}θ) · src_C
#   — transcribed INDEPENDENTLY below from spectre CircularOrbit.cpp:322-332 (`rotation * 0.5 *
#   r / M_PI`, then `* Σ/((r²+a²) sin^m θ)`; the extra `* Δ/(r²+a²)` there belongs to the
#   NON-horizon-penetrating r_* chart and does NOT apply to our r chart).
#
# GATE CALIBRATION (PRE-RED, evidence dir calibrate_effsource{,2,4,5,6,7}.jl +
# calibration*.log; never tuned post hoc). Three instrument facts drive every constant:
#   * BigFloat here costs ~1000× Float64 (26 µs vs 28 ms at 128 bits, 138 ms at 256), so the
#     BigFloat legs are small and budgeted;
#   * the direct assembly cancels FOUR orders of singularity (main.tex:488-490), so its
#     roundoff amplification is ~λ^{-2} — measured 1e3…2e6 on the near-puncture identity box
#     and 1e5…1.7e10 at λ ∈ [1e-2, 1e-1) on the fixture grid. A 1e-10 RELATIVE identity gate
#     is therefore IMPOSSIBLE in Float64 (measured 1e-4…1e-6); the full 20×20 grid is gated
#     in Float64 against the TERM SCALE (an assembly sign/factor error is O(1) there) and the
#     1e-10 RELATIVE gate is met at BigFloat(128) on a subgrid of the same box;
#   * the same cancellation bounds the fixture agreement: at λ ∈ [1e-2, 1e-1) BOTH the C
#     oracle and this implementation carry 1e-6…1e-3 relative error against BigFloat(256)
#     truth, so the DAG's "1e-6 at λ ≳ 1e-2" is unreachable by ANY Float64 direct evaluation.
#     The windows below are the measured envelope, and the BigFloat(256) ADJUDICATION testset
#     shows the excluded near-worldline region is the C's loss, not ours (≤ 1e-2 of the C's
#     error there) — the fixture header's own warning ("the direct Box[Phi^P] evaluation is
#     DESTROYED by cancellation near the worldline … trustworthy only for ϱ² ≳ 1e-5").
using Test
using ScalarSelfForce
using ScalarSelfForce: kerr_varphi_offset, kerr_varphi_offset_dr

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

# ---- NaN-sentinel read discipline (W2-D header contract, inherited from W6-C) ------------
# calc_m rows THIS node touches: meta 1-5 and the src pair 36-37. The 5 mixed second
# derivatives the C never provides are rows 18-21, 26-29, 32-33 and are never read (this node
# needs NO mixed second derivative: the flux form is diagonal in (r, z)).
const _ESD_READ_COLS = Set([1, 2, 3, 4, 5, 36, 37])
const _ESD_NAN_SLOTS = Set([18, 19, 20, 21, 26, 27, 28, 29, 32, 33])
function _esd_read(fix, col, k)
    @assert k in _ESD_READ_COLS "read of non-whitelisted calc_m slot $k"
    return fix[k, col]
end

# ---- in-test instruments (uniquely named — the package already has two `cdiff` in Main) ---
_esd_cd1(f, x, h) = (f(x + h) - f(x - h)) / (2 * h)
_esd_cd2(f, x, h) = (f(x + h) - 2 * f(x) + f(x - h)) / h^2
_esd_rich(f, x, h) = let d1 = _esd_cd1(f, x, h), d2 = _esd_cd1(f, x, h / 2),
                         d3 = _esd_cd1(f, x, h / 4)
    (16 * (4 * d3 - d2) / 3 - (4 * d2 - d1) / 3) / 15
end
# 8th-order central stencils (Float64 leg; h = 0.025λ is the measured optimum, part 4 (G))
const _ESD_W1 = (4 / 5, -1 / 5, 4 / 105, -1 / 280)
const _ESD_W2 = (8 / 5, -1 / 5, 8 / 315, -1 / 560)
const _ESD_W20 = -205 / 72
_esd_d1_8(f, x, h) =
    (_ESD_W1[1] * (f(x + h) - f(x - h)) + _ESD_W1[2] * (f(x + 2h) - f(x - 2h)) +
     _ESD_W1[3] * (f(x + 3h) - f(x - 3h)) + _ESD_W1[4] * (f(x + 4h) - f(x - 4h))) / h
_esd_d2_8(f, x, h) =
    (_ESD_W2[1] * (f(x + h) + f(x - h)) + _ESD_W2[2] * (f(x + 2h) + f(x - 2h)) +
     _ESD_W2[3] * (f(x + 3h) + f(x - 3h)) + _ESD_W2[4] * (f(x + 4h) + f(x - 4h)) +
     _ESD_W20 * f(x)) / h^2

# λ = ϱ, the puncture's local distance scale (main.tex:452)
_esd_lambda(o, dr, dth) = sqrt(o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dth^2)

# THE convention transformation chain, transcribed INDEPENDENTLY of the implementation from
# spectre CircularOrbit.cpp:322-332 + main.tex:210-262/414-460 (see the file header).
function _esd_chain(o, m, r, z)
    Sigma = r^2 + o.a^2 * z^2
    return cis(-m * kerr_varphi_offset(o.a, r)) * r * Sigma /
           (2 * oftype(r, pi) * (r^2 + o.a^2) * sinm_theta(m, z))
end

# in-test trapezoid of ∫R^n e^{−imΔφ}dΔφ (spectrally convergent, periodic integrand)
function _esd_Rn_quad(n::Integer, m::Integer, rho2::Float64, zc2::Float64; k::Int = 18)
    N = 2^k
    s = complex(0.0)
    for j in 0:(N - 1)
        phi = -pi + 2 * pi * j / N
        s += (rho2 + zc2 * sin(phi / 2)^2)^(n / 2) * cis(-m * phi)
    end
    return s * (2 * pi / N)
end

# the near-puncture identity box (worldline EXCLUDED): Δr ∈ [0.2, 1.2], z ∈ [0.01, 0.21],
# the same shape as spectre Test_CircularOrbit.cpp:33-42 (cosθ offset 0.1, width 0.2).
_esd_box_r(r0, i, n) = r0 + 0.2 + 1.0 * (i - 1) / (n - 1)
_esd_box_z(j, n) = 0.01 + 0.2 * (j - 1) / (n - 1)
const _ESD_IDENTITY_CFGS = ((0.0, 10.0), (0.5, 10.0), (-0.9, 10.0),
                            (0.0, 6.0), (0.5, 6.0), (-0.9, 6.0))

@testset "EffectiveSource: direct m-mode source (effective-source-direct)" begin

    orbits = Dict(c.idx => CircularOrbit(c.a, c.r0) for c in EFFSOURCE_CONFIGS)
    puncs = Dict(c.idx => Puncture(orbits[c.idx]) for c in EFFSOURCE_CONFIGS)
    phi_p = EFFSOURCE_PHIP        # 0.31, the frozen-grid particle BOYER-LINDQUIST azimuth

    @testset "NaN sentinel discipline: this node reads only meta + src columns" begin
        for fix in (EFFSOURCE_CALC_M_COMOVING, EFFSOURCE_CALC_M_ELLIPTIC)
            nanset = Set(k for k in 1:size(fix, 1) for c in 1:size(fix, 2) if isnan(fix[k, c]))
            @test nanset == _ESD_NAN_SLOTS
            @test isempty(intersect(_ESD_READ_COLS, nanset))
        end
    end

    @testset "mode_integral_Rn_d2rho2: ∂²J_n/∂(ϱ²)² = (n/2)((n−2)/2) J_{n−4}" begin
        # (a) vs a BigFloat(192) central difference of the SEALED W6-C first-derivative
        #     primitive (h = 1e-10 ϱ²: truncation ~1e-21, roundoff ~1e-48). Measured 1.1e-18.
        setprecision(BigFloat, 192) do
            worst = 0.0
            for be in (big"457.14285714285717", big"24.294546275416153"),
                lgr in (-6, -3, 0, 2), n in -7:2:1, m in (0, 1, 2, 5, 10, 20)

                u = be * big(10.0)^lgr
                a = mode_integral_Rn_d2rho2(n, m, u, be)
                b = _esd_cd1(x -> mode_integral_Rn_drho2(n, m, x, be), u, big"1e-10" * u)
                @test abs(a - b) <= big"1e-14" * abs(a)
                worst = max(worst, Float64(abs(a - b) / abs(a)))
            end
            println("effective-source-direct d2rho2 vs BigFloat FD of sealed drho2: worst = ",
                    worst)
        end
        # (b) n = −11 range closure: the SECOND derivative pushes `mode_integral_Rn` two orders
        #     below W6-C's n = −9 (itself two below W4-D's n = −7). Trapezoid identity at the
        #     W4-D noise floor (measured 9.9e-4 in gate units).
        for be in (457.14285714285717, 24.294546275416153), ratio in (0.1, 1.0, 10.0),
            m in (0, 2, 5)

            rho2 = ratio^2 * be
            q = _esd_Rn_quad(-11, m, rho2, be)
            floor0 = 1e-12 * abs(_esd_Rn_quad(-11, 0, rho2, be))
            v = mode_integral_Rn(-11, m, rho2, be)
            @test abs(v - real(q)) <= max(1e-9 * abs(q), floor0)
        end
        # (c) type inheritance + evenness in m. `mode_integral_Rn` reaches m < 0 through the
        #     reciprocal-falling Pochhammer + negative-order ₃P connection (its docstring's
        #     Γ-function identity), a DIFFERENT computational path — so evenness holds to the
        #     working precision's roundoff, not bitwise (measured 4.9e-78 absolute at 256 bits,
        #     eps(BigFloat) = 1.7e-77; the sealed W4-D gate is the Float64 analogue at 1e-13).
        setprecision(BigFloat, 256) do
            v = mode_integral_Rn_d2rho2(-7, 5, big"1.0", big"100.0")
            @test v isa BigFloat
            @test abs(mode_integral_Rn_d2rho2(-7, -5, big"1.0", big"100.0") - v) <=
                  10 * eps(BigFloat) * abs(v)
        end
    end

    @testset "puncture_d2phi_m: value + first-derivative slots == the sealed puncture_dphi_m" begin
        worst = 0.0
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10, 20),
            (dr, dt) in ((1e-4, 0.0), (0.01, 0.0), (0.3, 0.0), (-0.1, 0.0), (0.0, 0.01),
                         (0.0, 0.3), (0.01, 0.005), (1.0, 0.2), (3.0, -0.5))

            p = puncs[c.idx]
            a = puncture_d2phi_m(p, m, dr, dt)
            b = puncture_dphi_m(p, m, dr, dt)
            sc = max(abs(b.phi), abs(b.dr), abs(b.dtheta))
            # the value slot is BITWISE the sealed assembly (same c_j, same J, same order) —
            # measured on 324/324 calibration points
            @test a.phi == b.phi
            # the first derivatives regroup the ∂/∂ϱ² sum (the c_j ladder identity below), so
            # they agree to rounding only
            @test abs(a.dr - b.dr) <= 1e-13 * sc
            @test abs(a.dtheta - b.dtheta) <= 1e-13 * sc
            worst = max(worst, max(abs(a.dr - b.dr), abs(a.dtheta - b.dtheta)) / sc)
            # θ-evenness of I_m: the FIRST θ-derivative is an exact zero on the equator, the
            # SECOND is not (it carries the leading ∂²_θ of the source there).
            if dt == 0
                @test a.dtheta == 0
                @test a.dthetatheta != 0
            end
        end
        println("effective-source-direct d2phi_m vs sealed dphi_m 1st slots: worst = ", worst)
    end

    @testset "puncture_d2phi_m: ∂²I/∂Δr², ∂²I/∂Δθ² vs FD of the sealed first derivatives" begin
        # BigFloat leg — truncation-limited central difference of the SEALED W6-C
        # puncture_dphi_m; the primary correctness gate on the new second-derivative chain
        # (the C d2A groups + the c_j(ϱ²) ladder + the J power ladder). Measured 1.2e-18.
        setprecision(BigFloat, 192) do
            worst_r = 0.0
            worst_t = 0.0
            for c in (EFFSOURCE_CONFIGS[1], EFFSOURCE_CONFIGS[4], EFFSOURCE_CONFIGS[5])
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                for m in (0, 2, 20),
                    (dr, dt) in ((big"0.01", big(0)), (big"0.3", big"0.05"),
                                 (big(0), big"0.1"), (big(1), big"0.2"))

                    d = _esd_lambda(oB, dr, dt)
                    a = puncture_d2phi_m(pB, m, dr, dt)
                    frr = _esd_cd1(x -> puncture_dphi_m(pB, m, x, dt).dr, dr,
                                   big"1e-10" * d / sqrt(oB.rho2_dr2))
                    @test abs(a.drr - frr) <= big"1e-14" * abs(frr)
                    worst_r = max(worst_r, Float64(abs(a.drr - frr) / abs(frr)))
                    if dt != 0
                        ftt = _esd_cd1(x -> puncture_dphi_m(pB, m, dr, x).dtheta, dt,
                                       big"1e-10" * d / sqrt(oB.rho2_dtheta2))
                        @test abs(a.dthetatheta - ftt) <= big"1e-14" * abs(ftt)
                        worst_t = max(worst_t, Float64(abs(a.dthetatheta - ftt) / abs(ftt)))
                    end
                end
            end
            println("effective-source-direct I_rr/I_tt vs BigFloat FD of sealed dphi_m: ",
                    worst_r, " / ", worst_t)
        end
        # Float64 leg — the production path at the Float64 Richardson instrument's own limit
        # (measured 1.4e-9 / 3.6e-11).
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 20),
            (dr, dt) in ((0.01, 0.0), (0.1, 0.0), (0.3, 0.05), (-0.2, 0.0), (0.0, 0.1),
                         (1.0, 0.2))

            p = puncs[c.idx]
            o = p.orbit
            d = max(_esd_lambda(o, dr, dt), 1e-8)
            a = puncture_d2phi_m(p, m, dr, dt)
            frr = _esd_rich(x -> puncture_dphi_m(p, m, x, dt).dr, dr,
                            1e-2 * d / sqrt(o.rho2_dr2))
            @test abs(a.drr - frr) <= 1e-7 * abs(frr)
            if dt != 0
                ftt = _esd_rich(x -> puncture_dphi_m(p, m, dr, x).dtheta, dt,
                                1e-2 * d / sqrt(o.rho2_dtheta2))
                @test abs(a.dthetatheta - ftt) <= 1e-7 * abs(ftt)
            end
        end
    end

    @testset "convention chain: Ψ^P_m == the sealed Ψ-level API at the chart-corrected φ̄" begin
        # The ONLY difference between the chart-complete bundle and the sealed W4-D/W6-C
        # Ψ-level API is the FIELD-POINT rotation e^{−imχ(r)}: passing
        # varphi_p = kerr_varphi(a, r, φ_p) to the sealed calls reproduces psi (measured
        # 3.5e-16) and dr up to the exactly-known −im χ′(r) Ψ chain term, χ′ = a/Δ — spectre
        # CircularOrbit.cpp:338-341 adds precisely `singular_field/r − i a m/Δ · singular_field`
        # (the /r piece is its own r-scaling convention, not ours).
        wpsi = 0.0
        for c in EFFSOURCE_CONFIGS, m in (0, 1, 2, 5, 10),
            (dr, z) in ((0.3, 0.05), (1.0, 0.2), (-0.5, -0.1), (3.0, 0.4))

            p = puncs[c.idx]
            o = p.orbit
            r = o.r0 + dr
            vp = kerr_varphi(o.a, r, phi_p)
            b = puncture_psi_chart(p, m, r, z; phi_p = phi_p)
            ps = puncture_psi_m(p, m, r, z; varphi_p = vp)
            @test abs(b.psi - ps) <= 1e-14 * abs(ps)
            wpsi = max(wpsi, abs(b.psi - ps) / abs(ps))
            sealed = puncture_dpsi_m(p, m, r, z; varphi_p = vp)
            expect_dr = sealed.dr - im * m * kerr_varphi_offset_dr(o.a, r) * b.psi
            @test abs(b.dr - expect_dr) <= 1e-12 * abs(expect_dr)
            # ∂_z = −(1/sinθ)∂_θ (the locked z = cosθ chart)
            expect_dz = -sealed.dtheta / sintheta_of_z(z)
            @test abs(b.dz - expect_dz) <= 1e-12 * abs(expect_dz)
        end
        println("effective-source-direct chart Ψ^P_m vs sealed puncture_psi_m: worst = ", wpsi)
    end

    @testset "S^eff_m vs the baked effsource_calc_m src fixtures (both C variants)" begin
        # Region-restricted oracles (adopted decision 4). λ-windows are the PRE-RED measured
        # envelope (calibration7.log); λ < 1e-2 is excluded because BOTH sides lose 1e-4…1e0
        # there to the shared direct-□ cancellation — adjudicated in the next testset.
        # The elliptic (kerr-circular.c E/K-polynomial) variant is additionally restricted to
        # m ≤ 1 in the λ ≥ 1 window: main.tex:465-466's far-field ill-conditioning, which the
        # SECOND derivatives in `src` amplify far past the W4-D PhiS_m/first-derivative case
        # (measured 2.0e-5 at m = 2 and O(1) at m ≥ 5 there).
        for (variant, fix, gates, restrict) in
            (("comoving", EFFSOURCE_CALC_M_COMOVING, (1e-7, 1e-4, 1e-2), false),
             ("elliptic", EFFSOURCE_CALC_M_ELLIPTIC, (1e-8, 1e-4, 2e-2), true))

            worst = 0.0
            worst_norm = 0.0
            nskip = 0
            nwin = [0, 0, 0]
            for col in 1:size(fix, 2)
                cfg = Int(_esd_read(fix, col, 1))
                m = Int(_esd_read(fix, col, 2))
                o = orbits[cfg]
                p = puncs[cfg]
                r = o.r0 + _esd_read(fix, col, 4)     # the capture driver's effective offsets
                th = pi / 2 + _esd_read(fix, col, 5)
                z = cos(th)
                lam = _esd_lambda(o, r - o.r0, th - pi / 2)
                if lam < 1e-2 || (restrict && lam >= 1.0 && m >= 2)
                    nskip += 1
                    continue
                end
                win = lam >= 1.0 ? 1 : (lam >= 1e-1 ? 2 : 3)
                nwin[win] += 1
                srcC = complex(_esd_read(fix, col, 36), _esd_read(fix, col, 37))
                theirs = srcC * _esd_chain(o, m, r, z)
                b = effective_source_bundle(p, m, r, z; phi_p = phi_p)
                rel = abs(b.source - theirs) / max(abs(theirs), abs(b.source))
                @test rel <= gates[win]
                worst = max(worst, rel)
                worst_norm = max(worst_norm, rel / (eps(Float64) * b.cancellation))
            end
            # the disagreement is ENTIRELY explained by the shared λ^{-2} roundoff
            # amplification of the direct assembly (measured 7.0e3 comoving / 2.2e5 elliptic)
            @test worst_norm <= (variant == "comoving" ? 1e5 : 1e8)
            println("effective-source-direct vs ", variant, " calc_m src: worst rel dev = ",
                    worst, " over ", sum(nwin), " rows ", nwin, " (", nskip, " excluded); ",
                    "worst rel/(eps·cancellation) = ", worst_norm)
            @test nwin == [variant == "comoving" ? 306 : 102, 288, 276]
            @test nskip > 0
        end
    end

    @testset "BigFloat(256) adjudication of the excluded near-worldline region" begin
        # WHY the λ < 1e-2 rows are excluded: there the C oracle's own direct evaluation is
        # destroyed by the cancellation (its fixture header says so), while ours is 2-5 orders
        # closer to BigFloat(256) truth. Sampled deterministically: the first row of each
        # (config, m) bucket with λ < 1e-3.
        setprecision(BigFloat, 256) do
            fix = EFFSOURCE_CALC_M_COMOVING
            seen = Set{Tuple{Int,Int}}()
            nchk = 0
            worst_ratio = 0.0
            for col in 1:size(fix, 2)
                cfg = Int(_esd_read(fix, col, 1))
                m = Int(_esd_read(fix, col, 2))
                cfg in (1, 4) || continue
                (cfg, m) in seen && continue
                o = orbits[cfg]
                p = puncs[cfg]
                r = o.r0 + _esd_read(fix, col, 4)
                th = pi / 2 + _esd_read(fix, col, 5)
                z = cos(th)
                _esd_lambda(o, r - o.r0, th - pi / 2) < 1e-3 || continue
                push!(seen, (cfg, m))
                nchk += 1
                oB = CircularOrbit(big(o.a), big(o.r0))
                pB = Puncture(oB)
                truth = effective_source_m(pB, m, big(r), big(z); phi_p = big(phi_p))
                ours = effective_source_m(p, m, r, z; phi_p = phi_p)
                theirs = complex(_esd_read(fix, col, 36), _esd_read(fix, col, 37)) *
                         _esd_chain(o, m, r, z)
                e_ours = Float64(abs(ours - truth) / abs(truth))
                e_c = Float64(abs(theirs - truth) / abs(truth))
                @test e_ours <= 1e-2 * e_c
                worst_ratio = max(worst_ratio, e_ours / e_c)
            end
            @test nchk == 12          # 2 configs × 6 m values
            println("effective-source-direct near-worldline adjudication: worst ",
                    "ours/C error ratio over ", nchk, " points = ", worst_ratio)
        end
    end

    @testset "operator↔source identity on the 20×20 near-puncture grid (spectre twin)" begin
        # −∂_i𝓕^i(Ψ^P_m) + β_mΨ^P_m + γ^i_m∂_iΨ^P_m == −S^eff_m, i.e. △_mΨ^P_m == −S^eff_m
        # (Test_CircularOrbit.cpp:99-112). Coefficients: the SEALED fluxform_callables in the
        # t region (H = 0, the worldtube slab). Derivatives of Ψ^P_m: an INDEPENDENT
        # finite-difference instrument applied to the chart-complete Ψ^P_m — never the source
        # formula. The nested flux divergence is expanded as A∂² + (∂A)∂ with ∂A itself
        # finite-differenced from the sealed coefficient callables.
        #
        # (a) Float64, 8th-order stencils at h = 0.025λ (the measured optimum), the FULL 20×20
        #     grid at all 6 (a, r₀) × 4 m = 9600 points, gated against the TERM SCALE: an
        #     assembly sign/factor error is O(1)·term_scale, while the cancellation (1e3…2e6
        #     here) forbids a relative Float64 gate. Measured worst 3.2e-8.
        gworst_ts = 0.0
        for (a, r0) in _ESD_IDENTITY_CFGS
            o = CircularOrbit(a, r0)
            p = Puncture(o)
            v = vtu_slicing(o, r0 - 1.0, r0 + 2.0)
            for m in (0, 1, 2, 5)
                fc = fluxform_callables(o, m, v, :t)
                psi(rr, zz) = puncture_psi_chart(p, m, rr, zz; phi_p = phi_p).psi
                for i in 1:20, j in 1:20
                    r = _esd_box_r(r0, i, 20)
                    z = _esd_box_z(j, 20)
                    d = _esd_lambda(o, r - o.r0, -asin(z))
                    hr = 0.025 * d / sqrt(o.rho2_dr2)
                    hz = 0.025 * min(d / sqrt(o.rho2_dtheta2), z / 2)
                    pr(x) = psi(x, z)
                    pz(x) = psi(r, x)
                    dpr = _esd_d1_8(pr, r, hr)
                    d2pr = _esd_d2_8(pr, r, hr)
                    dpz = _esd_d1_8(pz, z, hz)
                    d2pz = _esd_d2_8(pz, z, hz)
                    dArr = _esd_d1_8(x -> fc.A[1](x, z), r, hr)
                    dAzz = _esd_d1_8(x -> fc.A[2](r, x), z, hz)
                    lhs = -(fc.A[1](r, z) * d2pr + dArr * dpr) -
                          (fc.A[2](r, z) * d2pz + dAzz * dpz) +
                          fc.beta(r, z) * psi(r, z) +
                          fc.gamma[1](r, z) * dpr + fc.gamma[2](r, z) * dpz
                    b = effective_source_bundle(p, m, r, z; phi_p = phi_p)
                    res = abs(lhs + b.source) / b.term_scale
                    @test res <= 1e-6
                    gworst_ts = max(gworst_ts, res)
                end
            end
        end
        println("effective-source-direct identity (Float64, 9600 pts): worst ",
                "|△_mΨ^P + S^eff| / term_scale = ", gworst_ts)
        # (b) BigFloat(128), 2nd-order stencils at h = 1e-9λ (truncation ~1e-19, roundoff
        #     ~1e-21), 3×3 subgrid of the SAME box at all 24 configurations = 216 points, the
        #     DAG's 1e-10 RELATIVE gate. Measured 2.2e-12 on a 2×3 subgrid.
        setprecision(BigFloat, 128) do
            gworst = 0.0
            npt = 0
            for (a, r0) in _ESD_IDENTITY_CFGS
                oB = CircularOrbit(big(a), big(r0))
                pB = Puncture(oB)
                v = vtu_slicing(oB, oB.r0 - 1, oB.r0 + 2)
                for m in (0, 1, 2, 5)
                    fc = fluxform_callables(oB, m, v, :t)
                    psi(rr, zz) = puncture_psi_chart(pB, m, rr, zz; phi_p = big(phi_p)).psi
                    for i in 1:3, j in 1:3
                        r = big(r0) + big"0.2" + big(1) * (i - 1) / 2
                        z = big"0.01" + big"0.2" * (j - 1) / 2
                        d = _esd_lambda(oB, r - oB.r0, -asin(z))
                        hr = big"1e-9" * d / sqrt(oB.rho2_dr2)
                        hz = big"1e-9" * min(d / sqrt(oB.rho2_dtheta2), z / 2)
                        pr(x) = psi(x, z)
                        pz(x) = psi(r, x)
                        dpr = _esd_cd1(pr, r, hr)
                        d2pr = _esd_cd2(pr, r, hr)
                        dpz = _esd_cd1(pz, z, hz)
                        d2pz = _esd_cd2(pz, z, hz)
                        dArr = _esd_cd1(x -> fc.A[1](x, z), r, hr)
                        dAzz = _esd_cd1(x -> fc.A[2](r, x), z, hz)
                        lhs = -(fc.A[1](r, z) * d2pr + dArr * dpr) -
                              (fc.A[2](r, z) * d2pz + dAzz * dpz) +
                              fc.beta(r, z) * psi(r, z) +
                              fc.gamma[1](r, z) * dpr + fc.gamma[2](r, z) * dpz
                        S = effective_source_m(pB, m, r, z; phi_p = big(phi_p))
                        res = Float64(abs(lhs + S) / max(abs(S), abs(lhs)))
                        @test res <= 1e-10
                        gworst = max(gworst, res)
                        npt += 1
                    end
                end
            end
            @test npt == 216
            println("effective-source-direct identity (BigFloat(128), ", npt,
                    " pts): worst RELATIVE residual = ", gworst)
        end
    end

    @testset "conjugation, m = 0 reality, BigFloat(256) inheritance, allocations" begin
        for c in EFFSOURCE_CONFIGS
            p = puncs[c.idx]
            o = p.orbit
            for m in (1, 2, 5, 10),
                (dr, z) in ((0.3, 0.05), (1.0, 0.2), (-0.5, -0.1), (3.0, 0.4))

                r = o.r0 + dr
                bp = effective_source_bundle(p, m, r, z; phi_p = phi_p)
                bm = effective_source_bundle(p, -m, r, z; phi_p = phi_p)
                # S is cancellation-limited ⇒ compare against the TERM SCALE (measured 2.0e-13)
                @test abs(bm.source - conj(bp.source)) <= 1e-11 * bp.term_scale
                @test abs(bm.term_scale - bp.term_scale) <= 1e-11 * bp.term_scale
                cp = puncture_psi_chart(p, m, r, z; phi_p = phi_p)
                cm = puncture_psi_chart(p, -m, r, z; phi_p = phi_p)
                for k in (:psi, :dr, :dz, :drr, :dzz)
                    @test abs(getfield(cm, k) - conj(getfield(cp, k))) <=
                          1e-11 * abs(getfield(cp, k))
                end
            end
            # m = 0: the rotation is exactly 1 and every coefficient is real ⇒ EXACT reality
            for (dr, z) in ((0.3, 0.05), (1.0, 0.2), (-0.5, -0.1))
                @test imag(effective_source_m(p, 0, o.r0 + dr, z; phi_p = phi_p)) == 0
            end
        end
        # BigFloat(256) inheritance + the Float64 loss PREDICTED by the exposed diagnostic
        setprecision(BigFloat, 256) do
            worst = 0.0
            worst_pred = 0.0
            for c in (EFFSOURCE_CONFIGS[1], EFFSOURCE_CONFIGS[2], EFFSOURCE_CONFIGS[5])
                p = puncs[c.idx]
                o = p.orbit
                oB = CircularOrbit(big(c.a), big(c.r0))
                pB = Puncture(oB)
                for m in (0, 2, 10), (dr, z) in ((0.3, 0.05), (1.0, 0.2))
                    r = o.r0 + dr
                    b = effective_source_bundle(p, m, r, z; phi_p = phi_p)
                    sB = effective_source_m(pB, m, big(r), big(z); phi_p = big(phi_p))
                    @test sB isa Complex{BigFloat}
                    lost = Float64(abs(b.source - sB) / abs(sB))
                    # the node-37 switch criterion: relative Float64 loss ≲ eps · cancellation
                    # (main.tex:488-490); measured worst ratio 2.0e2
                    @test lost <= 1e4 * eps(Float64) * b.cancellation
                    worst = max(worst, lost)
                    worst_pred = max(worst_pred, lost / (eps(Float64) * b.cancellation))
                end
            end
            println("effective-source-direct Float64-vs-BigFloat(256) loss: worst = ", worst,
                    " ; worst loss/(eps·cancellation) = ", worst_pred)
        end
        setprecision(BigFloat, 256) do
            oB = CircularOrbit(big"0.5", big"10.0")
            pB = Puncture(oB)
            b = puncture_psi_chart(pB, 2, big"10.3", big"0.05")
            for k in (:psi, :dr, :dz, :drr, :dzz)
                @test getfield(b, k) isa Complex{BigFloat}
            end
            @test puncture_d2phi_m(pB, 2, big"0.3", big"0.05").drr isa BigFloat
            eb = effective_source_bundle(pB, 2, big"10.3", big"0.05")
            @test eb.source isa Complex{BigFloat}
            @test eb.term_scale isa BigFloat
        end
        # allocation sanity (Float64 production path; measured 0 bytes)
        let p = puncs[2]
            effective_source_m(p, 5, 10.3, 0.03; phi_p = phi_p)
            effective_source_bundle(p, 5, 10.3, 0.03; phi_p = phi_p)
            puncture_psi_chart(p, 5, 10.3, 0.03; phi_p = phi_p)
            @test (@allocated effective_source_m(p, 5, 10.3, 0.03; phi_p = phi_p)) <= 64
            @test (@allocated effective_source_bundle(p, 5, 10.3, 0.03; phi_p = phi_p)) <= 64
            @test (@allocated puncture_psi_chart(p, 5, 10.3, 0.03; phi_p = phi_p)) <= 64
        end
    end

    @testset "cancellation diagnostic: the node-37 switch data (main.tex:488-490)" begin
        # Four orders of singularity cancel in the assembly ⇒ the roundoff amplification grows
        # like λ^{-2}. The diagnostic must be monotone towards the worldline and follow that
        # law (measured ratios 5.4e3…1.8e4 over two decades of λ).
        o = CircularOrbit(0.5, 10.0)
        p = Puncture(o)
        for m in (0, 1, 2, 5, 10, 20)
            prev = 0.0
            for lg in (-1, -2, -3, -4)
                dr = exp10(float(lg)) / sqrt(o.rho2_dr2)
                b = effective_source_bundle(p, m, o.r0 + dr, 0.0; phi_p = phi_p)
                @test b.cancellation > prev
                @test b.term_scale >= abs(b.source)
                prev = b.cancellation
            end
            c1 = effective_source_bundle(p, m, o.r0 + 1e-1 / sqrt(o.rho2_dr2), 0.0;
                                         phi_p = phi_p).cancellation
            c3 = effective_source_bundle(p, m, o.r0 + 1e-3 / sqrt(o.rho2_dr2), 0.0;
                                         phi_p = phi_p).cancellation
            @test 1e3 <= c3 / c1 <= 1e5
        end
    end
end
