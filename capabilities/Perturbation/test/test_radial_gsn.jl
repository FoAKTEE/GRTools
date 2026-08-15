# :GSN backend (GSN-2) — grounded characterization vs the BAKED vendored-oracle fixtures
# (test/fixtures/gsn_oracle.jl: GeneralizedSasakiNakamura.jl v0.8.5 commit 982b087, solved with
# its own env at abstol=reltol=1e-12; provenance: progress/emri-disk/DG_SPIN_GSN_BENCHMARKS.md),
# PLUS the independent :MST referee, the radial_sn (s=-2 Up) cross-backend check, the paper
# Eq.-(42) identity, Wronskian constancy, the DP5-swap tolerance study, and the 3-part BigFloat
# acceptance gate (GSN_PLAN §3a/§3c/§4).
#
# Trust floor: the oracle's solved rows carry ~1e-12 ODE tolerance ⇒ solved-quantity agreement
# is gated at 1e-9 relative (headroom over the 1e-10 floor); closed-form rows at 1e-13 relative.

module TestRadialGSN

using Test
using Perturbation
import Perturbation:
    _gsn_rstar_from_r, _gsn_r_from_rstar, _gsn_rp,
    _gsn_sF, _gsn_sU, _gsn_VT, _gsn_VSN, _gsn_eta, _gsn_eta_prime, _gsn_eta_primeprime,
    _gsn_alpha, _gsn_beta, _gsn_eta_coefficient,
    _gsn_Btrans, _gsn_Binc, _gsn_Bref, _gsn_Ctrans, _gsn_Cinc, _gsn_Cref,
    _gsn_TS_abs_Csq, _gsn_TS_Bsq,
    _gsn_coeff_in_inf, _gsn_coeff_out_inf, _gsn_coeff_in_hor, _gsn_coeff_out_hor,
    _gsn_Xin_bc, _gsn_Xup_bc,
    _gsn_solve_X, _gsn_X_at, _gsn_sanity_identity, _gsn_wronskian,
    _gsn_teukolsky_radial, _gsn_value,
    SpinWeightedSpheroidalEigenvalue

include(joinpath(@__DIR__, "fixtures", "gsn_oracle.jl"))

relerr(x, y) = abs(x - y) / max(abs(y), one(abs(y)) * 1e-30)
mixerr(x, y) = abs(x - y) / max(one(abs(y)), abs(y))   # abs below 1, rel above

const RSIN_FIX, RSOUT_FIX = -20.0, 250.0   # the fixture capture domain
const GSNCFG = (rsin=RSIN_FIX, rsout=RSOUT_FIX, autowiden=false)

@testset "Radial :GSN backend (GSN-2)" begin

    # -------- 1. tortoise coordinates (33 rows) + BigFloat roundtrip ------------------------------
    @testset "gsn_coordinates vs oracle + BigFloat" begin
        for (a, r, rstar) in FIX_RSTAR
            @test relerr(_gsn_rstar_from_r(a, r), rstar) < 1e-13
        end
        for (a, rstar, r) in FIX_R_FROM_RSTAR
            @test relerr(_gsn_r_from_rstar(a, rstar), r) < 1e-12
        end
        # BigFloat roundtrip at 256 bits — the audited (0, 1.4) Float64-bracket blocker is gone
        setprecision(BigFloat, 256) do
            for (a, rs) in ((big"0.7", big"-20.0"), (big"0.7", big"35.0"), (big"0.99", big"-5.0"))
                r = _gsn_r_from_rstar(a, rs)
                @test r isa BigFloat
                @test abs(_gsn_rstar_from_r(a, r) - rs) < 1e-70 * (1 + abs(rs))
            end
        end
    end

    # -------- 2. potentials / transformation (56 rows) -------------------------------------------
    @testset "gsn_potentials + gsn_transform vs oracle" begin
        for (s, l, m, a, omega, lam, r, sf, su, et, al, be, vt) in FIX_POTENTIALS
            @test mixerr(_gsn_sF(s, m, a, omega, lam, r), sf) < 1e-13
            @test mixerr(_gsn_sU(s, m, a, omega, lam, r), su) < 1e-13
            @test mixerr(_gsn_eta(s, m, a, omega, lam, r), et) < 1e-13
            @test mixerr(_gsn_alpha(s, m, a, omega, lam, r), al) < 1e-13
            @test mixerr(_gsn_beta(s, m, a, omega, lam, r), be) < 5e-13
            @test mixerr(_gsn_VT(s, m, a, omega, lam, r), vt) < 1e-13
        end
        # internal identity: V_SN reconstructed from η stays finite and s=0 ⇒ sF ≡ 0 typed
        z = _gsn_sF(0, 2, 0.7, 0.5, 4.0, 10.0)
        @test z == 0 && z isa Complex
        @test isfinite(_gsn_VSN(-2, 2, 0.7, 0.5, 1.69660940163534, 10.0))
    end

    # -------- 3. conversion factors + Teukolsky-Starobinsky (23 rows) ----------------------------
    @testset "gsn_conversions vs oracle + TS constants" begin
        for (s, l, m, a, omega, lam, bt, ct, bi, ci) in FIX_CONVFACTORS
            @test mixerr(_gsn_Btrans(s, m, a, omega, lam), bt) < 1e-13
            @test mixerr(_gsn_Ctrans(s, m, a, omega, lam), ct) < 1e-13
            @test mixerr(_gsn_Binc(s, m, a, omega, lam), bi) < 1e-13
            @test mixerr(_gsn_Cinc(s, m, a, omega, lam), ci) < 1e-13
            # alias identities Bref ≡ Ctrans, Cref ≡ Btrans [GSN ConversionFactors.jl l.145/150]
            @test _gsn_Bref(s, m, a, omega, lam) == _gsn_Ctrans(s, m, a, omega, lam)
            @test _gsn_Cref(s, m, a, omega, lam) == _gsn_Btrans(s, m, a, omega, lam)
        end
        for (s, l, m, a, omega, lam, val) in FIX_TS_CONSTANTS
            got = abs(s) == 2 ? _gsn_TS_abs_Csq(s, m, a, omega, lam) : _gsn_TS_Bsq(s, m, a, omega, lam)
            @test mixerr(got, val) < 1e-12
        end
        # BigFloat: sqrt(2) fix — value must sharpen beyond Float64 when precision doubles
        setprecision(BigFloat, 256) do
            b256 = _gsn_Btrans(-2, 2, big"0.7", big"0.5", big"1.69660940163534")
            b64 = _gsn_Btrans(-2, 2, 0.7, 0.5, 1.69660940163534)
            @test abs(b256 - b64) < 1e-14                       # same value...
            @test !(Complex{Float64}(b256) === b256)            # ...at genuine BigFloat type
        end
    end

    # -------- 4. asymptotic-expansion coefficients (144 rows: closed-form AND recurrence regimes)
    @testset "gsn_asymptotics vs oracle (orders 0..8)" begin
        fn = Dict(1 => _gsn_coeff_in_inf, 2 => _gsn_coeff_out_inf,
                  3 => _gsn_coeff_in_hor, 4 => _gsn_coeff_out_hor)
        for (kind, s, m, a, omega, lam, order, coeff) in FIX_ASYMP_COEFFS
            # recurrence-only port must reproduce upstream's hand-coded ∞-orders 1..3 too
            # (this is the gate that licenses dropping those closed forms, GSN_PLAN §1)
            @test mixerr(fn[kind](s, m, a, omega, lam, order, Float64), coeff) < 2e-11
        end
        # cache-poisoning regression (audited blocker §3c-2): Float64 call, then BigFloat call of
        # the SAME mode must yield full-precision BigFloat coefficients (upstream downcast here)
        setprecision(BigFloat, 256) do
            c64 = _gsn_coeff_in_inf(-2, 2, 0.7, 0.5, 1.69660940163534, 6, Float64)
            cbig = _gsn_coeff_in_inf(-2, 2, big"0.7", big"0.5", big"1.69660940163534", 6, BigFloat)
            @test cbig isa Complex{BigFloat}
            @test abs(cbig - c64) < 1e-12          # same value
            # different BigFloat precisions must not alias in the cache
            c512 = setprecision(BigFloat, 512) do
                _gsn_coeff_in_inf(-2, 2, big"0.7", big"0.5", big"1.69660940163534", 6, BigFloat)
            end
            @test precision(real(c512)) == 512
        end
    end

    # -------- 5. boundary values (18 rows) --------------------------------------------------------
    @testset "gsn_boundary vs oracle" begin
        for (which, s, l, m, a, omega, lam, rsb, order, X, dX) in FIX_INITIAL_CONDITIONS
            got = which == 1 ? _gsn_Xup_bc(s, m, a, omega, lam, rsb, order, Float64) :
                               _gsn_Xin_bc(s, m, a, omega, lam, rsb, order, Float64)
            @test mixerr(got[1], X) < 1e-11
            @test mixerr(got[2], dX) < 1e-11
        end
    end

    # -------- 6. end-to-end: amplitudes (48) + Teukolsky profiles (192) + GSN profiles (192) ------
    @testset "end-to-end :GSN vs oracle net" begin
        # index fixtures by mode-config
        amps = Dict{Tuple,Tuple}()
        for (s, l, m, a, omega, bc, lam, inc, ref) in FIX_TEUK_AMPS
            amps[(s, l, m, a, omega, bc)] = (lam, inc, ref)
        end
        profs = Dict{Tuple,Vector{Tuple}}()
        for (s, l, m, a, omega, bc, r, R, dR) in FIX_TEUK_PROFILES
            push!(get!(profs, (s, l, m, a, omega, bc), Tuple[]), (r, R, dR))
        end
        gprofs = Dict{Tuple,Vector{Tuple}}()
        for (s, l, m, a, omega, bc, rs, X, dX) in FIX_GSN_PROFILES
            push!(get!(gprofs, (s, l, m, a, omega, bc), Tuple[]), (rs, X, dX))
        end
        modes = sort(collect(Set((k[1], k[2], k[3], k[4], k[5]) for k in keys(amps))))
        worst_amp = 0.0; worst_prof = 0.0; worst_gsn = 0.0; worst_sanity = 0.0
        for (s, l, m, a, omega) in modes
            lam_fix = amps[(s, l, m, a, omega, 1)][1]
            # OQ-5 λ-source gate: GrTools Harmonics vs the upstream SWSH.jl per-row value
            lam_h = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
            @test abs(real(lam_h) - lam_fix) <= 1e-11 * max(1.0, abs(lam_fix))
            R = TeukolskyRadial(s, l, m, a, omega; method=:GSN, nu=0.0 + 0.0im, gsn=GSNCFG)
            for (bci, RF) in ((1, R.In), (2, R.Up))
                lam, inc, ref = amps[(s, l, m, a, omega, bci)]
                e1 = mixerr(RF.amplitudes.Incidence, inc); e2 = mixerr(RF.amplitudes.Reflection, ref)
                worst_amp = max(worst_amp, e1, e2)
                @test e1 < 1e-9
                @test e2 < 1e-9
                @test RF.amplitudes.Transmission == 1
                for (r, Rv, dRv) in profs[(s, l, m, a, omega, bci)]
                    ep = max(mixerr(RF(r), Rv), mixerr(teukolsky_radial_deriv(RF, r), dRv))
                    worst_prof = max(worst_prof, ep)
                    @test ep < 3e-9
                end
                sol = RF.cache
                for (rs, Xv, dXv) in gprofs[(s, l, m, a, omega, bci)]
                    X, dX = _gsn_X_at(sol, rs)
                    eg = max(mixerr(X, Xv), mixerr(dX, dXv))
                    worst_gsn = max(worst_gsn, eg)
                    @test eg < 3e-9
                end
            end
            # internal gate per mode: Wronskian constancy, CONDITION-AWARE — at potential-
            # barrier modes (high l) the Wronskian is the difference of two terms that each
            # dwarf it (|Xin·Xup′| ≫ |W|), so achievable constancy is limited by the local
            # cancellation scale, not by |W| itself. Measured 2026-07-03 after the quintic-
            # Hermite fix (scratchpad/gsn-close/measure_wronskian.jl over the whole 24-mode
            # net): worst pure relative drift 1.26e-10, and the cancellation-floor term is
            # never the binding one (worst (drift−5e-9|w0|)/cancel = −1.7e-13 < 0). Gate at
            # 5e-9 relative + 1e-10·cancellation floor = ≥40x measured headroom, kept
            # condition-aware so barrier modes stay honestly gated on other platforms.
            wcond(rs) = begin
                Xi, dXi = _gsn_X_at(R.In.cache, rs); Xu, dXu = _gsn_X_at(R.Up.cache, rs)
                abs(Xi*dXu) + abs(dXi*Xu)
            end
            w0 = _gsn_wronskian(R.In.cache, R.Up.cache, -5.0)
            for rs in (0.0, 30.0, 120.0)
                tolW = 5e-9*max(1.0, abs(w0)) + 1e-10*(wcond(rs) + wcond(-5.0))
                @test abs(_gsn_wronskian(R.In.cache, R.Up.cache, rs) - w0) < tolW
            end
            # Eq.-(42) residual at the FIXED fixture orders (3/6) is BC-truncation-limited and
            # not asserted here (the identity is promised only on the auto-widen path — tested
            # in its own testset below); collect the worst relative residual as info.
            actual, expected = _gsn_sanity_identity(R.In.cache, R.Up.cache)
            worst_sanity = max(worst_sanity, Float64(abs(actual - expected)/max(1.0, abs(expected))))
        end
        @info "GSN end-to-end gate maxima" worst_amp worst_prof worst_gsn worst_sanity
    end

    # -------- 7. :MST as independent referee + radial_sn cross-backend ----------------------------
    @testset "Eq.-(42) sanity via auto-widen" begin
        # The identity Binc/Cinc = p·c0/(ω·η(r₊)) [arXiv:2306.16469 Eq. (42)] is promised where
        # the auto-widen loop enforces it (relative criterion — see _gsn_homogeneous).
        for (s, l, m, a, omega) in ((-2, 2, 2, 0.7, 0.5), (0, 2, 2, 0.99, 0.3))
            lam = real(SpinWeightedSpheroidalEigenvalue(s, l, m, a*omega))
            Xin, Xup = Perturbation._gsn_homogeneous(s, l, m, a, omega, lam;
                                                     rsin=RSIN_FIX, rsout=RSOUT_FIX, autowiden=true)
            actual, expected = _gsn_sanity_identity(Xin, Xup)
            @test abs(actual - expected) <= 1e-6*max(1, abs(expected))
        end
    end

    @testset ":MST referee + radial_sn cross-check" begin
        for (s, l, m, a, omega) in ((-2, 2, 2, 0.7, 0.5), (0, 2, 2, 0.7, 0.5), (2, 2, 2, 0.0, 0.5))
            G = TeukolskyRadial(s, l, m, a, omega; method=:GSN, gsn=GSNCFG)
            M = TeukolskyRadial(s, l, m, a, omega; method=:MST)
            @test relerr(G.In.amplitudes.Incidence, M.In.amplitudes.Incidence) < 1e-8
            @test mixerr(G.In.amplitudes.Reflection, M.In.amplitudes.Reflection) < 1e-8
            for r in (6.0, 12.0)
                @test relerr(G.In(r), M.In(r)) < 1e-8
                @test relerr(G.Up(r), M.Up(r)) < 1e-8
            end
        end
        # radial_sn (s=-2 Up only, MST-anchored) against GSN at its one supported point
        G = TeukolskyRadial(-2, 2, 2, 0.7, 0.5; method=:GSN, gsn=GSNCFG)
        S = TeukolskyRadial(-2, 2, 2, 0.7, 0.5; method=:SasakiNakamura)
        for r in (8.0, 20.0)
            @test relerr(G.Up(r), S.Up(r)) < 1e-7
        end
    end

    # -------- 8. ODE-swap validation: DP5 tolerance study (GSN_PLAN §3a) --------------------------
    @testset "DP5 tolerance study (vs :MST referee)" begin
        M = TeukolskyRadial(-2, 2, 2, 0.7, 0.5; method=:MST)
        errs = Float64[]
        for tol in (1e-6, 1e-8, 1e-10, 1e-12)
            G = _gsn_teukolsky_radial(-2, 2, 2, 0.7, 0.5, 1.69660940163534, 0.0 + 0.0im;
                                      rsin=RSIN_FIX, rsout=RSOUT_FIX, autowiden=false, tol=tol)
            push!(errs, relerr(G.In.amplitudes.Incidence, M.In.amplitudes.Incidence))
        end
        @info "GSN DP5 tolerance study |ΔBinc|/|Binc| vs MST at tol ∈ {1e-6,1e-8,1e-10,1e-12}" errs
        @test errs[end] < 5e-9                    # converged to the referee
        @test errs[end] <= errs[1] + 1e-12        # no worse than the loose-tol solve
        @test errs[2] <= 10 * errs[1] + 1e-12     # sane ordering with headroom for noise floors
        # both integration forms agree
        GR = _gsn_solve_X(-2, 2, 2, 0.7, 0.5, 1.69660940163534, :In;
                          rsin=RSIN_FIX, rsout=RSOUT_FIX, form=:riccati)
        GL = _gsn_solve_X(-2, 2, 2, 0.7, 0.5, 1.69660940163534, :In;
                          rsin=RSIN_FIX, rsout=RSOUT_FIX, form=:linear)
        @test relerr(GR.inc_SN, GL.inc_SN) < 1e-9
    end

    # -------- 9. BigFloat acceptance gate (3-part, GSN_PLAN §3c as amended by deviation D-GBS) --
    @testset "BigFloat gate" begin
        s, l, m = -2, 2, 2
        lamf = SpinWeightedSpheroidalEigenvalue(s, l, m, 0.7 * 0.5)
        # (iii-a) closed-form/coefficient/boundary layers: 256-vs-512 agreement genuinely BEYOND
        # Float64 (these layers have no integrator in the loop)
        v256 = setprecision(BigFloat, 256) do
            (_gsn_Btrans(-2, 2, big"0.7", big"0.5", big(lamf)),
             _gsn_sU(-2, 2, big"0.7", big"0.5", big(lamf), big"10.0"),
             _gsn_coeff_in_inf(-2, 2, big"0.7", big"0.5", big(lamf), 8, BigFloat),
             _gsn_Xin_bc(-2, 2, big"0.7", big"0.5", big(lamf), big"-20.0", 3, BigFloat)[1])
        end
        v512 = setprecision(BigFloat, 512) do
            (_gsn_Btrans(-2, 2, big"0.7", big"0.5", big(lamf)),
             _gsn_sU(-2, 2, big"0.7", big"0.5", big(lamf), big"10.0"),
             _gsn_coeff_in_inf(-2, 2, big"0.7", big"0.5", big(lamf), 8, BigFloat),
             _gsn_Xin_bc(-2, 2, big"0.7", big"0.5", big(lamf), big"-20.0", 3, BigFloat)[1])
        end
        for (x256, x512) in zip(v256, v512)
            @test x256 isa Complex{BigFloat}
            @test abs(x256 - x512) < big"1e-60" * max(1, abs(x512))   # far beyond Float64 eps
        end
        # (i)+(ii)+(iii-b) solve tier: the type-generic DP5 mesh at the BigFloat default
        # tol 1e-12 (deviation D-GBS: the GBS phase tier has a measured accuracy-biasing
        # extrapolation floor; see gsn_solve.jl). Threading + cross-precision determinism.
        Gf = _gsn_teukolsky_radial(s, l, m, 0.7, 0.5, lamf, 0.0 + 0.0im;
                                   rsin=-10.0, rsout=60.0, autowiden=false)
        inc256 = setprecision(BigFloat, 256) do
            G = _gsn_teukolsky_radial(s, l, m, big"0.7", big"0.5", big(lamf),
                                      Complex{BigFloat}(0);
                                      rsin=big"-10.0", rsout=big"60.0", autowiden=false)
            G.In.amplitudes.Incidence
        end
        inc512 = setprecision(BigFloat, 512) do
            G = _gsn_teukolsky_radial(s, l, m, big"0.7", big"0.5", big(lamf),
                                      Complex{BigFloat}(0);
                                      rsin=big"-10.0", rsout=big"60.0", autowiden=false)
            G.In.amplitudes.Incidence
        end
        @test inc256 isa Complex{BigFloat}                                  # (i) threading
        @test abs(inc256 - Gf.In.amplitudes.Incidence) < 1e-8 * abs(inc256) # (ii) Float64 agreement
        # (iii-b) cross-precision determinism: identical DP5 algorithm/truncation at
        # 256 vs 512 bits — the difference reflects only round-off propagation
        @test abs(inc256 - inc512) < 1e-10 * abs(inc256)
    end

    # -------- 10. quarantines + error surface ------------------------------------------------------
    @testset "quarantines (OQ-4/D-dispositions)" begin
        # complex ω: explicitly quarantined to GSN-3 (ruling OQ-4)
        @test_throws ArgumentError _gsn_solve_X(-2, 2, 2, 0.7, 0.5 + 0.1im, 1.7, :In;
                                                rsin=-20.0, rsout=250.0)
        # only :In/:Up on the solve surface (Out/Down internal-only, ruling OQ-2)
        @test_throws ArgumentError _gsn_solve_X(-2, 2, 2, 0.7, 0.5, 1.7, :Down;
                                                rsin=-20.0, rsout=250.0)
        # dispatcher surface: unknown method still errors, :GSN accepted
        @test_throws ErrorException TeukolskyRadial(-2, 2, 2, 0.7, 0.5; method=:NoSuch)
    end
end

end # module
