# Grounding for the :SasakiNakamura radial backend (radial_sn.jl, s=-2 Up only). CHARACTERIZED vs the
# LIGHT live oracle `Teukolsky` + cross-checked against the already-green :MST Up:
#   (1) _sn_outer_c 13-term OuterBCs recurrence == legacy Teukolsky._sn_outer_c at 1e-12 (identical
#       analytic recurrence — a pure algebraic characterization);
#   (2) public TeukolskyRadial(-2,l,m,a,ω;method=:SasakiNakamura).Up VALUE + DERIV == legacy SN at the
#       legacy SN rtol (both integrate the SN ODE inward; GrTool uses the in-firewall NumericsCore DP5);
#   (3) STRONG cross-check :SN Up ≈ already-green :MST Up (same physical solution after scale-matching)
#       + an r-independent Wronskian W(In_MST, Up_SN);
#   (4) NotYetPorted throw for s≠-2 and for the In BC, matching the legacy faithful restriction.
# Replaces the :SasakiNakamura @test_throws that used to live in test_radial_carrier.jl.

include("fixtures/leg_radial_sn.jl")  # REF_SN_LEG / REF_SN_RESTRICT (baked Teukolsky SN oracle + inputs, Phase 7-B)
using Perturbation: TeukolskyRadial, teukolsky_radial_deriv

# s=-2 Up modes (the legacy REF_SN grid). λ,ν pinned via the legacy oracle so GrTool and legacy are
# apples-to-apples (identical inputs → the SN coefficients + outer BC are the same port).
const _SN_MODES = [
    (a=0.0, w=0.1), (a=0.5, w=0.1), (a=0.5, w=0.3), (a=0.9, w=0.1),
]

@testset "Radial: _sn_outer_c OuterBCs recurrence == legacy Teukolsky._sn_outer_c (1e-12, baked, captured once)" begin
    for md in _SN_MODES
        leg = REF_SN_LEG[(md.a, md.w)]
        lam = leg.lam
        cg = Perturbation._sn_outer_c(lam, 1, md.a, md.w)
        cl = leg.snc
        @test length(cg) == length(cl)
        for (a, b) in zip(cg, cl)
            @test isapprox(a, b; rtol=1e-12, atol=1e-14)
        end
    end
end

@testset "Radial: public :SasakiNakamura Up value+deriv == legacy Teukolsky :SN (baked, captured once)" begin
    for md in _SN_MODES
        leg = REF_SN_LEG[(md.a, md.w)]
        lam = leg.lam
        nu  = leg.nu
        Rg = TeukolskyRadial(-2, 2, 1, md.a, md.w; method=:SasakiNakamura, lambda=lam, nu=nu)
        @test Rg.Up.method === :SasakiNakamura
        @test Rg.Up.amplitudes.Transmission ≈ 1
        for (j, rr) in enumerate((2.5, 3.0, 5.0, 10.0, 30.0))
            # both integrate the SN equation inward from rout=1000 (GrTool via NumericsCore DP5,
            # legacy via its inline DP5) then scale-match to their own MST Up — same physical solution.
            @test isapprox(Rg.Up(rr), leg.up[j]; rtol=1e-6, atol=1e-8)
            @test isapprox(teukolsky_radial_deriv(Rg.Up, rr), leg.upd[j]; rtol=1e-6, atol=1e-8)
        end
    end
end

@testset "Radial: STRONG cross-check :SN Up ≡ :MST Up + r-independent Wronskian" begin
    for md in _SN_MODES
        leg = REF_SN_LEG[(md.a, md.w)]
        lam = leg.lam
        nu  = leg.nu
        Rs = TeukolskyRadial(-2, 2, 1, md.a, md.w; method=:SasakiNakamura, lambda=lam, nu=nu)
        Rm = TeukolskyRadial(-2, 2, 1, md.a, md.w; method=:MST, lambda=lam, nu=nu)
        for rr in (2.5, 3.0, 5.0, 10.0, 30.0)
            (md.a == 0.0 || rr > 1 + sqrt(1 - md.a^2)) || continue
            rtolUp = rr ≤ 2.5 ? 1e-5 : 1e-6                 # r=2.5 limited by MST near-horizon Up
            @test isapprox(Rs.Up(rr), Rm.Up(rr); rtol=rtolUp, atol=1e-8)
        end
        # Wronskian W = In_MST·Up_SN' - In_MST'·Up_SN of the homogeneous s=-2 pair: Δ^{s+1}·W = Δ^{-1}·W
        # is r-independent. In from MST (regular at horizon), Up from SN (regular at infinity).
        wron(r) = begin
            Δ = r^2 - 2r + md.a^2
            Δ^(-2 + 1) * (Rm.In(r) * teukolsky_radial_deriv(Rs.Up, r) -
                          teukolsky_radial_deriv(Rm.In, r) * Rs.Up(r))
        end
        @test isapprox(wron(6.0), wron(12.0); rtol=1e-5)
    end
end

@testset "Radial: :SasakiNakamura faithful restriction (s≠-2 / In throw NotYetPorted)" begin
    lam0 = REF_SN_RESTRICT.lam0
    nu0  = REF_SN_RESTRICT.nu0
    Rwrong_s = TeukolskyRadial(0, 2, 1, 0.5, 0.2; method=:SasakiNakamura, lambda=lam0, nu=nu0)
    @test_throws Perturbation.NotYetPorted Rwrong_s.Up(5.0)     # s≠-2 → NotYetPorted
    lam2 = REF_SN_RESTRICT.lam2
    nu2  = REF_SN_RESTRICT.nu2
    Rs2 = TeukolskyRadial(-2, 2, 1, 0.5, 0.2; method=:SasakiNakamura, lambda=lam2, nu=nu2)
    @test_throws Perturbation.NotYetPorted Rs2.In(5.0)          # In BC → NotYetPorted
    @test isfinite(Rs2.Up(5.0))                                 # s=-2 Up works
end
