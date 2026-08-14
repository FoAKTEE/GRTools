# RT-A1 full reduced-ODE operator assembly via the test-function method (Modes.MasterPDE).
# GROUNDING (NON-CIRCULAR, self-contained): the master-PDE A_m (leading-order equatorial reduction
# of eq:MasterEquations, extracted by applying L to radial test functions 1,(r-r0),(r-r0)²) is
# cross-checked against the INDEPENDENT exact Abel A_m = −d/dr log𝒲 (Modes.ReducedODE.A_m, grounded
# EXACTLY by the closed-form Wronskian). A leading-order construction matches the exact value only
# STRUCTURALLY — same sign and same order of magnitude (factor ≤10), all coeffs finite. This is
# verbatim the legacy DiscModes cov_heightavg_pde `all_finite + sign_ok + max_ratio<10` assertion.
using Test
using AccretionDisc
using AccretionDisc.Modes: RelativisticDisc, master_radial_coeffs, A_m

@testset "Modes.MasterPDE (RT-A1 reduction)" begin
    d = RelativisticDisc(a = 0.6, h0 = 0.1); rp = 15.0
    all_finite = true; sign_ok = true; max_ratio = 1.0
    for (r, m) in [(10.0, 2), (12.0, 2), (20.0, 2), (11.0, 3), (25.0, 3)]
        c = master_radial_coeffs(d, r, rp, m)
        ap = real(c.A_m); aw = A_m(d, r, rp, m)
        @test isfinite(ap) && isfinite(real(c.B_m)) && isfinite(real(c.P))
        all_finite &= isfinite(ap) && isfinite(real(c.B_m))
        sign_ok &= (sign(ap) == sign(aw))
        ratio = abs(ap / aw); max_ratio = max(max_ratio, max(ratio, 1 / ratio))
    end
    @test all_finite                                 # finite A_m, B_m
    @test sign_ok                                    # A_m(PDE) sign == A_m(𝒲) sign
    @test max_ratio < 10                             # same order of magnitude (leading-order reduction)

    # BigFloat element-type inheritance: the assembly is Complex{BigFloat}.
    db = RelativisticDisc(a = big"0.6", h0 = big"0.1", r_I = big"1.0", eden_I = big"1.0")
    cb = master_radial_coeffs(db, big"12.0", big"15.0", 2)
    @test cb.A_m isa Complex{BigFloat}
    @test isfinite(real(cb.A_m)) && isfinite(real(cb.B_m))
end
