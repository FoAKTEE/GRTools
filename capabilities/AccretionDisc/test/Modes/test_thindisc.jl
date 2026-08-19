# Module — Modes/ThinDisc.jl + the Forcing Newtonian-amplitude leaves. PROMOTED from the legacy
# figure-local library `results/emri-disk/figures/2510.02433/lindblad_torque_lib.jl` (arXiv:2510.02433
# RT5/RT6 + tab:FitParameters) and `repro_discProfile.jl` (RT4), which had NO package home
# (`grep -i 'thorne|shakura' GrTools/capabilities/*/src` = 0 hits before this node).
#
# GROUNDING (every number below is either a CLOSED-FORM identity, a PUBLISHED anchor, or a ONE-SHOT
# capture of the legacy oracle):
#   • legacy capture: `julia --project=results/emri-disk/figures <scratch>/capture_legacy_thindisc.jl`
#     includes lindblad_torque_lib.jl VERBATIM and prints the literals baked below (17 sig digits).
#   • published anchors: Hirata 1010.0759 RT-N1 Newtonian-limit strengths N = −2.36 (m=2 ILR),
#     −7.50 (m=3 ILR), +0.19 (m=1 OLR), +8.62 (m=2 OLR); arXiv:2510.02433 h(a=0, r=1000M) ≈ 0.0033.
#   • closed-form identities: Newtonian limit A,B,C,D,E→1, Q→y ⇒ Σ_NT→Σ_SS ∝ r^{3/2}, h_NT→h_SS ∝ 1/r;
#     the specific epicyclic impedance Z = κ u^t r²/Δ (Modes._kZ) equals the legacy metric form
#     √(C₀₂₀/C₀₀₂) — REUSED, not re-implemented (dedup: one home).
using Test
using AccretionDisc
using AccretionDisc.Modes
const MT = AccretionDisc.Modes

@testset "Modes.ThinDisc" begin

    # ── legacy one-shot captures (see header) ────────────────────────────────────────────────────
    legacy_factors = Dict(
        (0.0, 10.0)   => (A=1.0, B=1.0, C=0.69999999999999996, D=0.80000000000000004, E=1.0, Q=0.09510593383579144),
        (0.9, 10.0)   => (A=1.00972, B=1.0284604989415154, C=0.7569209978830308, D=0.80810000000000004, E=1.02935683, Q=0.36012207527247098),
        (0.0, 1000.0) => (A=1.0, B=1.0, C=0.997, D=0.998, E=1.0, Q=0.87858738994860397),
        (0.9, 1000.0) => (A=1.0000008116199999, B=1.0000284604989416, C=0.99705692099788301, D=0.99800080999999996, E=1.0000032367619682, Q=0.91059550362223596),
        (0.3, 25.0)   => (A=1.0001555199999999, B=1.0024, C=0.88480000000000003, D=0.92014399999999996, E=1.0005530222079999, Q=0.39825211390279164))
    legacy_sigma = Dict((0.0,10.0)=>139094.96518023298, (0.9,10.0)=>41953.907661489502,
                        (0.0,100.0)=>766558.25792169757, (0.9,100.0)=>665676.75169179239,
                        (0.0,1000.0)=>17969360.344188981, (0.9,1000.0)=>17339727.176367544)
    legacy_h     = Dict((0.0,10.0)=>0.037340367124365545, (0.9,10.0)=>0.13252413222287035,
                        (0.0,100.0)=>0.02423719922498439, (0.9,100.0)=>0.027959762344051727,
                        (0.0,1000.0)=>0.0033, (0.9,1000.0)=>0.0034200207367766377)
    legacy_sigma_ss = Dict(10.0=>17969.360344188983, 100.0=>568241.06783944403, 1000.0=>17969360.344188981)
    legacy_h_ss     = Dict(10.0=>0.32999999999999996, 100.0=>0.033000000000000002, 1000.0=>0.0033)
    legacy_pcut = Dict(0.0=>1.0, 0.5=>0.37754066879814541, 1.0=>0.10757656854799805, 3.0=>0.0025635607123009074)
    legacy_b    = Dict((1,0.5)=>0.55586619792668102, (2,0.63)=>0.36537190748965526,
                       (3,0.7)=>0.28512795003441788, (1,0.1)=>0.10037736097554596)
    legacy_db   = Dict((1,0.5)=>1.3795088245938218, (2,0.63)=>1.4601790493166851,
                       (3,0.7)=>1.541843434153916, (1,0.1)=>1.0113683960248638)
    legacy_N    = Dict((1,0.0,50.0,:ILR)=>-0.23149152841386023, (2,0.0,50.0,:ILR)=>-3.413593422964512,
                       (2,0.0,50.0,:OLR)=>9.8598854793631894,  (1,0.0,50.0,:OLR)=>0.22705867015167266,
                       (3,0.9,20.0,:ILR)=>-13.594113058346728, (2,0.9,20.0,:OLR)=>11.694980673184427)
    legacy_Nlim = Dict((2,:ILR)=>-2.3564658896560582, (3,:ILR)=>-7.5046537302809746,
                       (1,:OLR)=>0.19219467525449724, (2,:OLR)=>8.6149743714028251)
    legacy_TLR_NT = Dict((0.0,10.0)=>-14891722512.296196, (0.0,30.0)=>30665905549.235046,
                         (0.9,10.0)=>385456850.98508561, (0.9,30.0)=>12865201476.307467)
    legacy_TLR_SS = Dict(10.0=>37306276.302393027, 30.0=>4357667476.6736422)

    @testset "laplace_b_deriv — analytic derivative vs FD of laplace_b + legacy capture" begin
        for ((m, β), val) in legacy_db
            @test laplace_b_deriv(1//2, m, β) ≈ val rtol = 1e-9
            @test laplace_b(1//2, m, β) ≈ legacy_b[(m, β)] rtol = 1e-11
        end
        # independent oracle: central finite difference of laplace_b
        h = 1e-5
        for (m, β) in ((1, 0.5), (2, 0.63), (3, 0.7))
            fd = (laplace_b(1//2, m, β + h) - laplace_b(1//2, m, β - h)) / (2h)
            @test laplace_b_deriv(1//2, m, β) ≈ fd rtol = 1e-8
        end
        @test laplace_b_deriv(1//2, 1, 1e-6) ≈ 1.0 rtol = 1e-9   # b^{(1)′}(0) = 1 (Hirata l.916)
        @test laplace_b_deriv(1//2, 1, big"0.5") isa BigFloat     # BigFloat inheritance
    end

    @testset "epicyclic_impedance — Schwarzschild closed form + legacy metric form" begin
        d0 = RelativisticDisc(a = 0.0)
        for r in (8.0, 14.5, 33.444855735549922)
            @test epicyclic_impedance(d0, r) ≈ (1/(r-2)) * sqrt((r-6)/(r*(r-3))) rtol = 1e-13
        end
        # legacy Zimp = √(C₀₂₀/C₀₀₂) at the m=2 ILR of r_p=50 (a=0) and (3, a=0.9, r_p=20)
        @test epicyclic_impedance(d0, 33.444855735549922) ≈ 0.0052210667240685455 rtol = 1e-7
        @test epicyclic_impedance(RelativisticDisc(a = 0.9), 15.959686605566404) ≈ 0.016623713545128146 rtol = 1e-7
    end

    @testset "lindblad_detuning / _deriv — root identity + legacy D′" begin
        d0 = RelativisticDisc(a = 0.0); rp = 50.0
        rI, rO = lindblad_locations(d0, rp, 2)
        @test rI ≈ 33.444855735549922 rtol = 1e-9
        @test rO ≈ 64.826565242937761 rtol = 1e-9
        @test abs(lindblad_detuning(d0, rI, rp, 2; kind = :ILR)) < 1e-11
        @test abs(lindblad_detuning(d0, rO, rp, 2; kind = :OLR)) < 1e-11
        @test lindblad_detuning_deriv(d0, rI, 2; kind = :ILR) ≈ -0.00026901705971367484 rtol = 1e-6
        @test lindblad_detuning_deriv(d0, rO, 2; kind = :OLR) ≈ -0.00012945630078626166 rtol = 1e-6
    end

    @testset "lindblad_strength — hybrid kernel vs legacy N_signed" begin
        for ((m, a, rp, kind), val) in legacy_N
            @test lindblad_strength(RelativisticDisc(a = a), rp, m; kind = kind) ≈ val rtol = 1e-6
        end
        # sign convention: ILR < 0 (inward), OLR > 0 (outward)
        @test lindblad_strength(RelativisticDisc(a = 0.0), 50.0, 2; kind = :ILR) < 0
        @test lindblad_strength(RelativisticDisc(a = 0.0), 50.0, 2; kind = :OLR) > 0
    end

    @testset "lindblad_strength_newtonian — PUBLISHED Hirata RT-N1 constants + r_p-independence" begin
        # published 1010.0759 values (l.1008 + eq:smolr with the corrected ς⁺ argument)
        @test lindblad_strength_newtonian(2; kind = :ILR) ≈ -2.36 rtol = 3e-3
        @test lindblad_strength_newtonian(3; kind = :ILR) ≈ -7.50 rtol = 3e-3
        @test lindblad_strength_newtonian(1; kind = :OLR) ≈ 0.19 rtol = 2e-2
        @test lindblad_strength_newtonian(2; kind = :OLR) ≈ 8.62 rtol = 3e-3
        # legacy capture = the HYBRID kernel evaluated at r_p = 1e6 (its numerical Newtonian limit).
        # It agrees with the analytic closed form only to ~2e-4: the legacy D′ uses a FIXED absolute
        # step h=1e-6 which loses ~4 digits at r₁ ≈ 6e5 (this promotion uses h = 1e-6·max(1,r)).
        for ((m, kind), val) in legacy_Nlim
            @test lindblad_strength_newtonian(m; kind = kind) ≈ val rtol = 3e-4
            @test lindblad_strength(RelativisticDisc(a = 0.0), 1.0e6, m; kind = kind) ≈ val rtol = 3e-4
        end
        # the hybrid kernel converges onto the closed form as r_p → ∞
        for (m, kind) in ((2, :ILR), (3, :ILR), (1, :OLR), (2, :OLR))
            @test lindblad_strength(RelativisticDisc(a = 0.0), 1.0e6, m; kind = kind) ≈
                  lindblad_strength_newtonian(m; kind = kind) rtol = 1e-4
        end
        @test lindblad_strength_newtonian(2; kind = :ILR) isa Float64
    end

    @testset "Page–Thorne factors + NT/SS profiles vs legacy capture" begin
        for ((a, r), f) in legacy_factors
            g = nt_factors(a, r)
            for k in (:A, :B, :C, :D, :E, :Q)
                @test getfield(g, k) ≈ getfield(f, k) rtol = 1e-12
            end
        end
        nt0 = NovikovThorneDisc(a = 0.0); nt9 = NovikovThorneDisc(a = 0.9)
        for ((a, r), val) in legacy_sigma
            nt = NovikovThorneDisc(a = a)
            @test nt_surface_density(nt, r) ≈ val rtol = 1e-11
            @test nt_aspect_ratio(nt, r) ≈ legacy_h[(a, r)] rtol = 1e-11
        end
        for (r, val) in legacy_sigma_ss
            @test ss_surface_density(nt0, r) ≈ val rtol = 1e-11
            @test ss_aspect_ratio(nt0, r) ≈ legacy_h_ss[r] rtol = 1e-11
        end
        # PUBLISHED anchor: h(a=0, 1000M) = 0.0033
        @test nt_aspect_ratio(nt0, 1000.0) ≈ 0.0033 rtol = 1e-12
        # Newtonian limit: 𝒜,ℬ,𝒞,𝒟,ℰ→1 and 𝒬→y ⇒ the SS SCALINGS Σ ∝ r^{3/2}, h ∝ 1/r. (The SS
        # overlay is anchored AT r_norm=1000, where 𝒬 is still 12% short of its limit, so the
        # RATIO Σ_NT/Σ_SS tends to a constant ≠ 1 — the log-slope is the invariant statement.)
        logslope(f, r) = (log(f(r*1.001)) - log(f(r/1.001))) / (log(r*1.001) - log(r/1.001))
        for nt in (nt0, nt9)
            @test logslope(r -> nt_surface_density(nt, r), 1e8) ≈ 1.5 rtol = 1e-3
            @test logslope(r -> nt_aspect_ratio(nt, r), 1e8) ≈ -1.0 rtol = 1e-3
        end
        @test logslope(r -> ss_surface_density(nt0, r), 1e3) ≈ 1.5 rtol = 1e-12
        @test logslope(r -> ss_aspect_ratio(nt0, r), 1e3) ≈ -1.0 rtol = 1e-12
        # the overlay coincides with the a=0 NT curve exactly at the anchor radius
        @test ss_surface_density(nt0, nt0.r_norm) ≈ nt_surface_density(nt0, nt0.r_norm) rtol = 1e-14
        @test ss_aspect_ratio(nt0, nt0.r_norm) ≈ nt_aspect_ratio(nt0, nt0.r_norm) rtol = 1e-14
        # ISCO delegation (no second kerr_isco)
        @test nt_isco(nt9) == MT.kerr_isco(0.9)
    end

    @testset "pressure_cutoff — eq:TPressureCutoff" begin
        for (x, val) in legacy_pcut
            @test pressure_cutoff(x) ≈ val rtol = 1e-13
        end
        @test pressure_cutoff(0.0) == 1.0
        @test pressure_cutoff(10.0) < 1e-4          # kills m ≫ 1/h
        @test pressure_cutoff(big"1.0") isa BigFloat
    end

    @testset "total_lindblad_torque — signed assembly vs legacy TLR_NT / TLR_SS" begin
        for ((a, rs), val) in legacy_TLR_NT
            @test total_lindblad_torque(NovikovThorneDisc(a = a), rs) ≈ val rtol = 1e-6
        end
        for (rs, val) in legacy_TLR_SS
            @test total_lindblad_torque(NovikovThorneDisc(a = 0.0), rs; profile = :ss) ≈ val rtol = 1e-6
        end
        # sign structure: ILR-dominated (negative) close in, OLR-dominated (positive) further out
        @test total_lindblad_torque(NovikovThorneDisc(a = 0.0), 10.0) < 0
        @test total_lindblad_torque(NovikovThorneDisc(a = 0.0), 30.0) > 0
        # an enhancement callable multiplies each mode
        base = total_lindblad_torque(NovikovThorneDisc(a = 0.0), 30.0)
        @test total_lindblad_torque(NovikovThorneDisc(a = 0.0), 30.0;
                                    enhancement = (m, r, kind) -> 2.0) ≈ 2 * base rtol = 1e-12
    end

    @testset "powerlaw_fit — exact on a synthetic power law" begin
        rs = collect(range(10.0, 100.0; length = 25))
        f = powerlaw_fit(rs, 2.5 .* (rs ./ 10) .^ 3.4)
        @test f.A ≈ 2.5 rtol = 1e-10
        @test f.n ≈ 3.4 rtol = 1e-10
        @test f.A ≈ 2.5000000000000493 rtol = 1e-12   # legacy fit_single capture
        @test f.n ≈ 3.3999999999999879 rtol = 1e-12
    end

    @testset "c13_alpha — Hirata amplitude bundle, eq:cplex1 closure" begin
        # |c₁₃|² + α = (2ω)⁶/|C|² |ℵ|²  (1010.0759 eq:cplex1) with the Teukolsky–Press ε=Γ/2 factor
        for (l, m, a, w) in ((2, 2, 0.5, 0.12), (2, 1, 0.5, 0.08), (3, 2, 0.7, 0.15),
                             (2, 2, 0.0, 0.10), (4, 3, 0.9, 0.20), (3, 3, 0.3, 0.18))
            cc = c13_alpha(l, m, a, w; Nm = 300)
            rhs = (2w)^6 / cc.Csq * abs2(cc.aleph)
            @test abs2(cc.c13) + cc.alpha ≈ rhs rtol = 1e-11
        end
        # the (c13, alpha) pair is the one lindblad_torque consumes (captured DiscModes oracle)
        cc = c13_alpha(2, 2, 0.0, 2 * MT._kOmega(0.0, 20.0); Nm = 300)
        @test cc.c13 ≈ (-23.58767132477094 - 100.88873420202292im) rtol = 1e-8
        # NOTE the horizon factor α is CORRECTED relative to the captured DiscModes literal
        # (2.8512019997500306e-5): that value carries Γ = 4√(1−a²)/(r₊²+a²), i.e. Hirata's
        # AS-PRINTED bracket, which fails the eq:cplex1 closure by up to 9.1e-5. With the
        # Teukolsky–Press Γ = 2√(1−a²)/(r₊²+a²) the closure holds to ≤6e-13 (gate above) and
        # α is 15.88× smaller. lindblad_torque shifts by only 2.0e-8 (α ≪ |c₁₃|²).
        @test cc.alpha ≈ 1.7953595800524933e-6 rtol = 1e-8
    end

end
