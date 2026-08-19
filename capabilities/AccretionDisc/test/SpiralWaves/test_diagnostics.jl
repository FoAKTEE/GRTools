# test_diagnostics.jl — EGD::sw-diagnostics: the AMF / torque diagnostics (69)-(76) of arXiv:2608.17003
# (src/SpiralWaves/diagnostics.jl). Gates, in order:
#   (a) √−g = e^{ν+ψ+μ} equals r for the Kerr (A1) potentials, and equals the CAS fixture column;
#   (b) THE FIG-3 GATE (M:898-908): dF_J/dr == d𝔗/dr pointwise away from corotation (|r − r_c| > 3σ_ε),
#       and the deviation converges at SECOND order (the scheme's own truncation error — the identity
#       itself is an exact on-shell symbolic result, [K71] of EGD::deriv-diagnostics);
#   (c) the (75) split T_L = T − (Δ_cT − Δ_cF_J) reproduces F_J(r_out) − F_J(r_in) — the whole chain
#       (solve → AMF → torque density → cumulative torque) has to be consistent for this to hold;
#   (d) exact scalings: T ∝ q² (m•) and T ∝ ε₀, T/F_{J,0} invariant, F_{J,0}, F_{J,1} of (76);
#   (e) T_L > 0 for the nominal disc (angular momentum carried outward by the waves);
#   (f) m_cut on a synthetic sequence;
#   (g) BigFloat path, and internal consistency (cumulative_torque ↔ torque_density, total_torque).
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using Spacetime: kerr_equatorial_potentials

include(joinpath(@__DIR__, "fixtures", "r2_diag_fixture.jl"))

_bg_ref(::Type{T} = Float64; χ = 0, aEMRI = 15, h0 = 0.05, 𝔭 = 0, 𝔮 = 0, ε0 = 1) where {T} =
    DiscBackground{T}(; χ = T(χ), aEMRI = T(aEMRI), h0 = T(h0), 𝔭 = T(𝔭), 𝔮 = T(𝔮), ε0 = T(ε0))

@testset "SpiralWaves diagnostics: AMF (72a), torque (72b), (74), (75), cutoffs (M:988)" begin

    @testset "(a) √−g = e^{ν+ψ+μ}" begin
        for χ in (-0.9, 0.0, 0.51, 0.99), r in (6.5, 15.0, 102.9)
            bg = DiscBackground(; χ = χ, aEMRI = 200.0, h0 = 0.05)
            r > r_isco(bg) || continue
            @test isapprox(sqrt_minus_g(bg, r), r; rtol = 1e-14)     # Kerr (A1): √−g ≡ r
        end
        worst = 0.0
        for pt in R2_DIAG_FIXTURE
            bg = DiscBackground(; χ = parse(Float64, pt.chi), aEMRI = parse(Float64, pt.aEMRI),
                                h0 = parse(Float64, pt.h0), 𝔭 = parse(Float64, pt.pdisk), 𝔮 = parse(Float64, pt.qdisk))
            r = parse(Float64, pt.r)
            worst = max(worst, abs(sqrt_minus_g(bg, r) - parse(Float64, pt.sqrtmg)) / abs(parse(Float64, pt.sqrtmg)))
        end
        @info "√−g vs the CAS fixture (50 points): worst rel dev = $worst"
        @test worst < 1e-14
    end

    bg = _bg_ref()
    m = 3
    g1 = build_grid(bg, m; nL = 30, nR = 30, nC = 600, n_anc = 4000)
    sol1 = solve_mode(bg, m; grid = g1)
    F1 = amf(bg, sol1)
    D1 = torque_density(bg, sol1)

    @testset "(b) FIG-3 GATE: dF_J/dr == d𝔗/dr away from corotation, 2nd-order convergent" begin
        devs = Float64[]; meds = Float64[]
        for (g, sol, F) in ((g1, sol1, F1),
                            let g = build_grid(bg, m; nL = 60, nR = 60, nC = 1200, n_anc = 4000)
                                s = solve_mode(bg, m; grid = g); (g, s, amf(bg, s))
                            end)
            rm, lhs, rhs, rel = balance_residual(bg, sol; flux = F)
            away = findall(x -> abs(x - g.r_c) > 3 * bg.σε, rm)
            @test !isempty(away)
            push!(devs, maximum(rel[away]))
            push!(meds, sort(rel[away])[cld(length(away), 2)])
            @info "Fig-3 balance, N = $(length(g.r)) nodes: max |dF_J/dr − d𝔗/dr| / max|d𝔗/dr| away from corotation (|r−r_c| > 3σ_ε = $(3*bg.σε)) = $(devs[end]), median = $(meds[end]); over the WHOLE grid (corotation included) = $(maximum(rel))"
            # the LOCAL relative form the paper's figure shows: the two curves coincide where the
            # torque density is not passing through zero
            big = [i for i in away if abs(rhs[i]) > 0.05 * maximum(abs, view(rhs, away))]
            @test maximum(i -> abs(lhs[i] - rhs[i]) / abs(rhs[i]), big) < 3e-3
        end
        @test devs[1] < 1e-4
        @test devs[2] < 3e-5
        @test 3.0 < devs[1] / devs[2] < 5.0          # second order: halving Δr quarters the deviation
        @test 3.0 < meds[1] / meds[2] < 5.0
    end

    @testset "(c) the (75) corotation split closes: T − (Δ_cT − Δ_cF_J) == F_J(r_out) − F_J(r_in)" begin
        𝔗 = cumulative_torque(bg, sol1; density = D1)
        @test total_torque(bg, sol1; density = D1) ≈ 𝔗[end] - 𝔗[1]
        clos = Float64[]
        for (g, sol, F, D) in ((g1, sol1, F1, D1),
                               let g = build_grid(bg, m; nL = 60, nR = 60, nC = 1200, n_anc = 4000)
                                   s = solve_mode(bg, m; grid = g); (g, s, amf(bg, s), torque_density(bg, s))
                               end)
            sp = corotation_split(bg, sol; flux = F, density = D)          # default window 3σ_ε
            push!(clos, abs(sp.T - sp.T_c - sp.T_L) / abs(sp.T_L))
            # the sinh corotation BLOCK (C4)-(C5) is 100× narrower than the Landau scale σ_ε and does
            # NOT contain the whole deposition — the split then stalls at ≈ 0.7 % [FINDING]
            spb = corotation_split(bg, sol; halfwidth = g.r_Cout - g.r_c, flux = F, density = D)
            @info "(75) split, N = $(length(g.r)): T = $(sp.T), T_L = $(sp.T_L), T_c = $(sp.T_c) (Δ_c𝔗 = $(sp.Δ𝔗), Δ_cF_J = $(sp.ΔF)), closure = $(clos[end]); with the (C4) block as the window: T_c = $(spb.T_c), closure = $(abs(spb.T - spb.T_c - spb.T_L) / abs(spb.T_L))"
            @test sp.window isa UnitRange
            @test isapprox(sp.T_L, F[end] - F[1]; rtol = 1e-14)
        end
        @test clos[1] < 1e-3
        @test clos[2] < 3e-4
        @test 3.0 < clos[1] / clos[2] < 5.0            # second order, as the balance law itself
    end

    @testset "(d) exact scalings and the (76) units" begin
        # T ∝ q² (m• enters the source linearly, every diagnostic is quadratic)
        s2 = solve_mode(bg, m; grid = g1, mEMRI = 2)
        @test isapprox(total_torque(bg, s2), 4 * total_torque(bg, sol1; density = D1); rtol = 1e-10)
        @test isapprox(lindblad_torque(bg, s2), 4 * lindblad_torque(bg, sol1; flux = F1); rtol = 1e-10)
        # T ∝ ε₀ (the master system, hence y, is independent of ε₀; the diagnostics carry one power)
        bg2 = _bg_ref(; ε0 = 2)
        g2 = build_grid(bg2, m; nL = 30, nR = 30, nC = 600, n_anc = 4000)
        s3 = solve_mode(bg2, m; grid = g2)
        @test g2.r == g1.r                                            # ε₀ does not move the grid
        @test s3.y ≈ sol1.y rtol = 1e-12
        @test isapprox(total_torque(bg2, s3), 2 * total_torque(bg, sol1; density = D1); rtol = 1e-10)
        # (76): F_{J,0} = ε₀ a•⁴ Ω•² q² / h₀³, F_{J,1} = F_{J,0}/a•
        u = torque_units(bg; q = 3e-5)
        @test u[1] ≈ bg.ε0 * bg.aEMRI^4 * Omega_emri(bg)^2 * (3e-5)^2 / bg.h0^3
        @test u[2] ≈ u[1] / bg.aEMRI
        @test torque_units(bg2)[1] ≈ 2 * torque_units(bg)[1]
    end

    @testset "(e) the Lindblad torque is positive for the nominal disc" begin
        TL = lindblad_torque(bg, sol1; flux = F1)
        T = total_torque(bg, sol1; density = D1)
        @info "nominal disc (χ=0, a•=15, h0=0.05, 𝔭=𝔮=0, m=3): T = $T, T_L = $TL, F_J(r_in) = $(F1[1]), F_J(r_out) = $(F1[end]), T/(F_{J,0}h₀) = $(T / (F_J0(bg) * bg.h0))"
        @test TL > 0
        @test F1[1] > 0 && F1[end] > 0            # AM flows outward across BOTH edges
        @test F1[end] > F1[1]
        @test T > 0                               # net positive torque ON THE DISC (M:1014)
        @test all(isfinite, F1) && all(isfinite, D1)
        # with the LITERAL Im[k] < 0 boundary branch (ERRATUM E-M4) the wave torque flips sign
        solp = solve_mode(bg, m; grid = g1, bc = :paper)
        TLp = lindblad_torque(bg, solp)
        @info "same mode with the printed Im[k] < 0 branch: T_L = $TLp (physical branch: $TL)"
        @test TLp < 0
    end

    @testset "(f) m_cut (M:988-1000)" begin
        Tm = [4.0, 3.0, 2.0, 1.0]                  # total 10
        @test m_cut(Tm, 0.4) == 1
        @test m_cut(Tm, 0.7) == 2
        @test m_cut(Tm, 0.9) == 3
        @test m_cut(Tm, 1.0) == 4
        @test m_cut(fill(1.0, 20), 0.45) == 9      # ⌈0.45·20⌉
        @test m_cut(fill(1.0, 20), 0.95) == 19
        @test m_cut(-Tm, 0.7) == 2                 # sign-symmetric (torque-reversal discs)
        @test_throws ArgumentError m_cut(Tm, 0.0)
        @test_throws ArgumentError m_cut([1.0, -1.0], 0.5)
    end

    @testset "(g) internal consistency, per-mode weight, and the BigFloat path" begin
        # cumulative_torque is the trapezoid of torque_density on the SAME grid
        𝔗 = cumulative_torque(bg, sol1; density = D1)
        @test 𝔗[1] == 0
        @test 𝔗[2] ≈ (g1.r[2] - g1.r[1]) * (D1[1] + D1[2]) / 2
        @test total_torque(bg, sol1; density = D1) == 𝔗[end]
        # the (73) weight of a single m ≥ 1 mode is 4π (an INTEGRAL over φ, not an average)
        hm = singular_field_mode(bg, m, g1.r[10])
        y = (sol1.y[1, 10], sol1.y[2, 10])
        @test amf(bg, sol1)[10] ≈ amf_mode(bg, m, g1.r[10], y; hmodes = hm, rc = g1.r_c)
        @test torque_density(bg, sol1)[10] ≈ torque_density_mode(bg, m, g1.r[10], y; hmodes = hm, rc = g1.r_c)
        setprecision(BigFloat, 256) do
            bgB = _bg_ref(BigFloat)
            gB = build_grid(bgB, m; nL = 8, nR = 8, nC = 60, n_anc = 400)
            solB = solve_mode(bgB, m; grid = gB, hmodes_N = 512)
            FB = amf(bgB, solB); DB = torque_density(bgB, solB)
            @test eltype(FB) === BigFloat && eltype(DB) === BigFloat
            @test all(isfinite, FB) && all(isfinite, DB)
            @test lindblad_torque(bgB, solB; flux = FB) > 0
            @info "BigFloat(256) diagnostics: T_L = $(Float64(lindblad_torque(bgB, solB; flux = FB))), T = $(Float64(total_torque(bgB, solB; density = DB)))"
        end
    end
end
