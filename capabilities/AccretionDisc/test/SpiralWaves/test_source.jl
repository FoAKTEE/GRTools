# SpiralWaves/test_source.jl — EGD::sw-singular-source (arXiv:2608.17003 Sec. II.D, main.tex:653-711,
# appendix.tex:27-53 eqs (52)-(60), (B1)-(B2)). Tests FIRST for `source.jl`: the leading-order singular
# field h_μν = (2m•/s0)(g_μν(w) + 2v_μv_ν) (52) of the circular equatorial EMRI worldline (53)-(55), the
# s0² approximation (58) in its printed equatorial closed form (B1), the Plummer softening (59), the
# cosine Fourier projection (60) by the periodic trapezoid rule, the analytic ∂_r of every mode, both
# index positions (covariant = the master-system convention, ERRATUM E-M2; contravariant via the background
# metric at the field point), and the grid-vectorised `source_modes`. GROUNDING:
#   (i)   s0² (B1) vs a DIRECT evaluation of (58) with (56)-(57): Jacobian ∂x/∂y at the worldline point,
#         [g + v v]_μν from `Spacetime.kerr_metric` (the dense BL 4×4 — an independent metric route),
#         Δy from (57) at RANDOM (t, φ) — 1e-13 relative at 20 random points;
#   (ii)  Newtonian limit a• ≫ M (χ = 0): ĥ_tt,m/(−2) → the softened Plummer potential mode
#         Φ̂_m(r) = −(1/2π)∫cos(mΨ)dΨ/√(r² + a•² − 2ra• cosΨ + ε²), evaluated INDEPENDENTLY through the
#         existing `Modes.laplace_b` (the exact substitution β + 1/β = (r² + a•² + ε²)/(ra•) maps the
#         softened kernel onto the Laplace coefficient b_{1/2}^{(m)}(β)): agreement O(M/a•), ratio ∝ 1/a•
#         between a• = 1e4 and 1e6; sign convention Φ̂_m < 0 (attractive), ĥ_tt > 0;
#   (iii) m-decay: |ĥ_m| strictly decreasing for m beyond a•/ε_soft (softened kernel);
#   (iv)  trapezoid N-convergence: doubling N changes every mode by < 1e-13 RELATIVE for b ≥ 0.3
#         (the strip half-width of the integrand is δ ≈ ε_soft/a• = h0 b at r = a•; the default node
#         count is derived from the exact nearest complex zero of s0² — b → 0 would need N → ∞);
#   (v)   BigFloat(256) path on a tiny case, agreeing with Float64 to 1e-13 and with itself at doubled N;
#   (vi)  analytic ∂_r of every mode vs central finite differences (Richardson) to 1e-8 relative;
#   (vii) worldline: v·v = −1, v_t = −E, v_φ = L (Bardeen–Press–Teukolsky), dt/dτ = (55);
#   (viii) nonzero pattern of (52) and the per-component constant tensor K_μν = 2(g_μν(w) + 2v_μv_ν);
#   (ix)  `source_modes` == pointwise `singular_field_mode` bitwise, threaded == serial bitwise;
#   (x)   contravariant modes: ĥ^{μν} = g^{μα}g^{νβ}ĥ_αβ vs a dense `Spacetime.kerr_metric_inv`
#         contraction and FD of the r-derivatives;
#   (xi)  end-to-end: the relativistic master source S (master.jl, route R1) fed with the singular-field
#         modes agrees with the Newtonian S (50f-g) with Φ̂ = −ĥ_tt/2 to O(M/a•) at a• = 1e4.
using Test
using LinearAlgebra
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.Modes: laplace_b
using Spacetime: kerr_metric, kerr_metric_inv, kerr_circular_dtdtau, kerr_equatorial_potentials
using Spacetime: kerr_equatorial_potentials_dr
import Spacetime

const _SRC_KEYS = (:htt, :htt_r, :htphi, :htphi_r, :hphiphi, :hphiphi_r, :hrr, :hrr_r, :hzz, :hzz_r)

# --- independent route for (58): y^μ = (t, r cosφ, r sinφ, z) (56), Jacobian ∂x^μ/∂y^α at w (53),
#     [g_μν + v_μ v_ν](w) from the dense BL Kerr metric, Δy^α (57) ----------------------------------
function _s0sq_direct(χ, a, h0, b, r, φ, t)
    ΩE = 1 / (a^(3 / 2) + χ)
    τ = (a^(3 / 2) + χ) / (a^(3 / 4) * sqrt(sqrt(a) * (a - 3) + 2χ))          # (55) literally
    gBL = kerr_metric(χ, a, π / 2)                                            # order (t, r, θ, φ)
    perm = [1, 2, 4, 3]                                                        # → (t, r, φ, z≡θ-dir)
    g = gBL[perm, perm]
    g[4, 4] = 1.0                                                              # dz² on the equator (z = r cosθ → dz = −r dθ... g_zz = 1 by (4))
    vup = τ .* [1.0, 0.0, ΩE, 0.0]                                             # (54)
    vdn = g * vup
    P = g + vdn * vdn'
    φw = ΩE * t
    # x(y) = (Y0, √(Y1²+Y2²), atan(Y2, Y1), Y3); ∂x/∂y at (t, a cosφw, a sinφw, 0)
    J = [1.0 0.0 0.0 0.0;
         0.0 cos(φw) sin(φw) 0.0;
         0.0 -sin(φw)/a cos(φw)/a 0.0;
         0.0 0.0 0.0 1.0]
    Δy = [0.0, r * cos(φ) - a * cos(φw), r * sin(φ) - a * sin(φw), 0.0]
    ε = a * h0 * b
    dot(Δy, (J' * P * J) * Δy) + ε^2
end

# softened Plummer potential mode via the Laplace coefficient (independent route; Modes.laplace_b)
function _plummer_mode(m, r, a, ε)
    R2 = r^2 + a^2 + ε^2
    x = R2 / (r * a)                                   # β + 1/β = x > 2
    β = (x - sqrt(x^2 - 4)) / 2
    -sqrt(β / (r * a)) * laplace_b(0.5, m, β; N = 4096) / 2
end

# Richardson-extrapolated central difference of a scalar function
function _fd(f, r; h = 1e-4 * r)
    d1 = (f(r + h) - f(r - h)) / (2h)
    d2 = (f(r + h / 2) - f(r - h / 2)) / h
    (4d2 - d1) / 3
end

@testset "SpiralWaves source: singular field (52)-(60), (B1)" begin
    SW = AccretionDisc.SpiralWaves
    bg = DiscBackground(; χ = 0.9, aEMRI = 8.0, h0 = 0.1, b = 0.6)

    @testset "(vii) worldline (53)-(55): normalisation, E, L, dt/dτ" begin
        for (χ, a) in ((0.9, 8.0), (0.0, 10.0), (-0.7, 12.0), (0.99, 3.0))
            bgw = DiscBackground(; χ = χ, aEMRI = a, h0 = 0.1)
            w = worldline_velocity(bgw)
            τ = (a^(3 / 2) + χ) / (a^(3 / 4) * sqrt(sqrt(a) * (a - 3) + 2χ))
            @test w.ut ≈ τ rtol = 1e-14
            @test w.ut == kerr_circular_dtdtau(χ, a)
            @test w.uphi ≈ Omega_emri(bgw) * w.ut rtol = 1e-14
            @test w.ut * w.vt + w.uphi * w.vphi ≈ -1 rtol = 1e-13            # v^μ v_μ = −1
            E = (a^(3 / 2) - 2 * sqrt(a) + χ) / (a^(3 / 4) * sqrt(a^(3 / 2) - 3 * sqrt(a) + 2χ))
            L = (a^2 - 2χ * sqrt(a) + χ^2) / (a^(3 / 4) * sqrt(a^(3 / 2) - 3 * sqrt(a) + 2χ))
            @test w.vt ≈ -E rtol = 1e-13
            @test w.vphi ≈ L rtol = 1e-13
        end
    end

    @testset "(i) s0² (B1) == direct (58) with (56)-(57), 20 random points" begin
        seed = 20260819                      # deterministic pseudo-random points via a simple LCG
        lcg() = (seed = (1103515245 * seed + 12345) % 2^31; seed / 2^31)
        maxrel = 0.0
        for _ in 1:20
            a = 4 + 296 * lcg(); χ = -0.9 + 1.8 * lcg(); r = a * (0.7 + 0.6 * lcg())
            φ = 2π * lcg(); t = 100 * lcg(); h0 = 0.02 + 0.2 * lcg(); b = 0.1 + 0.6 * lcg()
            a > 1.01 * Spacetime.kerr_isco(χ) || continue
            bgi = DiscBackground(; χ = χ, aEMRI = a, h0 = h0, b = b)
            ΩE = Omega_emri(bgi)
            Ψ = φ - ΩE * t
            s_b1 = s0_squared(bgi, r, Ψ)
            s_pr = s0_squared(bgi, r, Ψ; form = :printed)              # the literal regrouped (B1)
            s_dir = _s0sq_direct(χ, a, h0, b, r, φ, t)
            maxrel = max(maxrel, abs(s_b1 - s_dir) / abs(s_dir), abs(s_pr - s_dir) / abs(s_dir))
            @test s_b1 ≈ s_dir rtol = 1e-13
            @test s_pr ≈ s_dir rtol = 1e-13
            @test s0_squared_dr(bgi, r, Ψ) ≈ s0_squared_dr(bgi, r, Ψ; form = :printed) rtol = 1e-10
            @test s_b1 > 0
        end
        println("  (i) s0²: max relative |(B1) − direct (58)| over the random points = ", maxrel)
        # evenness (B2): s0²(−Ψ) = s0²(Ψ); 2π-periodicity
        @test s0_squared(bg, 9.0, 0.7) == s0_squared(bg, 9.0, -0.7)
        @test s0_squared(bg, 9.0, 0.7) ≈ s0_squared(bg, 9.0, 0.7 + 2π) rtol = 1e-14
        # ε_soft enters as +ε² (59)
        bg0 = DiscBackground(; χ = 0.9, aEMRI = 8.0, h0 = 0.1, b = 0.0)
        @test s0_squared(bg, 8.0, 0.0) ≈ s0_squared(bg0, 8.0, 0.0) + softening_length(bg)^2 atol = 1e-13
        @test s0_squared(bg0, 8.0, 0.0) ≈ 0 atol = 1e-10                         # on the worldline, unsoftened
        # analytic ∂_r s0² vs FD
        for r in (6.5, 8.0, 11.0), Ψ in (0.0, 1.0, 2.5)
            @test s0_squared_dr(bg, r, Ψ) ≈ _fd(rr -> s0_squared(bg, rr, Ψ), r) rtol = 1e-9
        end
        # near the worldline the :stable form keeps eps relative accuracy, the :printed regrouping loses
        # ~(a•/ε_soft)² eps (BigFloat reference; a• = 20, h0 b = 0.015 → ε = 0.3)
        bgs = DiscBackground(; χ = -0.5, aEMRI = 20.0, h0 = 0.05, b = 0.3)
        bgsB = DiscBackground{BigFloat}(; χ = -1 // 2, aEMRI = 20, h0 = 5 // 100, b = 3 // 10)
        for Ψ in (0.0, 0.01, 0.05)
            ref = Float64(s0_squared(bgsB, big(20), big(Ψ)))
            @test abs(s0_squared(bgs, 20.0, Ψ) - ref) ≤ 8 * eps(ref)
            @test abs(s0_squared(bgs, 20.0, Ψ; form = :printed) - ref) ≤ 1e-12 * ref   # (a/ε)² eps ≈ 1e-12
        end
        @test abs(s0_squared(bgs, 20.0, 0.0; form = :printed) - 0.09) > 10 * eps(0.09)    # the loss is real
        @test s0_squared(bgs, 20.0, 0.0) == 0.3^2                                           # ε² exactly
        @test s0_squared_dr(bgs, 20.0, 0.0) == 0.0
    end

    @testset "(viii) (52): nonzero pattern, constant tensor K_μν, pointwise field" begin
        K = singular_field_coeffs(bg)
        w = worldline_velocity(bg)
        pot = kerr_equatorial_potentials(bg.χ, bg.aEMRI)
        gtt = -exp(2pot.ν) + exp(2pot.ψ) * pot.ω̃^2
        gtφ = -exp(2pot.ψ) * pot.ω̃
        gφφ = exp(2pot.ψ)
        grr = exp(2pot.μ)
        @test K.Ktt ≈ 2 * (gtt + 2 * w.vt^2) rtol = 1e-14
        @test K.Ktphi ≈ 2 * (gtφ + 2 * w.vt * w.vphi) rtol = 1e-14
        @test K.Kphiphi ≈ 2 * (gφφ + 2 * w.vphi^2) rtol = 1e-14
        @test K.Krr ≈ 2 * grr rtol = 1e-14
        @test K.Kzz == 2
        # against the dense BL metric (independent route; order t,r,θ,φ)
        g = kerr_metric(bg.χ, bg.aEMRI, π / 2)
        @test K.Ktt ≈ 2 * (g[1, 1] + 2 * w.vt^2) rtol = 1e-13
        @test K.Kphiphi ≈ 2 * (g[4, 4] + 2 * w.vphi^2) rtol = 1e-13
        @test K.Krr ≈ 2 * g[2, 2] rtol = 1e-13
        @test K.εsoft == softening_length(bg)
        # pointwise h_μν(r, Ψ) (52), m• = 1: K_μν / s0
        h = singular_field(bg, 9.0, 0.4)
        s0 = sqrt(s0_squared(bg, 9.0, 0.4))
        @test h.htt ≈ K.Ktt / s0 rtol = 1e-14
        @test h.htphi ≈ K.Ktphi / s0 rtol = 1e-14
        @test h.hphiphi ≈ K.Kphiphi / s0 rtol = 1e-14
        @test h.hrr ≈ K.Krr / s0 rtol = 1e-14
        @test h.hzz ≈ 2 / s0 rtol = 1e-14
        @test h.htt > 0                                   # g_tt + 2v_t² > 0 (Φ = −h_tt/2 < 0, attractive)
        # mass scaling: per unit mass; q enters linearly
        @test singular_field(bg, 9.0, 0.4; mEMRI = 0.5).htt ≈ h.htt / 2 rtol = 1e-15
        @test singular_field_mode(bg, 3, 9.0; mEMRI = 2.0).htt ≈ 2 * singular_field_mode(bg, 3, 9.0).htt rtol = 1e-15
    end

    @testset "(ii) Newtonian limit: ĥ_tt,m/(−2) → softened Plummer mode via laplace_b" begin
        res = Dict{Float64,Float64}()
        for a in (1e4, 1e6)
            bgN = DiscBackground(; χ = 0.0, aEMRI = a, h0 = 0.1, b = 0.6)
            ε = softening_length(bgN)
            worst = 0.0
            for m in (0, 1, 2, 5, 10, 20), rfac in (0.8, 0.97, 1.0, 1.03, 1.5)
                r = rfac * a
                hm = singular_field_mode(bgN, m, r)
                Φm = _plummer_mode(m, r, a, ε)
                @test Φm < 0
                @test hm.htt > 0
                rel = abs(hm.htt / (-2) - Φm) / abs(Φm)
                worst = max(worst, rel)
                # derivative: ∂_r of the Plummer mode by FD of the laplace route
                dΦm = _fd(rr -> _plummer_mode(m, rr, a, ε), r)
                @test abs(hm.htt_r / (-2) - dΦm) ≤ 2e-3 * (abs(dΦm) + abs(Φm) / r)   # O(M/a•) of the derivative scale
            end
            res[a] = worst
            println("  (ii) Newtonian limit a• = $a: max relative |ĥ_tt,m/(−2) − Φ̂_m(Plummer, laplace_b)| = ", worst)
            @test worst < 2e-3 * (1e4 / a)   # O(M/a•)
        end
        # the residual is O(M/a•): ratio ≈ 100 between a• = 1e4 and 1e6 (allow a generous factor)
        @test 20 < res[1e4] / res[1e6] < 500
        # other components at leading order: h_rr = h_zz = h_tt, h_φφ = a² h_tt, h_tφ = O(√(M/a•)) (CAS report)
        bgN = DiscBackground(; χ = 0.0, aEMRI = 1e6, h0 = 0.1, b = 0.6)
        hm = singular_field_mode(bgN, 2, 1.1e6)
        @test hm.hrr ≈ hm.htt rtol = 1e-5
        @test hm.hzz ≈ hm.htt rtol = 1e-5
        @test hm.hphiphi ≈ (1e6)^2 * hm.htt rtol = 1e-5
        @test abs(hm.htphi) / (1e6 * hm.htt) < 5e-3            # |h_tφ|/(a h_tt) ~ 2√(M/a) ≈ 2e-3
    end

    @testset "(iii) m-decay beyond m ≈ a•/ε_soft" begin
        for (χ, a, h0, b) in ((0.9, 8.0, 0.1, 0.6), (0.0, 10.0, 0.05, 0.3))
            bgd = DiscBackground(; χ = χ, aEMRI = a, h0 = h0, b = b)
            mstar = ceil(Int, a / softening_length(bgd))
            for r in (a, 1.15a, 0.9a)
                vals = [singular_field_mode(bgd, m, r).htt for m in mstar:2:3mstar]
                @test all(vals .> 0)                                  # Laplace-like: positive cosine modes
                @test all(diff(abs.(vals)) .< 0)                      # strictly decreasing
                # e^{−mδ}-type decay: between m* ≈ 1/δ and 3m* the modes drop by ≳ e^{−2} (≈ 0.09 at r = a•,
                # faster away from it where the strip is wider)
                @test abs(vals[end]) < 0.15 * abs(vals[1])
            end
        end
    end

    @testset "(iv) trapezoid N-convergence (default N vs 2N) for b ≥ 0.3" begin
        for (χ, a, h0, b) in ((0.9, 8.0, 0.1, 0.6), (0.0, 10.0, 0.1, 0.3), (-0.5, 20.0, 0.05, 0.3))
            bgc = DiscBackground(; χ = χ, aEMRI = a, h0 = h0, b = b)
            for m in (0, 1, 7, 40), r in (a, 1.02a, 0.8a, 2.5a)
                N0 = trapezoid_nodes(bgc, m, r)
                @test N0 ≥ max(256, 8m + 64)
                @test iseven(N0)
                h1 = singular_field_mode(bgc, m, r)
                h2 = singular_field_mode(bgc, m, r; N = 2N0)
                h0 = singular_field_mode(bgc, 0, r)                    # the m = 0 mode sets the roundoff scale
                @test h1.N == N0
                for k in _SRC_KEYS
                    # spectral: 1e-13 relative to the mode itself, down to the roundoff floor of the N-term sum
                    # (4 eps of the largest mode of the same component)
                    @test abs(getproperty(h1, k) - getproperty(h2, k)) ≤ 1e-13 * abs(getproperty(h2, k)) + 4e-15 * abs(getproperty(h0, k))
                end
            end
        end
        # the strip half-width is a property of s0² only: δ ≈ ε/a• at r = a•, larger away from it;
        # halving b halves δ (→ N grows ∝ 1/b: b → 0 breaks the spectral convergence)
        bg1 = DiscBackground(; χ = 0.9, aEMRI = 8.0, h0 = 0.1, b = 0.6)
        bg2 = DiscBackground(; χ = 0.9, aEMRI = 8.0, h0 = 0.1, b = 0.3)
        δ1 = strip_halfwidth(bg1, 8.0); δ2 = strip_halfwidth(bg2, 8.0)
        @test 0.8 * 0.06 < δ1 < 0.06
        @test δ2 / δ1 ≈ 0.5 rtol = 0.05
        @test strip_halfwidth(bg1, 12.0) > 3δ1
        @test trapezoid_nodes(bg2, 10, 8.0) > trapezoid_nodes(bg1, 10, 8.0)
        # under-resolved N is visibly worse (the default is not accidental)
        hN = singular_field_mode(bg2, 5, 8.0; N = 64)
        hR = singular_field_mode(bg2, 5, 8.0)
        @test abs(hN.htt - hR.htt) > 1e-6 * abs(hR.htt)
    end

    @testset "(vi) analytic ∂_r of every mode vs FD (Richardson)" begin
        for (χ, a, h0) in ((0.9, 8.0, 0.1), (-0.5, 15.0, 0.05))
            bgd = DiscBackground(; χ = χ, aEMRI = a, h0 = h0, b = 0.6)
            for m in (0, 3, 12), r in (0.85a, a, 1.3a)
                hm = singular_field_mode(bgd, m, r)
                for (k, kr) in ((:htt, :htt_r), (:htphi, :htphi_r), (:hphiphi, :hphiphi_r), (:hrr, :hrr_r), (:hzz, :hzz_r))
                    fd = _fd(rr -> getproperty(singular_field_mode(bgd, m, rr; N = hm.N), k), r)
                    @test abs(getproperty(hm, kr) - fd) ≤ 1e-8 * abs(fd) + 1e-14 * abs(getproperty(hm, k) / r)
                end
                @test hm.I_r ≈ _fd(rr -> singular_field_mode(bgd, m, rr; N = hm.N).I, r) rtol = 1e-8
            end
        end
    end

    @testset "(v) BigFloat(256) path" begin
        setprecision(BigFloat, 256) do
            bgB = DiscBackground{BigFloat}(; χ = 0.9, aEMRI = 8, h0 = 1 // 10, b = 6 // 10)
            hB = singular_field_mode(bgB, 2, big(9))
            hF = singular_field_mode(bg, 2, 9.0)
            @test hB.htt isa BigFloat && hB.htt_r isa BigFloat && hB.I isa BigFloat
            for k in _SRC_KEYS
                @test Float64(getproperty(hB, k)) ≈ getproperty(hF, k) rtol = 1e-13
            end
            hB2 = singular_field_mode(bgB, 2, big(9); N = 2hB.N)
            @test abs(hB2.htt - hB.htt) < big(10)^-60 * abs(hB.htt)
            @test s0_squared(bgB, big(9), big(4) / 10) isa BigFloat
            @test abs(s0_squared(bgB, big(9), big(4) / 10) - s0_squared(bg, 9.0, 0.4)) < 1e-14 * s0_squared(bg, 9.0, 0.4)
            @test worldline_velocity(bgB).ut isa BigFloat
            @test raise_modes(bgB, big(9), hB).htt isa BigFloat
        end
    end

    @testset "(ix) source_modes: grid == pointwise bitwise; threads == serial bitwise" begin
        rgrid = collect(range(5.0, 20.0, length = 257))
        m = 5
        S1 = source_modes(bg, m, rgrid)
        @test length(S1) == length(rgrid)
        @test eltype(S1) <: NamedTuple
        for (i, r) in enumerate(rgrid)
            p = singular_field_mode(bg, m, r)
            for k in _SRC_KEYS
                @test getproperty(S1[i], k) === getproperty(p, k)      # bitwise
            end
        end
        S2 = source_modes(bg, m, rgrid; threads = true)
        for i in eachindex(rgrid), k in _SRC_KEYS
            @test getproperty(S2[i], k) === getproperty(S1[i], k)
        end
        # in-place variant, allocation-free per point
        out = similar(S1)
        source_modes!(out, bg, m, rgrid)
        @test out == S1
        alloc = @allocated source_modes!(out, bg, m, rgrid)
        @test alloc ≤ 256                                               # no per-point allocation
        # the modes feed master_inputs directly (same property names as the generated kernels)
        inp = master_inputs(bg, m, rgrid[100]; hmodes = S1[100], landau = false)
        @test inp.htt == S1[100].htt && inp.hzz == S1[100].hzz
    end

    @testset "(x) contravariant modes via the background metric at the field point" begin
        for r in (6.0, 9.5, 14.0), m in (0, 4)
            hm = singular_field_mode(bg, m, r)
            up = raise_modes(bg, r, hm)
            gi = kerr_metric_inv(bg.χ, r, π / 2)                   # (t, r, θ, φ); z-direction: g^{zz} = 1
            # covariant 4×4 (t, r, φ, z) of the modes
            H = [hm.htt 0 hm.htphi 0; 0 hm.hrr 0 0; hm.htphi 0 hm.hphiphi 0; 0 0 0 hm.hzz]
            G = [gi[1, 1] 0 gi[1, 4] 0; 0 gi[2, 2] 0 0; gi[4, 1] 0 gi[4, 4] 0; 0 0 0 1.0]
            Hup = G * H * G
            @test up.htt ≈ Hup[1, 1] rtol = 1e-12
            @test up.htphi ≈ Hup[1, 3] rtol = 1e-12
            @test up.hphiphi ≈ Hup[3, 3] rtol = 1e-12
            @test up.hrr ≈ Hup[2, 2] rtol = 1e-12
            @test up.hzz ≈ Hup[4, 4] rtol = 1e-12
            # r-derivatives of the raised modes vs FD of the whole map
            for (k, kr) in ((:htt, :htt_r), (:htphi, :htphi_r), (:hphiphi, :hphiphi_r), (:hrr, :hrr_r), (:hzz, :hzz_r))
                fd = _fd(rr -> getproperty(raise_modes(bg, rr, singular_field_mode(bg, m, rr; N = hm.N)), k), r)
                @test abs(getproperty(up, kr) - fd) ≤ 1e-7 * abs(fd) + 1e-14 * abs(getproperty(up, k) / r)
            end
            # an even-in-Ψ raising preserves the cosine projection: raising the pointwise field and
            # projecting equals raising the projected modes (the background does not depend on Ψ)
            N = hm.N
            acc = 0.0
            for k in 0:N-1
                Ψ = 2π * k / N
                hp = singular_field(bg, r, Ψ)
                Hp = [hp.htt 0 hp.htphi 0; 0 hp.hrr 0 0; hp.htphi 0 hp.hphiphi 0; 0 0 0 hp.hzz]
                acc += (G * Hp * G)[1, 1] * cos(m * Ψ)
            end
            @test acc / N ≈ up.htt rtol = 1e-12
        end
    end

    @testset "(xi) end-to-end: relativistic S with singular-field modes → Newtonian S (50f-g) as a• → ∞" begin
        for (a, tol) in ((1e4, 5e-3), (1e6, 5e-5))
            bgN = DiscBackground(; χ = 0.0, aEMRI = a, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0, b = 0.6)
            for m in (2, 6), rfac in (0.9, 1.3)
                r = rfac * a
                hm = singular_field_mode(bgN, m, r)
                C, S = master_matrices(bgN, m, r; hmodes = hm, landau = false)
                CN, SN = master_matrices_newtonian(bgN, m, r, -hm.htt / 2, -hm.htt_r / 2; landau = false)
                @test abs(S[1] - SN[1]) ≤ tol * abs(SN[1])
                @test abs(S[2] - SN[2]) ≤ tol * abs(SN[2])
            end
        end
    end
end
