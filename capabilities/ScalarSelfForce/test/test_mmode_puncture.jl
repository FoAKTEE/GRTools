# test_mmode_puncture.jl — ESSF node `mmode-puncture-integrals` (schedule W4-D), written
# BEFORE the implementation (TDD). Target: Puncture/mmode_integrals.jl — the m-mode puncture
# Ψ^P_m from the 4D Φ^P via the CRP ₃P route (arXiv:2606.04998 main.tex:459-486; Bourg
# arXiv:2404.10082 Eq. (115)). Gates (GLOBAL_DAG #34 oracle column + worker contract,
# calibrated by the pre-test prototype measurements recorded in the node evidence):
#   * closed-form identity ∫R^n e^{-imΔφ}dΔφ vs in-test trapezoid quadrature, n odd ∈ [-7,9],
#     m ∈ {0,1,2,5,10,20}, ϱ/z_c ∈ {1e-2,0.1,1,10} × two z_c² scales; Float64 with an
#     explicit trapezoid noise floor 1e-12·|I_0| (skip rule m·μ ≤ 14, μ = 2asinh(ϱ/z_c),
#     the e^{-mμ} Fourier-decay scale), BigFloat(256) spot checks for the high-m corners
#     the Float64 quadrature cannot see;
#   * characterization vs BOTH baked effsource_PhiS_m fixtures (test/fixtures/effsource/),
#     ALL 6 frozen configs incl. a = -0.5 / -0.998, at the capture driver's effective offsets:
#     - comoving (kerr-circular-comoving.c — the SAME ₃P/2F1 algorithm, independent C/GSL
#       implementation): rtol 1e-10 away from the worldline, 1e-9 near (C1 = ϱ²/z_c² < 1e-8),
#       ALL fixture m ∈ {-2,-1,0,1,2,5,10,20};
#     - elliptic (kerr-circular.c — INDEPENDENT Thornburg E/K-polynomial algorithm): rtol
#       3e-9 in the header's overlap window C1 ∈ [1e-8,1e-2); 5e-3 near the worldline
#       (header-stated elliptic E/K digit loss, measured 1.2e-3); far field (C1 ≥ 1e-2)
#       1e-10 for m ≤ 1 only — for m ≥ 2 the ELLIPTIC path itself loses digits far out
#       (measured 1.2e-5 at m=2, total loss by m=5; the CRP's stated reason for switching
#       routes, main.tex:465-466), so far-field m ≥ 2 is adjudicated against quadrature:
#       ours beats the C polynomial route (GLOBAL_DAG #34), Float64 trapezoid for m ∈ {2,5},
#       BigFloat(256) trapezoid of the BigFloat 4D puncture for m ∈ {10,20} (the integral is
#       e^{-mμ}-small there, below the Float64 quadrature noise floor);
#   * brute-force cross-check at ONE moderate point per config (worker contract): adaptive
#     trapezoid refinement of the 4D puncture to 1e-11 relative stall, all m ∈ {0,1,2,5,10,20};
#   * Ψ^P_m assembly: exact relation to the bare integral through the LOCKED ansatz layer
#     (mode_amplitude), slicing-time invariance of the co-rotating assembly, fixture bridge,
#     conjugate symmetry Ψ^P_{-m} = conj(Ψ^P_m), m-evenness of the bare integrals;
#   * A_ijk Q-power regrouping vs the solid 4D puncture_phi (same 50-coefficient knowledge);
#   * BigFloat(256) smoke: near-worldline point (worst-case cancellation candidate — measured
#     loss is BENIGN there, ~1e-16: α ≪ β kills the alternating (α/β)-weighted sums) + the
#     measured WORST Float64-loss grid point (far field, extreme C1 ~ 120: cfg4 dr=10 m=20,
#     measured 5.7e-11) within documented loss; full BigFloat type inheritance.
using Test
using ScalarSelfForce
using ScalarSelfForce.BigSpecial: legendre_P3, pochhammer

# fixture constants may already be loaded by earlier test files under Pkg.test; guard the
# includes so this file also runs standalone.
if !@isdefined(EFFSOURCE_CONFIGS)
    include(joinpath(@__DIR__, "fixtures", "effsource", "effsource_configs.jl"))
end
if !@isdefined(EFFSOURCE_PHIS_M_COMOVING)
    include(joinpath(@__DIR__, "fixtures", "effsource", "phis_m_comoving.jl"))
end
if !@isdefined(EFFSOURCE_PHIS_M_ELLIPTIC)
    include(joinpath(@__DIR__, "fixtures", "effsource", "phis_m_elliptic.jl"))
end

# ---- in-test brute-force oracles (unique _mmq_ prefix; no new deps) ----------------------

# Trapezoid refinement of ∫_{-π}^{π} f(Δφ) e^{-imΔφ} dΔφ for periodic f — spectrally
# convergent; refined until successive doublings agree to `tol` relative (the worker-contract
# 1e-11 stall criterion) or kmax is hit (Float64 noise-floor cases assert looser gates).
function _mmq_trapz(::Type{T}, f, m; tol, kmax) where {T}
    prev = complex(zero(T))
    for k in 6:kmax
        N = 2^k
        s = complex(zero(T))
        for j in 0:(N - 1)
            phi = -T(pi) + 2 * T(pi) * j / N
            s += f(phi) * cis(-m * phi)
        end
        cur = s * (2 * T(pi) / N)
        k > 6 && abs(cur - prev) <= tol * abs(cur) && return cur
        prev = cur
    end
    return prev
end

# ∫ R^n e^{-imΔφ} dΔφ by quadrature, R² = ϱ² + z_c² sin²(Δφ/2)
_mmq_Rn_quad(::Type{T}, n, m, rho2, zc2; tol, kmax = 19) where {T} =
    _mmq_trapz(T, phi -> (rho2 + zc2 * sin(phi / 2)^2)^(T(n) / 2), m; tol, kmax)

# ∫ Φ^P e^{-imΔφ} dΔφ by quadrature of the SOLID 4D puncture (puncture-field-4d node)
_mmq_punc_quad(::Type{T}, p, m, dr, dth; tol, kmax = 19) where {T} =
    _mmq_trapz(T, phi -> puncture_phi(p, dr, dth, phi), m; tol, kmax)

# capture-driver effective offsets (generator/capture.c:290: x = (r0 + dr, π/2 + dθ, φ_p)),
# same diagnosis as test_puncture_field4d.jl
_mmq_eff(r0, dr, dth) = ((r0 + dr) - r0, (pi / 2 + dth) - pi / 2)

# fixture row → (cfg, m, dr, dθ, complex oracle);  columns per the fixture header
_mmq_row(fix, col) = (Int(fix[1, col]), Int(fix[2, col]), fix[4, col], fix[5, col],
                      complex(fix[6, col], fix[7, col]))

@testset "Puncture: m-mode integrals + Ψ^P_m (mmode-puncture-integrals)" begin

    orbits = Dict(c.idx => CircularOrbit(c.a, c.r0) for c in EFFSOURCE_CONFIGS)
    puncs = Dict(c.idx => Puncture(orbits[c.idx]) for c in EFFSOURCE_CONFIGS)
    r0s = Dict(c.idx => c.r0 for c in EFFSOURCE_CONFIGS)

    @testset "mode_integral_Rn: Bourg Eq.(115) closed form vs trapezoid (n odd ∈ [-7,9])" begin
        # Float64 leg over m·μ ≤ 14 (μ = 2asinh(ϱ/z_c), the Fourier decay rate of the
        # R² = ϱ²+z_c²sin²(Δφ/2) singularity structure); z_c² of configs 1 (a=0) and 4
        # (a=0.998 ISCO). Gate: rel 1e-9 OR the Float64-trapezoid noise floor 1e-12·|I_0|
        # (Fourier coefficients decay below eps·integrand-scale for smooth R^{n>0} at
        # high m; roundoff of the N-term sum ~ eps·√N·|I_0| ≤ 7e-14·|I_0| measured).
        worst = 0.0
        ncombo = 0
        for zc2 in (457.14285714285717, 24.294546275416153),
            ratio in (1e-2, 0.1, 1.0, 10.0), n in -7:2:9

            rho2 = ratio^2 * zc2
            floor0 = 1e-12 * abs(_mmq_Rn_quad(Float64, n, 0, rho2, zc2; tol = 1e-12))
            for m in (0, 1, 2, 5, 10, 20)
                m * 2 * asinh(ratio) <= 14 || continue
                ncombo += 1
                v = mode_integral_Rn(n, m, rho2, zc2)
                q = _mmq_Rn_quad(Float64, n, m, rho2, zc2; tol = 1e-12)
                @test abs(v - real(q)) <= max(1e-9 * abs(q), floor0)
                @test abs(imag(q)) <= max(1e-9 * abs(q), floor0)  # genuinely real
                worst = max(worst, abs(v - real(q)) / max(abs(q), floor0))
            end
        end
        println("mmode-puncture Rn identity [F64]: ", ncombo, " combos, worst rel dev = ", worst)

        # BigFloat(256) leg: the high-m corners invisible to Float64 quadrature (integral
        # e^{-mμ}-small, incl. the two m·μ > 14 skips and the worst F64-floor combo).
        # Gate 1e-20: bounded by the BigFloat(256) quadrature floor eps·|I_0|/|I_m| ≤
        # 1e-77·1e52 ≈ 1e-25 (measured worst 2.2e-23), 10+ orders below any Float64 scale.
        setprecision(BigFloat, 256) do
            zc2 = big"457.14285714285717"
            for (n, m, ratio) in ((-1, 20, 1.0), (-7, 10, 10.0), (1, 20, 10.0), (9, 20, 0.1))
                rho2 = big(ratio)^2 * zc2
                v = mode_integral_Rn(n, m, rho2, zc2)
                q = _mmq_Rn_quad(BigFloat, n, m, rho2, zc2; tol = big"1e-30", kmax = 14)
                @test isapprox(v, real(q); rtol = big"1e-20")
            end
        end

        # evenness in m (theorem: P^{-m}/(1+n/2)_{-m} = P^m/(1+n/2)_m; the -m code path is
        # different arithmetic, so equality is roundoff-level, gate 1e-13)
        for n in -7:2:9, m in (1, 2, 5, 10, 20)
            v = mode_integral_Rn(n, m, 4.31, 457.14285714285717)
            @test isapprox(mode_integral_Rn(n, -m, 4.31, 457.14285714285717), v; rtol = 1e-13)
        end
    end

    @testset "A_ijk Q-power regrouping vs solid 4D puncture_phi (rtol 1e-13)" begin
        # Σ_q A_q(Δr,Δθ) Q^{2q} / R⁷ must reproduce the puncture-field-4d evaluation — the
        # two homes of the 50-coefficient polynomial knowledge stay consistent.
        for cfg in 1:6, (dr, dth, dphi) in ((0.01, 0.007, 0.3), (-0.05, 0.02, 1.1),
                                            (0.5, 0.2, 2.5), (-1e-6, 0.0, 0.05))
            p = puncs[cfg]
            o = orbits[cfg]
            A = ScalarSelfForce._puncture_A_groups(p, dr, dth)
            Q2 = sin(dphi / 2)^2
            R2 = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dth^2 + o.zc2 * Q2
            v = sum(A[q + 1] * Q2^q for q in 0:4) / (R2^3 * sqrt(R2))
            @test isapprox(v, puncture_phi(p, dr, dth, dphi); rtol = 1e-13)
        end
    end

    @testset "characterization vs comoving PhiS_m (same ₃P route, independent C/GSL impl)" begin
        # gate: 1e-10 away from the worldline / 1e-9 near (C1 = ϱ²/z_c² < 1e-8); measured
        # worst 2.8e-12 away, 1.9e-11 near. ALL fixture m incl. the conjugate rows -2, -1.
        worst_near = 0.0
        worst_away = 0.0
        for col in axes(EFFSOURCE_PHIS_M_COMOVING, 2)
            cfg, m, dr, dth, oracle = _mmq_row(EFFSOURCE_PHIS_M_COMOVING, col)
            p = puncs[cfg]
            o = orbits[cfg]
            dre, dthe = _mmq_eff(r0s[cfg], dr, dth)
            v = puncture_phi_m(p, m, dre, dthe) * cis(-m * EFFSOURCE_PHIP)
            C1 = (o.rho2_dr2 * dre^2 + o.rho2_dtheta2 * dthe^2) / o.zc2
            rd = abs(v - oracle) / abs(oracle)
            if C1 < 1e-8
                @test isapprox(v, oracle; rtol = 1e-9)
                worst_near = max(worst_near, rd)
            else
                @test isapprox(v, oracle; rtol = 1e-10)
                worst_away = max(worst_away, rd)
            end
        end
        println("mmode-puncture vs comoving: worst rel dev near = ", worst_near,
                ", away = ", worst_away)
    end

    @testset "characterization vs elliptic PhiS_m (independent E/K algorithm), by window" begin
        # windows per the effsource_configs.jl header conditioning notes + prototype
        # adjudication (see file header): overlap C1 ∈ [1e-8,1e-2) both oracles good →
        # 3e-9 (header gate); near C1 < 1e-8 elliptic loses digits → 5e-3; far C1 ≥ 1e-2
        # elliptic is trustworthy only for m ≤ 1 → 1e-10 there, m ≥ 2 handled by the
        # quadrature adjudication testset below.
        worst_mid = 0.0
        for col in axes(EFFSOURCE_PHIS_M_ELLIPTIC, 2)
            cfg, m, dr, dth, oracle = _mmq_row(EFFSOURCE_PHIS_M_ELLIPTIC, col)
            p = puncs[cfg]
            o = orbits[cfg]
            dre, dthe = _mmq_eff(r0s[cfg], dr, dth)
            v = puncture_phi_m(p, m, dre, dthe) * cis(-m * EFFSOURCE_PHIP)
            C1 = (o.rho2_dr2 * dre^2 + o.rho2_dtheta2 * dthe^2) / o.zc2
            if C1 < 1e-8
                @test isapprox(v, oracle; rtol = 5e-3)
            elseif C1 < 1e-2
                @test isapprox(v, oracle; rtol = 3e-9)
                worst_mid = max(worst_mid, abs(v - oracle) / abs(oracle))
            elseif abs(m) <= 1
                @test isapprox(v, oracle; rtol = 1e-10)
            end
        end
        println("mmode-puncture vs elliptic overlap window: worst rel dev = ", worst_mid)
    end

    @testset "far-field m ≥ 2: ours beats the C polynomial route (vs quadrature)" begin
        # GLOBAL_DAG #34: far from the worldline the elliptic E/K-polynomial fixture is the
        # NOISY one (the CRP's reason for the ₃P route, main.tex:465-466). Adjudicate at the
        # (dr = 10, dθ = 0) grid extreme of every config: |ours − quad| < |elliptic − quad|
        # strictly, plus a direct ours-vs-quad gate. m ∈ {2,5}: Float64 trapezoid (gates
        # 1e-10 / 5e-3 — the m=5 integral sits near the Float64 quadrature noise floor);
        # m ∈ {10,20}: BigFloat(256) trapezoid of the BigFloat 4D puncture at cfg ∈ {1,4}
        # (the integral is e^{-mμ}-small, invisible to Float64 quadrature), ours-vs-quad
        # gate 1e-9 (measured worst Float64 combination loss 5.7e-11 at cfg4).
        ell_at = Dict{Tuple{Int,Int},Float64}()
        for col in axes(EFFSOURCE_PHIS_M_ELLIPTIC, 2)
            cfg, m, dr, dth, oracle = _mmq_row(EFFSOURCE_PHIS_M_ELLIPTIC, col)
            (dr == 10.0 && dth == 0.0) || continue
            ell_at[(cfg, m)] = real(oracle * cis(m * EFFSOURCE_PHIP))  # bare real integral
        end
        for cfg in 1:6, m in (2, 5)
            p = puncs[cfg]
            v = puncture_phi_m(p, m, 10.0, 0.0)
            q = real(_mmq_punc_quad(Float64, p, m, 10.0, 0.0; tol = 1e-11))
            @test isapprox(v, q; rtol = m == 2 ? 1e-10 : 5e-3)
            @test abs(v - q) < abs(ell_at[(cfg, m)] - q)
        end
        setprecision(BigFloat, 256) do
            for cfg in (1, 4), m in (10, 20)
                c = EFFSOURCE_CONFIGS[cfg]
                pb = Puncture(CircularOrbit(big(c.a), big(c.r0)))
                q = real(_mmq_punc_quad(BigFloat, pb, m, 10.0, 0.0;
                                        tol = big"1e-25", kmax = 14))
                v = puncture_phi_m(puncs[cfg], m, 10.0, 0.0)
                @test isapprox(v, Float64(q); rtol = 1e-9)
                @test abs(v - q) < abs(ell_at[(cfg, m)] - q)
            end
        end
    end

    @testset "brute-force cross-check: one moderate point per config, all m" begin
        # worker-contract gate: trapezoid refinement of the 4D puncture to 1e-11 relative
        # stall at the diagonal grid point (0.5, 0.2); rtol 1e-8 (measured worst 1.5e-9 at
        # cfg4 m=20, where the Float64 quadrature is the accuracy limiter).
        for cfg in 1:6, m in (0, 1, 2, 5, 10, 20)
            p = puncs[cfg]
            v = puncture_phi_m(p, m, 0.5, 0.2)
            q = _mmq_punc_quad(Float64, p, m, 0.5, 0.2; tol = 1e-11)
            @test abs(v - q) <= 1e-8 * abs(q)       # complex compare: pins Im(q) ≈ 0 too
        end
    end

    @testset "Ψ^P_m assembly per the LOCKED mode convention" begin
        # (a) exact relation to the bare integral through the ansatz layer: Ψ^P_m =
        #     mode_amplitude(m, I_m e^{-imφ̄_p}/2π; ts = 0) — and the same through ts ≠ 0
        #     with the CO-ROTATING phase φ̄_p(ts) = φ̄_p + Ω ts (slicing-time invariance:
        #     e^{+imΩts} from eq:modeIntegral cancels the puncture's rigid rotation).
        for cfg in (2, 5), (dr, dth) in ((0.01, 0.0), (0.5, 0.2), (-0.1, -0.05)),
            m in (0, 1, 2, 5, 20)

            p = puncs[cfg]
            o = orbits[cfg]
            r = o.r0 + dr
            z = cos(pi / 2 + dth)
            psi = puncture_psi_m(p, m, r, z; varphi_p = 0.31)
            I = puncture_phi_m(p, m, dr, dth)
            @test isapprox(psi,
                           mode_amplitude(m, I * cis(-m * 0.31) / (2pi);
                                          Omega = o.Omega, ts = 0.0, r = r, z = z);
                           rtol = 1e-13)
            ts = 0.7
            @test isapprox(psi,
                           mode_amplitude(m, I * cis(-m * (0.31 + o.Omega * ts)) / (2pi);
                                          Omega = o.Omega, ts = ts, r = r, z = z);
                           rtol = 1e-12)
            # (b) conjugate symmetry Ψ^P_{-m} = conj(Ψ^P_m) (bare integral even in m + the
            #     |m| ansatz factors + conjugated phase)
            @test isapprox(puncture_psi_m(p, -m, r, z; varphi_p = 0.31), conj(psi);
                           rtol = 1e-12)
        end
        # (c) fixture bridge: the baked comoving Ψ^P_m rows ARE I_m e^{-imφ_p} at φ_p = 0.31,
        #     so mode_amplitude over fixture/2π equals puncture_psi_m(...; varphi_p = 0.31)
        #     up to the characterization tolerance (end-to-end factor placement pin).
        for col in axes(EFFSOURCE_PHIS_M_COMOVING, 2)
            cfg, m, dr, dth, oracle = _mmq_row(EFFSOURCE_PHIS_M_COMOVING, col)
            (cfg == 3 && dr in (0.001, -0.001) && dth == 0.0 && m in (0, 2, 20)) || continue
            o = orbits[cfg]
            dre, dthe = _mmq_eff(r0s[cfg], dr, dth)
            r = o.r0 + dre
            z = cos(pi / 2 + dthe)
            @test isapprox(puncture_psi_m(puncs[cfg], m, r, z; varphi_p = EFFSOURCE_PHIP),
                           mode_amplitude(m, oracle / (2pi);
                                          Omega = o.Omega, ts = 0.0, r = r, z = z);
                           rtol = 1e-9)
        end
        # (d) m = 0 with varphi_p = 0: Ψ^P_0 is real and positive near the particle
        psi0 = puncture_psi_m(puncs[1], 0, r0s[1] + 0.1, 0.0)
        @test imag(psi0) == 0
        @test real(psi0) > 0
    end

    @testset "BigFloat(256) smoke + type inheritance" begin
        c = EFFSOURCE_CONFIGS[3]                    # a = 0.9, r0 = isco + 0.1
        p64 = puncs[3]
        setprecision(BigFloat, 256) do
            pb = Puncture(CircularOrbit(big(c.a), big(c.r0)))
            # near-worldline point (dr = -1e-6, m = 20): the α ≪ β limit kills the
            # alternating (α/β)-weighted sums — measured Float64 loss ~4e-17, gate 1e-12
            v64 = puncture_phi_m(p64, 20, -1e-6, 0.0)
            vb = puncture_phi_m(pb, 20, big(-1e-6), big(0.0))
            @test vb isa BigFloat
            @test isapprox(v64, Float64(vb); rtol = 1e-12)
            # measured WORST Float64-loss grid point (far field at extreme C1 ~ 120:
            # cfg4 dr = 10, m = 20, measured 5.7e-11): documented-loss gate 1e-9
            c4 = EFFSOURCE_CONFIGS[4]
            pb4 = Puncture(CircularOrbit(big(c4.a), big(c4.r0)))
            v64w = puncture_phi_m(puncs[4], 20, 10.0, 0.0)
            vbw = puncture_phi_m(pb4, 20, big(10.0), big(0.0))
            @test isapprox(v64w, Float64(vbw); rtol = 1e-9)
            # Ψ assembly inherits Complex{BigFloat}
            psib = puncture_psi_m(pb, 2, big(c.r0) + 1 // 10, big(1) / 20;
                                  varphi_p = big(31) / 100)
            @test psib isa Complex{BigFloat}
            @test isapprox(Float64(real(psib)),
                           real(puncture_psi_m(p64, 2, c.r0 + 0.1, 0.05; varphi_p = 0.31));
                           rtol = 1e-10)
        end
        # promotion smoke: BigFloat offsets through a Float64 table stay BigFloat
        @test puncture_phi_m(p64, 1, big(1) / 100, big(0)) isa BigFloat
        # mode_integral_Rn type inheritance
        @test mode_integral_Rn(-7, 5, big(4.31), big(457.0)) isa BigFloat
    end
end
