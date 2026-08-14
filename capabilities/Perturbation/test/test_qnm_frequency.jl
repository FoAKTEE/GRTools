# test_qnm_frequency.jl — Q1b: Kerr QNM frequency solver ω_{s,l,m,n}(a) (Leaver radial CF + Q1a
# complex-γ spheroidal eigenvalue). Grounded vs the in-repo h5 fixture (a0.0_l2_m2.h5) and the
# authoritative Cook-Zalutskiy radial coefficients (ref-code/qnm/qnm/radial.py, the h5 generator's
# convention); the rest are analytic self-consistency identities (spin-weight independence,
# reflection symmetry, m-degeneracy, χ-continuity, BF-256 self-consistency). Tolerances MEASURED via
# scratchpad/qqnm_phase3/probe_*.log and stated with headroom.

# ── baked h5 fixture: TeukEvolution.jl/qnm/a0.0_l2_m2.h5, groups [n=0]…[n=10] (F6 extraction
#    pattern). sha256 20e85b0ee8dff677dc7ac075a9b9081b50b2173f2be1a9cb1daedf5ceec533fa. At a=0
#    λ = l(l+1)−s(s+1) = 4 (m-independent). These are the omega dataset values baked verbatim. ──
const QNM_H5_SHA256 = "20e85b0ee8dff677dc7ac075a9b9081b50b2173f2be1a9cb1daedf5ceec533fa"
const QNM_H5_OMEGA = ComplexF64[
    0.37367168441804177 - 0.08896231568893546im,   # n=0
    0.34671099687916285 - 0.27391487529123504im,   # n=1
    0.30105345461250427 - 0.47827698322307044im,   # n=2
    0.25150496222643026 - 0.7051482024421033im,    # n=3
    0.2075145798286828  - 0.9468448909071687im,    # n=4
    0.1692994030698982  - 1.1956080540907499im,    # n=5
    0.133252340198513   - 1.4479106262097563im,    # n=6
    0.09282233536257913 - 1.703841171811082im,     # n=7
    8.718542227771156e-12 - 2.032269255277643im,   # n=8  (algebraically special, Re≈0)
    0.06326546325133285 - 2.302645555505304im,     # n=9
    0.07655379946386853 - 2.5608261891388246im,    # n=10
]
# Per-overtone tolerance vs the h5. The h5 is Float64 and its OWN accuracy degrades with overtone
# (proven by gate (a3): our BF-256 root satisfies the CF to ~1e-75 yet differs from the h5 by these
# same amounts, while matching our Float64 root to ≤2e-13). So these tolerances bound the h5's
# accuracy, not the solver's. Measured |ω_solver − ω_h5| + headroom.
const QNM_H5_TOL = [1e-12, 1e-12, 1e-11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-8, 1e-3, 1e-5, 1e-5]

# Cook-Zalutskiy radial D-coefficients — transcribed VERBATIM from ref-code/qnm/qnm/radial.py
# (D_coeffs + sing_pt_char_exps). The double-entry oracle for Perturbation._qnm_leaver_c.
function _cz_Dcoeffs(ω, a, s, m, A)
    root = sqrt(1 - a * a); r_p = 1 + root
    σ_p = (2ω * r_p - m * a) / (2root); σ_m = (2ω * (1 - root) - m * a) / (2root)
    ζ = im * ω; ξ = -s - im * σ_p; η = -im * σ_m
    p = root * ζ
    α = 1 + s + ξ + η - 2ζ + s
    γ = 1 + s + 2η
    δ = 1 + s + 2ξ
    σ = A + a * a * ω * ω - 8ω * ω + p * (2α + γ - δ) + (1 + s - 0.5 * (γ + δ)) * (s + 0.5 * (γ + δ))
    return (δ, 4p - 2α + γ - δ - 2, 2α - γ + 2, α * (4p - δ) - σ, α * (α - γ + 1))
end

@testset "QNM frequency ω_{s,l,m,n}(a) (Q1b)" begin

    # ── (a) the in-repo h5 fixture: ω(2,2,0,a=0) + the 11 overtones ──────────────────────────────
    @testset "(a) a=0 h5 fixture — 11 overtones" begin
        # (a1) the fundamental to the published digits
        ω0 = qnm_frequency(-2, 2, 2, 0, 0.0; seed=QNM_H5_OMEGA[1] * (1 + 1e-3))
        @test isapprox(ω0, QNM_H5_OMEGA[1]; atol=1e-13)
        # (a2) all 11 overtones (ninv = n), each within the h5's own accuracy
        for n in 0:10
            ω, res = qnm_frequency(-2, 2, 2, n, 0.0; seed=QNM_H5_OMEGA[n+1] * (1 + 1e-3),
                                   return_residual=true)
            @test abs(ω - QNM_H5_OMEGA[n+1]) <= QNM_H5_TOL[n+1]
            @test res <= 1e-11                              # a genuine CF root
        end
    end

    # ── (a3) BF-256 self-consistency: PROVES the h5 residual is the h5's error, not ours ─────────
    @testset "(a3) BF-256 self-consistency (n=0,3,8)" begin
        setprecision(BigFloat, 256) do
            for n in (0, 3, 8)
                f64 = qnm_frequency(-2, 2, 2, n, 0.0; seed=QNM_H5_OMEGA[n+1] * (1 + 1e-3))
                ωbf, res = qnm_frequency(-2, 2, 2, n, big(0.0);
                                         seed=Complex{BigFloat}(f64), return_residual=true)
                @test ωbf isa Complex{BigFloat}
                @test abs(ComplexF64(ωbf) - f64) <= 1e-10   # Float64 root == BF root
                @test Float64(res) <= 1e-68                 # CF defect at the BF floor
            end
        end
    end

    # ── coefficient double-entry vs Cook-Zalutskiy + the documented paper c₄ typo ────────────────
    # CONVENTION BRIDGE: the paper (Ma-Yang) c₃ takes the Teukolsky spheroidal eigenvalue λ (= the
    # Harmonics public λ, verified as the physical insertion by gates (a)+Berti); Cook-Zalutskiy
    # D[3] takes the RAW separation constant A = λ + 2mγ − γ² (γ=aω). So paper-c₃(λ) == CZ-D[3](A).
    # c₀,c₁,c₂,c₄ are A-independent and match either way.
    @testset "Leaver coeffs == Cook-Zalutskiy (double-entry) + c₄ typo" begin
        for (s, m, a, ω, λ) in ((-2, 2, 0.7, 0.53 - 0.08im, 4.9 + 0.1im),
                                (-2, 3, 0.4, 0.70 - 0.09im, 8.5 - 0.2im),
                                (+2, 2, 0.9, 0.67 - 0.06im, 3.2 + 0.3im))
            γ = a * ω
            A_cz = λ + 2 * m * γ - γ^2                  # Teukolsky λ → CZ raw separation constant
            c = Perturbation._qnm_leaver_c(s, m, a, ω, λ)
            D = _cz_Dcoeffs(ω, a, s, m, A_cz)
            for j in 1:5
                @test isapprox(c[j], D[j]; rtol=1e-12, atol=1e-12)
            end
            # the paper's c₄ (tex:1172, as written) differs from the correct (CZ) c₄ by −iω(2s+3)
            b = sqrt(1 - a^2); rp = 1 + b; σp = (ω * rp - a * m / 2) / b
            c4_paper = s + 1 - (2s + 1) * im * ω - (8ω + 2im) * σp
            @test isapprox(c[5] - c4_paper, -im * ω * (2s + 3); rtol=1e-12, atol=1e-12)
        end
    end

    # ── (e) Schwarzschild m-degeneracy (a=0 frequency is m-independent) ──────────────────────────
    @testset "(e) Schwarzschild m-degeneracy" begin
        ω22 = qnm_frequency(-2, 2, 2, 0, 0.0; seed=QNM_H5_OMEGA[1])
        ω20 = qnm_frequency(-2, 2, 0, 0, 0.0; seed=QNM_H5_OMEGA[1])
        ω2m = qnm_frequency(-2, 2, -2, 0, 0.0; seed=QNM_H5_OMEGA[1])
        @test abs(ω22 - ω20) <= 1e-13
        @test abs(ω22 - ω2m) <= 1e-13
    end

    # ── λ-convention gate: spin-weight independence ω(s=-2) == ω(s=+2) ───────────────────────────
    @testset "spin-weight independence ω(s=-2)==ω(s=+2)" begin
        for a in (0.0, 0.5, 0.9)
            ωm = qnm_frequency(-2, 2, 2, 0, a; seed=0.5 - 0.08im)
            ωp = qnm_frequency(+2, 2, 2, 0, a; seed=ωm)
            @test abs(ωm - ωp) <= 1e-12
        end
    end

    # ── (b) reflection symmetry  ω^+(s,l,m) == −conj(ω^-(s,l,−m))  at 3 χ ─────────────────────────
    @testset "(b) reflection symmetry" begin
        for a in (0.3, 0.6, 0.9)
            ωp = qnm_frequency(-2, 2, 2, 0, a; seed=0.5 - 0.08im)
            ωmir = qnm_frequency(-2, 2, -2, 0, a; seed=-conj(ωp))
            @test abs(ωp - (-conj(ωmir))) <= 1e-12
        end
    end

    # ── (c) χ-continuity across χ∈{0,0.1,…,0.9,0.95} — no branch jumps ───────────────────────────
    @testset "(c) χ-continuity (2,2,0)" begin
        χs = vcat(0.0:0.1:0.9, 0.95)
        ωprev = QNM_H5_OMEGA[1]; ωs = ComplexF64[]
        for χ in χs
            ω, res = qnm_frequency(-2, 2, 2, 0, χ; seed=ωprev, return_residual=true)
            @test res <= 1e-12                      # genuine CF root at every χ
            push!(ωs, ω); ωprev = ω
        end
        dRe = diff(real.(ωs)); dIm = diff(imag.(ωs))
        @test all(dRe .> 0)                          # Re ω monotone ↑ (prograde 220)
        @test all(dIm .> 0)                          # |Im ω| monotone ↓ (decay slows)
        @test maximum(abs.(diff(ωs))) < 0.15         # no jumps (smooth continuation)
    end

    # ── (d)+(f) CF residual == 0 at BF-256 at spot points; {T}/BF-256 inheritance ────────────────
    @testset "(d,f) CF residual → BF floor + {T}/BF-256" begin
        setprecision(BigFloat, 256) do
            @test precision(BigFloat) == 256
            for a in (big"0.7", big"0.9")
                seed = Complex{BigFloat}(qnm_frequency(-2, 2, 2, 0, Float64(a); seed=0.5 - 0.08im))
                ω, res = qnm_frequency(-2, 2, 2, 0, a; seed=seed, return_residual=true)
                @test ω isa Complex{BigFloat}
                # independent residual recompute at the returned root with the coupled angular λ
                λ = Perturbation._qnm_angular_lambda(-2, 2, 2, a, ω)
                r2 = abs(qnm_cf_residual(-2, 2, 2, a, ω, λ; ninv=0))
                @test Float64(res) <= 1e-68
                @test Float64(r2) <= 1e-68
            end
        end
    end

    # ── the 5 campaign multiplets sweep to a valid (CF-residual-zero) root at every χ ────────────
    @testset "campaign multiplets (2,2,0)/(3,3,0)/(4,4,0)/(2,0,0)/(2,-2,0) sweep" begin
        χs = vcat(0.0:0.1:0.9, 0.95)
        for (l, m) in ((2, 2), (3, 3), (4, 4), (2, 0), (2, -2))
            seed = qnm_frequency(-2, l, m, 0, 0.0)   # derived eikonal default seed (no memory value)
            @test abs(qnm_cf_residual(-2, l, m, 0.0, seed,
                                      Perturbation._qnm_angular_lambda(-2, l, m, 0.0, complex(0.0)))) <= 1e-11
            for χ in χs[2:end]
                ω, res = qnm_frequency(-2, l, m, 0, χ; seed=seed, return_residual=true)
                @test res <= 1e-11
                seed = ω
            end
        end
    end
end
