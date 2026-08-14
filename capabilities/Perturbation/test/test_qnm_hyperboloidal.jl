# test_qnm_hyperboloidal.jl — gates G-Q2H-1..5 for the Q2-H hyperboloidal QNM backend.
#   G-Q2H-1  operator transcription (double-entry α/β + Δ̃ typo-guard + a=0 == cond_bf ½-form)
#   G-Q2H-2  DECISIVE: :hyperboloidal ω == :leaver over the 5×11 campaign, ≤1e-12 (two formulations)
#   G-Q2H-3  MMA cross-oracle: the (2,2,0,a=0.7) radial eigenfunction vs the baked fixture
#   G-Q2H-4  s=+2 Hertz spectrum == s=−2 (Teukolsky s-pairing, an unfakeable internal invariant)
#   G-Q2H-5  spectral-convergence curve (n ∈ {24,32,48,64}, BF-256) + {T}/BF eltype cleanliness
using Test
using Perturbation
using Harmonics: SpinWeightedSpheroidalEigenvalue
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "fixtures", "qnm_hphc_mma_ref.jl"))
const MMA = QNMHphcMMARef

setprecision(BigFloat, 256)

# The VALIDATED a=0 operator (cond_bf.jl coeffs, transcribed from QNM.wl ScriptCapitalM → R22 form).
# The Q2-H ½-form operator must equal HALF of this (== R22/2), element-wise (G-Q2H-1).
function _condbf_coeffs(n::Integer; s::Integer=-2, lam::Real=4.0)
    N = n - 1
    x = [cos(pi * j / N) for j in 0:N]
    c = [(j == 0 || j == N) ? 2.0 : 1.0 for j in 0:N] .* [(-1.0)^j for j in 0:N]
    D = zeros(N + 1, N + 1)
    for i in 0:N, j in 0:N
        i != j && (D[i+1, j+1] = c[i+1] / c[j+1] / (x[i+1] - x[j+1]))
    end
    for i in 0:N
        D[i+1, i+1] = -sum(D[i+1, :])
    end
    M = 1.0; rp = 2.0
    rho = @. (1 / rp) * (1 + x) / 2
    fac = 2 * rp
    D1 = ComplexF64.(fac .* D); D2 = D1 * D1
    Δ = @. 1 - 2 * M * rho; P = @. -(rho^2) * Δ
    A1 = @. 2im * (1 - 8 * rho^2); A0 = @. -2 * (1 + s) * rho + 2 * (s + 3) * rho^2
    B2 = @. -16 * (1 + 2 * rho); B1 = @. 4im * (s - 2 * (s + 2) * rho); B0 = @. 2 * (s + 1) * rho
    C2 = Diagonal(ComplexF64.(B2))
    C1 = Diagonal(ComplexF64.(A1)) * D1 + Diagonal(ComplexF64.(B1))
    C0 = Diagonal(ComplexF64.(P)) * D2 + Diagonal(ComplexF64.(A0)) * D1 + Diagonal(ComplexF64.(B0 .+ lam))
    return Matrix(C0), Matrix(C1), Matrix(C2)
end

@testset "qnm_hyperboloidal (Q2-H) — the third QNM backend" begin

    # ── G-Q2H-1 operator transcription (three independent transcriptions on record) ───────────────
    @testset "G-Q2H-1 operator transcription + Δ̃ typo-guard" begin
        i = im
        # (a) closed forms (#1) == the a=0 anchors (#2), s=±2, at sample R,ω
        for (s, ω, R) in Iterators.product((-2, 2), (0.4 - 0.09im, 1.0 + 0.3im), (0.1, 0.37, 0.5))
            αc = Perturbation._hphc_alpha(s, 2, 0.0, ω, R); βc = Perturbation._hphc_beta(s, 2, 0.0, ω, R)
            if s == -2
                αa = 1 - 4i * ω - 8 * ω^2 - (1 + 16 * ω^2) * R
                βa = i * ω + R + (1 - 8i * ω) * R^2
            else
                αa = -1 + 4i * ω - 8 * ω^2 + (3 - 16i * ω - 16 * ω^2) * R
                βa = i * ω - 3 * R + (5 - 8i * ω) * R^2
            end
            @test abs(αc - αa) ≤ 1e-13 * max(1, abs(αa))
            @test abs(βc - βa) ≤ 1e-13 * max(1, abs(βa))
        end
        # (b) closed forms (#1) == the general-a s=−2 pre-expanded form (#3, Q2H_operator.md L70-74)
        for (a, m, ω, R) in Iterators.product((0.3, 0.7, 0.9), (2, 0, -2), (0.5 - 0.1im,), (0.2, 0.44))
            αc = Perturbation._hphc_alpha(-2, m, a, ω, R)
            βc = Perturbation._hphc_beta(-2, m, a, ω, R)
            α3 = (0.5 * a^2 * ω^2 + a * m * ω - 8 * ω^2 - 4i * ω + 1) +
                 R * (-1 + i * a * m + i * a^2 * ω + 4 * a * m * ω + 4 * a^2 * ω^2 - 16 * ω^2) +
                 R^2 * (-a^2 + 6i * a^2 * ω + 8 * a^2 * ω^2)
            β3 = i * ω + R + R^2 * (1 + i * a * m - 8i * ω + i * a^2 * ω) + R^3 * (-2 * a^2 + 4i * a^2 * ω)
            @test abs(αc - α3) ≤ 1e-12 * max(1, abs(α3))
            @test abs(βc - β3) ≤ 1e-12 * max(1, abs(β3))
        end
        # (c) THE Δ̃ TYPO-GUARD: Δ̃(1/r₊)=0 EXACTLY (1−2MR+a²R²); the printed 1−2MR² does NOT vanish
        for a in (0.0, 0.5, 0.7, 0.9, 0.95)
            rp = 1 + sqrt(1 - a^2); Rp = 1 / rp
            Δ̃ = 1 - 2 * Rp + a^2 * Rp^2                    # correct Δ-form at R=1/r₊
            @test abs(Δ̃) ≤ 1e-14
            typo = 1 - 2 * Rp^2 + a^2 * Rp^2                # the printed supp:44 form
            @test abs(typo) > 1e-2                          # emphatically NOT zero at the horizon
        end
        # (d) DECISIVE: the assembled a=0 ½-form (C0,C1,C2) == ½·cond_bf (== R22/2), element-wise
        for nc in (24, 40)
            op = hphc_radial_operator(-2, 2, 0.0, 0.4 - 0.09im, 4.0; n=nc, T=Float64)
            c0, c1, c2 = _condbf_coeffs(nc + 1; s=-2, lam=4.0)     # cond_bf uses n = nodes = nc+1
            @test maximum(abs, op.C0 .- 0.5 .* c0) ≤ 1e-9 * maximum(abs, c0)
            @test maximum(abs, op.C1 .- 0.5 .* c1) ≤ 1e-11 * maximum(abs, c1)
            @test maximum(abs, op.C2 .- 0.5 .* c2) ≤ 1e-12 * maximum(abs, c2)
        end
        # (e) the a=0 eigenvalue lands on the validated qnm-package QNM (transcription is correct)
        ωbf = qnm_frequency(-2, 2, 2, 0, big"0.0"; backend=:hyperboloidal)
        @test abs(ComplexF64(ωbf) - (0.37367168441804177 - 0.08896231568893546im)) ≤ 1e-11
    end

    # ── G-Q2H-3 MMA cross-oracle (radial eigenfunction, a=0.7 (2,2,0)) ─────────────────────────────
    @testset "G-Q2H-3 MMA radial eigenfunction cross-oracle" begin
        # Julia SWSH λ == BHPT fixture λ (the angular-eigenvalue convention agrees, ~1e-14)
        λj = SpinWeightedSpheroidalEigenvalue(-2, 2, 2, 0.7 * MMA.OMEGA)
        @test abs(λj - MMA.LAMBDA) ≤ 1e-11
        # r₊ / R₊ agree with the fixture
        rf = qnm_radial_hphc(-2, 2, 2, 0, 0.7; n=64, T=Float64)
        @test abs(rf.rp - MMA.RPLUS) ≤ 1e-12
        @test abs(rf.Rp - MMA.RP_COMP) ≤ 1e-12
        # the frequency matches the fixture's Leaver crosscheck value (F64 tier — the companion floor
        # is ~1e-8; the BF tier reaches ~1e-17, gated in G-Q2H-2)
        @test abs(rf.omega - MMA.OMEGA_LEAVER) ≤ 2e-8
        # EIGENFUNCTION shape: build the operator AT THE FIXTURE (ω,λ) (isolates operator structure
        # from the ω-solve) and compare the normalization-invariant scri-ratio R(ρ)/R(scri) to the
        # fixture's ef(ρ)/ef(scri). The ½-form op == ½·QNM.wl ⇒ identical null eigenfunction.
        op = hphc_radial_operator(-2, 2, 0.7, MMA.OMEGA, MMA.LAMBDA; n=100, T=Float64)
        v, _ = Perturbation._hphc_nullvector(op.L)
        Rvals = v ./ v[op.iscri]                       # R(scri)=1 (in-house pin)
        w = Perturbation.bary_weights(op.R)
        ef_scri = MMA.EF[1]                            # fixture ef at scri (= I)
        for (ρ, ef) in zip(MMA.RHO, MMA.EF)
            mine = Perturbation.bary_eval(op.R, w, Rvals, ρ)     # my R(ρ)/R(scri)
            theirs = ef / ef_scri                                 # fixture ef(ρ)/ef(scri)
            @test abs(mine - theirs) ≤ 1e-7 * max(1, abs(theirs))
        end
        # and the scri amplitude ratio directly vs the fixture's I (its own 32-pt resolution ~1e-8)
        @test abs(rf.ratio_scri_hor - MMA.I_SCRI) ≤ 1e-6
    end

    # ── G-Q2H-4 s=+2 Hertz spectrum == s=−2 (Teukolsky s-pairing) ──────────────────────────────────
    @testset "G-Q2H-4 s=+2 == s=−2 spectrum (isospectrality)" begin
        # the s=+2 field has the R^{1+2s}=R^5 peeling ⇒ its collocation needs higher n (spectral):
        # n=48 pairs to ~7e-10, n=80 to ~7e-15. Gate the CONVERGED pairing at n=80.
        for (l, m, a) in ((2, 2, big"0.7"), (3, 3, big"0.5"))
            ωm = Perturbation._qnm_frequency_hyperboloidal(-2, l, m, 0, a; n_coll=80)
            ωp = Perturbation._qnm_frequency_hyperboloidal(2, l, m, 0, a; n_coll=80)
            @test Float64(abs(ωp - ωm)) ≤ 1e-11
            # λ(s=+2) consistency: distinct from λ(s=−2) but both land the SAME ω
            λp = SpinWeightedSpheroidalEigenvalue(2, l, m, a * ωp)
            λm = SpinWeightedSpheroidalEigenvalue(-2, l, m, a * ωm)
            @test abs(λp - λm) > 1e-3                    # genuinely different angular eigenvalues
        end
        # the mirror carry helper (PREP-B L3)
        @test hphc_mirror_carry(0.5 - 0.08im, 2) == (-0.5 - 0.08im, -2)
        @test hphc_mirror_carry(0.3 + 0.1im, -2) == (-0.3 + 0.1im, 2)
    end

    # ── G-Q2H-5 spectral convergence + eltype cleanliness ─────────────────────────────────────────
    @testset "G-Q2H-5 spectral convergence + {T}/BF eltype" begin
        # BF eltype inheritance: operator + frequency are Complex{BigFloat} on a BigFloat input
        op = hphc_radial_operator(-2, 2, big"0.5", Complex{BigFloat}(0.4, -0.09),
                                  Complex{BigFloat}(4); n=24, T=BigFloat)
        @test eltype(op.C0) == Complex{BigFloat} && eltype(op.L) == Complex{BigFloat}
        ωbf = qnm_frequency(-2, 2, 2, 0, big"0.5"; backend=:hyperboloidal)
        @test ωbf isa Complex{BigFloat}
        # MONOTONE spectral decay of |ω(n) − ω_ref| (ω_ref = the leaver BF value), a=0.7 (2,2,0)
        ωref = qnm_frequency(-2, 2, 2, 0, big"0.7"; backend=:leaver,
                             seed=Complex{BigFloat}(qnm_frequency(-2, 2, 2, 0, 0.7; backend=:leaver)))
        errs = Float64[]
        for nc in (24, 32, 48, 64)
            ω = Perturbation._qnm_frequency_hyperboloidal(-2, 2, 2, 0, big"0.7"; n_coll=nc)
            push!(errs, Float64(abs(ω - ωref)))
        end
        @printf("\n[G-Q2H-5] spectral convergence a=0.7 (2,2,0):  n=24:%.2e  n=32:%.2e  n=48:%.2e  n=64:%.2e\n",
                errs...)
        @test errs[2] < errs[1] && errs[3] < errs[2]     # monotone decay 24→32→48
        @test errs[3] ≤ 1e-12                            # converged by n=48
        @test errs[end] ≤ 1e-13
    end

    # ── G-Q2H-2 DECISIVE: the 5×11 campaign, :hyperboloidal == :leaver ≤1e-12 ──────────────────────
    @testset "G-Q2H-2 DECISIVE 5×11 campaign (hyperboloidal == leaver, ≤1e-12)" begin
        χs = vcat(0.0:0.1:0.9, 0.95)                    # 11 spins
        mults = ((2, 2), (3, 3), (4, 4), (2, 0), (2, -2))   # 5 multiplets (n=0, s=-2)
        maxΔ = 0.0
        println("\n[G-Q2H-2] 5×11 campaign  |ω_hphc(BF,n=56) − ω_leaver(BF)| :")
        @printf("  %-8s", "χ")
        for (l, m) in mults; @printf(" (%d,%2d)   ", l, m); end
        println()
        for χ in χs
            @printf("  a=%-5.2f", χ)
            for (l, m) in mults
                aBF = BigFloat(χ)
                ωl64 = qnm_frequency(-2, l, m, 0, Float64(χ); backend=:leaver)
                ωlBF = qnm_frequency(-2, l, m, 0, aBF; backend=:leaver, seed=Complex{BigFloat}(ωl64))
                ωh = Perturbation._qnm_frequency_hyperboloidal(-2, l, m, 0, aBF; n_coll=56)
                Δ = Float64(abs(ωh - ωlBF))
                maxΔ = max(maxΔ, Δ)
                @printf(" %.2e", Δ)
                @test Δ ≤ 1e-12
            end
            println()
        end
        @printf("[G-Q2H-2] max deviation over the 5×11 campaign = %.2e (≤ 1e-12 required)\n", maxΔ)
        @test maxΔ ≤ 1e-12
    end
end
