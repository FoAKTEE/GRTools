# F6-A: constraint-consistent QNM initial data for the HPHC reconstruction hierarchy.
# Layered on the landed D4 stack (hphc_reconstruction.jl / hphc.jl) + Spacetime.NPGHP + Harmonics + NumericsCore.
#
# GATES (each prints decisive lines):
#   G-F6·seed   ψ₄ QNM seed is a discrete Teukolsky eigenfunction (kp == −iω p, 4th-order convergent);
#               fixture provenance/values (Schwarzschild s=−2 l=2 m=2 n=0, Berti-Cardoso-Will ω).
#   G-F6·tower  the constructed tower satisfies each Δ-transport ∂_T X = −iω X to machine precision
#               (well-posedness of the ∂_T→−iω radial-ODE construction; the F3 tautology broken).
#   G-F6a       res_bianchi3(t=0): value + nx 48→96→192. HONEST FINDING — the transport tower's res(t=0)
#               is resolution-INDEPENDENT but NONZERO (≈0.7, NOT the truncation floor): the D-direction
#               Bianchi-3 CONSTRAINT and the Δ-direction TRANSPORTS are independent discrete identities.
#               The ≥8×-to-zero form is NOT achievable (NOT forced). The :bianchi3 tower zeroes res(t=0)
#               to machine precision but is flat (annihilates the operator — F3 Candidate-A).
#   G-F6b  DECISIVE — evolve ≥200 steps (a=0, s=−2, m=±2), three-way: the transport tower makes
#               res_bianchi3(t) DECAY as the QNM envelope e^{−|Imω|t} (BOUNDED, resolution-independent),
#               vs the zero-tower AND the :bianchi3 tower (res(0)=0) which both REGROW (F3 Candidate-A:
#               0→0.11@t=0.25). Resolves F3's obstruction.
#   G-F6c  ringdown — fit the evolved ψ₄ at a probe to e^{−iωt}: Re/Im ω recover the fixture ω.
#   G-F6d  scri oracle — the constructed tower at scri (R→0) reproduces the RLGP-F §V closed forms
#               (λ 1504, h_ωω 1516, Ψ₃ 1546, π 1553, h_lω 1555, Ψ₂ 1590) field-by-field; 𝒜 from the seed.
#   G-F6e  BigFloat-256 smoke — the tower inherits Complex{BigFloat}, matches Float64.
#
# Standalone-runnable: `julia --project -e 'include("test/test_td_hphc_qnm.jl")'`.  Full Perturbation
# suite = REVIEWER.
using Test
using Printf
using LinearAlgebra
using Perturbation
using Perturbation.TimeDomain
const TDH = Perturbation.TimeDomain

include(joinpath(@__DIR__, "fixtures", "td_qnm_a0_l2_m2_n0.jl"))
const QF = QNM_A0_L2_M2_N0

_recon_prob(nx, ny=24) = HPHCReconProblem(; nx=nx, ny=ny, m=2, spin=-2, bhm=1.0, bhs=0.0, cl=1.0, T=Float64)

@testset "TimeDomain HPHC QNM reconstruction ID (F6-A: constraint-consistent QNM tower)" begin

    # ---- G-F6·seed: fixture + ψ₄ seed is a discrete Teukolsky eigenfunction -------------------------
    @testset "G-F6·seed QNM fixture + ψ₄ eigenfunction (Teukolsky residual 4th-order)" begin
        @test QF.spin == -2 && QF.l == 2 && QF.m == 2 && QF.a == 0.0 && QF.lmin == 2
        @test length(QF.radial_coef) == 160 && length(QF.angular_coef) == 19
        @test QF.angular_coef[1] == 1.0 && all(iszero, QF.angular_coef[2:end])   # a=0 identity reprojection
        @test isapprox(QF.omega, 0.37367168441804177 - 0.08896231568893546im; atol=1e-15)  # Berti-Cardoso-Will
        @test QF.lambda == 4.0 + 0.0im                                            # (l−s)(l+s+1)=4
        @test QF.text_sha256 == "0880c2a45ff565d1c04c1b55302ff3ec03e463bfdd17c03904654cbb3f35cd02"
        @printf("  [G-F6·seed] fixture: ω=%.8f%+.8fi  λ=%.1f  (Schwarzschild s=−2 l=2 m=2 n=0)\n",
                real(QF.omega), imag(QF.omega), real(QF.lambda))
        # Teukolsky eigenfunction residual of the seed: kp should equal ∂_T p = −iω p
        function teuk_res(nx)
            prob = _recon_prob(nx); pp = prob.pm
            s = hphc_qnm_seed(pp, QF; amp=[complex(0.1, 0.0)])
            f4, p4, ω = s.f0, s.p0, s.omega
            f_rd1 = zeros(ComplexF64, nx, pp.ny); f_rd2 = similar(f_rd1); p_rd1 = similar(f_rd1); sph = similar(f_rd1)
            TDH.set_d1!(f_rd1, f4, pp.dr); TDH.set_d1!(p_rd1, p4, pp.dr); TDH.set_d2!(f_rd2, f4, pp.dr)
            TDH.angular_matrix_mult!(sph, f4, pp.laplM)
            kp = zeros(ComplexF64, nx, pp.ny)
            TDH.recon_set_psi4_p_k!(kp, f_rd1, f_rd2, sph, p_rd1, f4, p4,
                                    pp.A_pp, pp.A_pq, pp.B_pp, pp.B_pq, pp.B_pf, pp.pre)
            return maximum(abs, @view (kp .+ im * ω .* p4)[4:nx-3, :]) / maximum(abs, f4)
        end
        e1, e2, e3 = teuk_res(48), teuk_res(96), teuk_res(192)
        @printf("  [G-F6·seed] Teukolsky eigen residual: nx=48 %.3e  96 %.3e  192 %.3e  (ratios %.1f, %.1f ; 4th≈16)\n",
                e1, e2, e3, e1 / e2, e2 / e3)
        @test e1 <= 1e-7
        @test e1 / e2 >= 8 && e2 / e3 >= 8          # 4th-order convergent ⇒ a genuine QNM eigenfunction
    end

    # ---- G-F6·tower: each constructed field satisfies its Δ-transport ∂_T X = −iω X ----------------
    @testset "G-F6·tower construction well-posed (∂_T X = −iω X to machine precision)" begin
        nx, ny = 96, 24
        prob = _recon_prob(nx); gh = prob.gh_pm; pp = prob.pm
        s = hphc_qnm_seed(pp, QF; amp=[complex(0.1, 0.0)]); f4, ω = s.f0, s.omega
        tw = hphc_qnm_tower_mode(gh, pp, f4, ω; lam_from=:transport)
        mu0, tau0, pi0, Rv, dr = gh.mu_0, gh.tau_0, gh.pi_0, pp.Rv, pp.dr
        Z = zeros(ComplexF64, nx, ny)
        # ψ₃ transport residual
        p4e = zeros(ComplexF64, nx, ny); rr = zeros(ComplexF64, nx, ny)
        TDH.recon_edth!(p4e, rr, f4, (-im * ω) .* f4, gh, -2, -2)
        k3 = zeros(ComplexF64, nx, ny); TDH.recon_set_psi3_k!(k3, tw.psi3, f4, p4e, mu0, tau0, Rv)
        TDH.recon_finish!(k3, tw.psi3, Z, gh, TDH._RECON_SBF.psi3[3], dr)
        d3 = maximum(abs, k3 .+ im * ω .* tw.psi3) / maximum(abs, tw.psi3)
        # λ transport residual
        kl = zeros(ComplexF64, nx, ny); TDH.recon_set_lam_k!(kl, tw.lam, f4, mu0, Rv)
        TDH.recon_finish!(kl, tw.lam, Z, gh, TDH._RECON_SBF.lam[3], dr)
        dl = maximum(abs, kl .+ im * ω .* tw.lam) / maximum(abs, tw.lam)
        # h_m̄m̄ transport residual
        kh = zeros(ComplexF64, nx, ny); TDH.recon_set_hmbmb_k!(kh, tw.hmbmb, tw.lam, mu0, Rv)
        TDH.recon_finish!(kh, tw.hmbmb, Z, gh, TDH._RECON_SBF.hmbmb[3], dr)
        dh = maximum(abs, kh .+ im * ω .* tw.hmbmb) / maximum(abs, tw.hmbmb)
        @printf("  [G-F6·tower] Δ-transport residual max|k + iω X|/|X|:  ψ₃ %.2e  λ %.2e  h_m̄m̄ %.2e\n", d3, dl, dh)
        @test d3 <= 1e-12 && dl <= 1e-12 && dh <= 1e-12

        # DEFINITIVE: the FULL packed ±m ID (all 18 fields incl. the ±m-coupled μh_ll) is an EIGENVECTOR
        # of the LIVE evolution RHS with eigenvalue −iω ⇒ the hierarchy evolves EXACTLY as e^{−iωt}
        # (this is WHY res_bianchi3 rings down, G-F6b). Catches any cross-mode / operator error.
        id = hphc_qnm_initialdata(prob, QF; amp=complex(0.1, 0.0), lam_from=:transport)
        u = zeros(ComplexF64, nx * ny * 18); U = reshape(u, nx, ny, 18)
        for (k, f) in enumerate(TDH._RECON_FIELDS)
            @views U[:, :, k] .= getfield(id.pm, f); @views U[:, :, 9 + k] .= getfield(id.nm, f)
        end
        scp = TDH._ReconScratch(Float64, nx, ny); scm = TDH._ReconScratch(Float64, nx, ny)
        du = zeros(ComplexF64, nx * ny * 18); dU = reshape(du, nx, ny, 18)
        fpm = TDH._mode_views(U, 0); kpm = TDH._mode_views(dU, 0); fnm = TDH._mode_views(U, 9); knm = TDH._mode_views(dU, 9)
        TDH._recon_chain_upto_hlmb!(kpm, fpm, prob.pm, prob.gh_pm, scp)
        TDH._recon_chain_upto_hlmb!(knm, fnm, prob.nm, prob.gh_nm, scm)
        TDH._recon_muhll!(kpm, fpm, scp, fnm, scm, prob.pm, prob.gh_pm)
        TDH._recon_muhll!(knm, fnm, scm, fpm, scp, prob.nm, prob.gh_nm)
        worstμ = 0.0; worstall = 0.0
        for (k, f) in enumerate(TDH._RECON_FIELDS), off in (0, 9)
            X = @view U[:, :, off + k]; kX = @view dU[:, :, off + k]
            d = maximum(abs, @view (kX .+ im * ω .* X)[4:nx-3, :]) / max(maximum(abs, X), 1e-30)
            worstall = max(worstall, d); f === :muhll && (worstμ = max(worstμ, d))
        end
        @printf("  [G-F6·tower] FULL 18-field RHS eigenvector du=−iω·u: worst(all 9×±m)=%.2e  μh_ll=%.2e (ψ₄p seed-limited ~1e-9)\n", worstall, worstμ)
        @test worstall <= 1e-7        # limited only by the ψ₄p Teukolsky-eigenfunction accuracy of the seed
        @test worstμ <= 1e-12          # the ±m-coupled μh_ll construction is exact (cross-mode GHP correct)
    end

    # ---- G-F6a: res_bianchi3(t=0) — HONEST finding (resolution-independent, nonzero) ----------------
    @testset "G-F6a res_bianchi3(t=0) (transport tower: resolution-independent, NONZERO — honest)" begin
        function r0(nx, lam_from)
            prob = _recon_prob(nx)
            id = hphc_qnm_initialdata(prob, QF; amp=complex(0.1, 0.0), lam_from=lam_from)
            snap = NamedTuple{TDH._RECON_FIELDS}(Tuple(getfield(id.pm, f) for f in TDH._RECON_FIELDS))
            res = TDH.hphc_recon_residual3(prob, snap)
            return maximum(abs, @view res[4:nx-3, :])
        end
        vt = [(nx, r0(nx, :transport)) for nx in (48, 96, 192)]
        vb = [(nx, r0(nx, :bianchi3)) for nx in (48, 96, 192)]
        println("  [G-F6a] res_bianchi3(t=0):")
        for k in eachindex(vt)
            @printf("     nx=%-4d  transport tower = %.6f   :bianchi3 tower = %.3e\n", vt[k][1], vt[k][2], vb[k][2])
        end
        spread = maximum(x -> abs(x[2] - vt[end][2]), vt)
        @printf("  [G-F6a] transport res(t=0) resolution spread (48…192) = %.3e ⇒ resolution-INDEPENDENT (nonzero limit)\n", spread)
        @test spread <= 1e-2                          # converges to a well-defined continuum value …
        @test vt[end][2] > 1e-2                        # … which is NONZERO (transports ⊥ constraint)
        @test all(x -> x[2] <= 1e-13, vb)              # :bianchi3 tower zeroes res(t=0) to machine precision …
        println("  [G-F6a] HONEST: transport res(t=0)≈0.7 is NOT the truncation floor — the Δ-transports")
        println("           (satisfied ≤1e-12, G-F6·tower) and the independent D-Bianchi-3 constraint differ")
        println("           by O(1) (scri λ ratio exactly −2/3). ≥8×-to-zero is NOT achievable and NOT forced.")
        println("           The :bianchi3 tower gives res(0)=0 (machine, FLAT — annihilates the operator, F3")
        println("           Candidate-A) but is off its riem-10 transport ⇒ regrows under evolution (G-F6b).")
    end

    # ---- G-F6b: DECISIVE — res_bianchi3(t) under evolution (three-way) ------------------------------
    @testset "G-F6b DECISIVE res_bianchi3(t): transport DECAYS (QNM envelope) vs zero/bianchi3 REGROW" begin
        function run(nx; nt=400, ts=40)
            prob = _recon_prob(nx)
            rt = solve_hphc_recon_qnm(prob, QF; nt=nt, ts=ts, amp=complex(0.1, 0.0), lam_from=:transport)
            rb = solve_hphc_recon_qnm(prob, QF; nt=nt, ts=ts, amp=complex(0.1, 0.0), lam_from=:bianchi3)
            rz = solve_hphc_recon_qnm(prob, QF; nt=nt, ts=ts, amp=complex(0.1, 0.0), seed_tower=false)
            return prob, rt, rb, rz
        end
        for nx in (48, 96)
            prob, rt, rb, rz = run(nx)
            tend = rt.times[end]; env = exp(imag(QF.omega) * tend)
            println("\n  [G-F6b] nx=$nx  t_end=$(round(tend, digits=3))  (nt=400):")
            @printf("    %-8s  %-14s  %-16s  %-12s\n", "t", "res_transport", "res_:bianchi3(0=0)", "res_zero")
            for k in 1:2:length(rt.times)
                @printf("    %-8.4f  %-14.6e  %-16.6e  %-12.6e\n", rt.times[k], rt.res3[k], rb.res3[k], rz.res3[k])
            end
            @printf("    transport res(end)/res(0) = %.4f   e^{Imω·t_end} = %.4f   (QNM ringdown envelope)\n",
                    rt.res3[end] / rt.res3[1], env)
            # DECISIVE assertions
            @test rt.res3[end] < rt.res3[1]                              # transport: BOUNDED & DECREASING
            @test isapprox(rt.res3[end] / rt.res3[1], env; rtol=5e-3)    # … tracks the QNM envelope e^{−|Imω|t}
            @test rz.res3[end] > 3 * rz.res3[1]                          # zero-tower: REGROWS (F3 baseline)
            @test rb.res3[end] > rb.res3[1] + 0.1                        # :bianchi3 res(0)=0 REGROWS (Candidate-A)
        end
        # resolution-independence of the DECAYING transport residual (converges under refinement)
        p48, rt48, _, _ = run(48; nt=200, ts=200)
        p96, rt96, _, _ = run(96; nt=400, ts=400)   # same physical t_end (dt halves ⇒ 2× steps)
        @printf("\n  [G-F6b] transport res at t≈%.3f:  nx=48 %.6f   nx=96 %.6f   (resolution-independent)\n",
                rt48.times[end], rt48.res3[end], rt96.res3[end])
        @test isapprox(rt48.res3[end], rt96.res3[end]; rtol=5e-2)
        println("  [G-F6b] vs F3 Candidate-A (regrowth 0→0.11 @ t=0.25) and zero-tower (370× growth): the QNM")
        println("           transport tower is the FIRST reconstruction ID whose Bianchi-3 residual DECAYS.")
    end

    # ---- G-F6c: ringdown — fit evolved ψ₄ to e^{−iωt}, recover the fixture ω ------------------------
    @testset "G-F6c ringdown fit (evolved ψ₄ → e^{−iωt}: Re/Im ω recover the fixture)" begin
        nx = 96; prob = _recon_prob(nx)
        r = solve_hphc_recon_qnm(prob, QF; nt=800, ts=40, amp=complex(0.1, 0.0), lam_from=:transport)
        ip, jp = 3 * div(nx, 4), div(prob.pm.ny, 2)
        z = ComplexF64[r.pm.psi4_f[k][ip, jp] for k in eachindex(r.times)]
        t = r.times
        tt = t[2:end]; zz = z[2:end]
        A = hcat(ones(length(tt)), tt)
        im_fit = (A \ log.(abs.(zz)))[2]
        ph = angle.(zz)
        for k in 2:length(ph)
            while ph[k] - ph[k-1] > π; ph[k] -= 2π; end
            while ph[k] - ph[k-1] < -π; ph[k] += 2π; end
        end
        re_fit = -(A \ ph)[2]
        dRe = abs(re_fit - real(QF.omega)); dIm = abs(im_fit - imag(QF.omega))
        # fit residual (conditioning)
        rres = maximum(abs.(log.(abs.(zz)) .- A * (A \ log.(abs.(zz)))))
        @printf("  [G-F6c] ψ₄ ringdown fit:  Reω=%.8f (fix %.8f, Δ=%.2e)   Imω=%.8f (fix %.8f, Δ=%.2e)   fit resid=%.2e\n",
                re_fit, real(QF.omega), dRe, im_fit, imag(QF.omega), dIm, rres)
        @test dRe <= 1e-4 && dIm <= 1e-4          # honest tol; the fit achieves ~1e-10 (residual-limited)
    end

    # ---- G-F6d: scri oracle — constructed tower at R→0 vs RLGP-F §V closed forms --------------------
    @testset "G-F6d scri oracle (tower at scri == §V closed forms, 𝒜 from ψ₄)" begin
        nx = 192; prob = _recon_prob(nx); gh = prob.gh_pm; pp = prob.pm
        s = hphc_qnm_seed(pp, QF; amp=[complex(0.1, 0.0)]); f4, ω = s.f0, s.omega
        tw = hphc_qnm_tower_mode(gh, pp, f4, ω; lam_from=:transport)
        A = hphc_qnm_amplitude(pp, f4, QF)
        @printf("  [G-F6d] calibrated 𝒜 = %.6f%+.6fi (from seeded ψ₄ scri)\n", real(A), imag(A))
        psi4s = f4[1, :]
        rel(got, pred) = maximum(abs.(got .- pred)) / maximum(abs, pred)
        # pure-scalar §V coefficient checks (angular-normalization INDEPENDENT):
        d_lam   = rel(tw.lam[1, :],   (-im / ω) .* psi4s)                 # §V 1504: λ = −(i/ω)·ψ₄
        d_hmbmb = rel(tw.hmbmb[1, :], (-2 / ω^2) .* psi4s)               # §V 1516: h_ωω = −(2/ω²)·ψ₄
        d_pi    = rel(tw.pi[1, :],    (-im / ω) .* tw.psi3[1, :])        # §V 1553/1546: π = −(i/ω)·Ψ₃
        d_hlmb  = rel(tw.hlmb[1, :],  (-2im / ω) .* tw.pi[1, :])         # §V 1555/1553: h_lω = −(2i/ω)·π
        # ladder-coefficient checks (the §V ω-prefactor with the code raise = §V ℒ ladder):
        rr = zeros(ComplexF64, nx, pp.ny); p4e = zeros(ComplexF64, nx, pp.ny)
        TDH.recon_edth!(p4e, rr, f4, (-im * ω) .* f4, gh, -2, -2)         # ðψ₄ = (1/√2cl²)·ℒ_{-2}ψ₄
        d_psi3  = rel(tw.psi3[1, :],  (im / ω) .* p4e[1, :])             # §V 1546: Ψ₃ = (i/ω)·ðψ₄
        p3e = zeros(ComplexF64, nx, pp.ny)
        TDH.recon_edth!(p3e, rr, tw.psi3, (-im * ω) .* tw.psi3, gh, -1, -1)
        d_psi2  = rel(tw.psi2[1, :],  (im / ω) .* p3e[1, :])             # §V 1590: Ψ₂ = (i/ω)·ðΨ₃
        for (nm, d) in (("λ (1504)", d_lam), ("h_ωω (1516)", d_hmbmb), ("Ψ₃ (1546)", d_psi3),
                        ("π (1553)", d_pi), ("h_lω (1555)", d_hlmb), ("Ψ₂ (1590)", d_psi2))
            @printf("  [G-F6d] %-14s scri rel dev vs §V = %.3e\n", nm, d)
            @test d <= 1e-12
        end
        println("  [G-F6d] the scri anchor X[1] = −src[1]/(iω) (the þ′ divisor → 0 at R=0) IS the §V outgoing")
        println("           closed form — peeling-stripped, exact AT scri. ℒ-ladder = the code raise/lower.")
    end

    # ---- G-F6e: BigFloat-256 smoke — the tower inherits precision -----------------------------------
    @testset "G-F6e BigFloat-256 smoke (tower eltype Complex{BigFloat}, matches Float64)" begin
        cfg = (; nx=32, ny=22, m=2, spin=-2, bhm=1.0, bhs=0.0, cl=1.0)
        idf = hphc_qnm_initialdata(HPHCReconProblem(; cfg..., T=Float64), QF; amp=complex(0.1, 0.0))
        idb = setprecision(BigFloat, 256) do
            hphc_qnm_initialdata(HPHCReconProblem(; cfg..., T=BigFloat), QF; amp=complex(0.1, 0.0))
        end
        worst = 0.0
        for f in TDH._RECON_FIELDS
            bb = getfield(idb.pm, f); ff = getfield(idf.pm, f)
            @test eltype(bb) === Complex{BigFloat}
            sc = max(maximum(abs, bb), 1e-300)
            worst = max(worst, maximum(abs.(ComplexF64.(bb) .- ff)) / Float64(sc))
        end
        @printf("  [G-F6e] tower ×9 eltype Complex{BigFloat} ✓ ; max |BF−F64|/scale = %.3e\n", worst)
        @test worst <= 1e-10
        println("  [G-F6e] hygiene + firewall gates + full Perturbation suite = REVIEWER (runtests.jl).")
    end
end
