# F7-1: the second-order Teukolsky source (𝔰_d, 𝔰_t, 𝒮) bilinear assembler.
# Layered on the landed D4 stack (hphc_reconstruction.jl / hphc_qnm_id.jl) + Spacetime.NPGHP + Harmonics.
#
# SOURCE OF TRUTH: arXiv:2010.00162 numerics_description.tex:924-981 (eq:schematic_source +
# eq:def_mathfrak_s), transcription DOUBLE-ENTRY-verified here; the owner's H_± channel/support
# formalism teukolsky-category.tex:1271-1341 + coeff/spin2-minus-minus-*.generated.tex (F7-α structural).
#
# There is NO time-domain 𝒮-VALUE oracle (upstream TeukEvolution.jl only DECLARES scd_f/scd_p and never
# assembles the source; the owner's tables are frequency-domain separated, a different basis). So this
# round gates the assembler STRUCTURALLY (transcription double-entry, H_± grading / Π_{--}=0, bilinearity,
# QNM smoke + time-scaling) — the F7-α program. Ψ₄⁽²⁾-solve value validation vs the QQNM amplitudes is
# F7-2 [FUTURE].
#
# GATES (each prints decisive lines):
#   G-F7-1a  per-term transcription DOUBLE-ENTRY: scd_sd!/scd_st! == an INDEPENDENT tex re-transcription
#            (bracket-grouped) on seeded-random bundle data; assembled 𝒮 == d₄𝔰_d+d₃𝔰_t; and the source
#            GHP applicators === the locked recon_* applicators where the ladders coincide.
#   G-F7-1b  STRUCTURAL H_± grading (owner's slot proofs): supp_{H±}(𝒮)⊂{(2,0),(1,1)}; Π_{--}𝒮=0 exactly;
#            per-term bidegree table (𝔰_d all-(2,0) except the h_lm term; 𝔰_t all-(1,1)); no H_- in the
#            Weyl slot. Realized via the ±m structure (H_+ = +m tower content, H_- = the mirror-conj).
#   G-F7-1c  bilinearity 𝒮[c₁B₁,c₂B₂]=c₁c₂𝒮[B₁,B₂] + the exchange/one-mixed pin ([AXIOM-F7]).
#   G-F7-1d  QNM smoke on the F6-A (2,2) tower: finite; nonzero in the m=4 & m=0 sectors; Π_{--}=0;
#            𝒮 scales as e^{−2iωt} EXACTLY (two snapshot assemblies).
#   G-F7-1e  BigFloat-256 eltype-clean (Complex{BigFloat}, matches Float64).
#
# Standalone-runnable: `julia --project -e 'include("test/test_td_hphc_scd.jl")'`. Full Perturbation
# suite (+ hygiene + firewall) = REVIEWER (runtests.jl auto-includes this file).
using Test
using Printf
using LinearAlgebra
# (no Random: the Pkg.test sandbox carries only [targets] deps — seeded LCG per campaign convention)
using Perturbation
using Perturbation.TimeDomain
const TDH = Perturbation.TimeDomain

include(joinpath(@__DIR__, "fixtures", "td_qnm_a0_l2_m2_n0.jl"))
const QF = QNM_A0_L2_M2_N0

_prob(nx, ny=24; T=Float64) = HPHCReconProblem(; nx=nx, ny=ny, m=2, spin=-2, bhm=1.0, bhs=0.0, cl=1.0, T=T)

# ---- test helpers ---------------------------------------------------------------------------------

# a seeded pseudo-random bundle (complex nx×ny fields) at (ω, m_ang) — for the algebra/grading gates.
# Deterministic literal-int LCG (campaign convention; the gates are same-data identities, so the
# stream itself is arbitrary — no baked fixture depends on it, and Random is not a test target).
function rand_bundle(::Type{T}, nx, ny, ω, m_ang; seed=1) where {T}
    s = UInt64(0x9e3779b97f4a7c15) ⊻ (UInt64(seed) * 0x2545f4914f6cdd1d)
    nextu() = (s = s * 0x5851f42d4c957f2d + 0x14057b7ef767814f; s)
    # uniform in (-1, 1), full 52-bit mantissa entropy
    nextr() = (T(reinterpret(Float64, (nextu() >> 12) | 0x3ff0000000000000) - 1) * 2 - 1)
    F() = Complex{T}[Complex{T}(nextr(), nextr()) for _ in 1:nx, _ in 1:ny]
    SCDBundle{T}(F(), F(), F(), F(), F(), F(), F(), F(), F(), F(), F(), Complex{T}(ω), Int(m_ang))
end

scale_bundle(B::SCDBundle{T}, c) where {T} =
    SCDBundle{T}((c .* getfield(B, f) for f in TDH._SCD_BUNDLE_FIELDS)..., B.omega, B.m_ang)

# split a bundle into its H_+ part (Weyl + H_+ metric; H_- = 0) and H_- part (only h_lm,h_mm,π̄).
function split_pm(B::SCDBundle{T}) where {T}
    Z = zeros(Complex{T}, size(B.psi4))
    Hp = SCDBundle{T}(B.psi4, B.psi3, B.psi2, B.hll, B.hlmb, B.hmbmb, B.lam, B.pi,
                      copy(Z), copy(Z), copy(Z), B.omega, B.m_ang)
    Hm = SCDBundle{T}(copy(Z), copy(Z), copy(Z), copy(Z), copy(Z), copy(Z), copy(Z), copy(Z),
                      B.hlm, B.hmm, B.pibar, B.omega, B.m_ang)
    return Hp, Hm
end

relerr(A, B) = maximum(abs, A .- B) / max(maximum(abs, B), eps(real(eltype(B))))

# INDEPENDENT re-transcription of eq:def_mathfrak_s (bracket-grouped, distinct accumulation order from
# the src module's flat term list) — the second reading of the double-entry (G-F7-1a).
function sd_reference(B1::SCDBundle{T}, B2::SCDBundle{T}, ops) where {T}
    nx, ny = size(B1.psi4); a1, a2 = B1.m_ang, B2.m_ang; ω1, ω2 = B1.omega, B2.omega
    mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; tau0 = ops.gh.tau_0
    cpi = conj.(pi0); ctau = conj.(tau0)
    W = () -> zeros(Complex{T}, nx, ny); DR = W(); rr = W(); lw = W(); o = W()
    # operator shorthands (DT = −iω·field, the QNM fact)
    Δ(f, ω, fall)    = (oo = W(); TDH.scd_thornprime!(oo, f, (-im * ω) .* f, DR, ops, fall); oo)
    ð(f, ω, s, b, m) = (oo = W(); TDH.scd_edth!(oo, rr, f, (-im * ω) .* f, ops, s, b, m); oo)
    ð′(f, ω, s, b, m)= (oo = W(); TDH.scd_edth_prime!(oo, lw, f, (-im * ω) .* f, ops, s, b, m); oo)
    Ψ4, Ψ3, Ψ2 = B1.psi4, B1.psi3, B1.psi2
    hll, hlmb, hmbmb, hlm, λ, π1 = B2.hll, B2.hlmb, B2.hmbmb, B2.hlm, B2.lam, B2.pi
    bracket1 = (1 // 2) .* (ð(hlmb, ω2, -1, 1, a2) .+ (cpi .+ 2 .* tau0) .* hlmb) .+
               (Δ(hll, ω2, 3) .+ (-mu .+ cmu) .* hll) .-
               (1 // 2) .* (ð′(hlm, ω2, 1, 1, a2) .+ (-5 .* pi0 .- 4 .* ctau) .* hlm)
    bracket2 = (ð(hmbmb, ω2, -2, 0, a2) .+ (cpi .+ tau0) .* hmbmb) .+
               (Δ(hlmb, ω2, 2) .+ (-2 .* mu .+ cmu) .* hlmb)
    bracket3 = hlmb .* Δ(Ψ3, ω1, 2) .- (1 // 2) .* hmbmb .* ð(Ψ3, ω1, -1, -1, a1) .- 4 .* π1 .* Ψ3
    return (1 // 2) .* hll .* (Δ(Ψ4, ω1, 1) .+ mu .* Ψ4) .+
           Ψ4 .* bracket1 .- (1 // 2) .* Ψ3 .* bracket2 .- bracket3 .- 3 .* λ .* Ψ2
end

function st_reference(B1::SCDBundle{T}, B2::SCDBundle{T}, ops) where {T}
    nx, ny = size(B1.psi4); a1, a2 = B1.m_ang, B2.m_ang; ω1, ω2 = B1.omega, B2.omega
    mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; tau0 = ops.gh.tau_0; ctau = conj.(tau0)
    W = () -> zeros(Complex{T}, nx, ny); DR = W(); lw = W()
    Δ(f, ω, fall)     = (oo = W(); TDH.scd_thornprime!(oo, f, (-im * ω) .* f, DR, ops, fall); oo)
    ð′(f, ω, s, b, m) = (oo = W(); TDH.scd_edth_prime!(oo, lw, f, (-im * ω) .* f, ops, s, b, m); oo)
    Ψ4 = B1.psi4; hlm, hmm, pibar = B2.hlm, B2.hmm, B2.pibar
    return .-hlm .* (Δ(Ψ4, ω1, 1) .+ (mu .+ 2 .* cmu) .* Ψ4) .+
           (1 // 2) .* hmm .* ð′(Ψ4, ω1, -2, -2, a1) .+
           Ψ4 .* (pibar .- Δ(hlm, ω2, 2) .+ (ð′(hmm, ω2, 2, 0, a2) .+ (-(1 // 2) .* pi0 .- (1 // 2) .* ctau) .* hmm))
end

@testset "TimeDomain HPHC 2nd-order source 𝔰_d/𝔰_t/𝒮 (F7-1: bilinear assembler, F7-α structural)" begin

    # ---- G-F7-1a: transcription DOUBLE-ENTRY + applicator grounding --------------------------------
    @testset "G-F7-1a per-term transcription double-entry (vs independent tex reading)" begin
        nx, ny = 32, 24
        prob = _prob(nx, ny); ops = SCDSourceOps(prob)
        # SEEDED random bundles (two distinct seeds) at ω from the fixture, azimuthal (+2, +2).
        ω = Complex{Float64}(QF.omega)
        B1 = rand_bundle(Float64, nx, ny, ω, 2; seed=11)
        B2 = rand_bundle(Float64, nx, ny, ω, 2; seed=23)
        sc = TDH._SCDScratch(Float64, nx, ny)
        sd = zeros(ComplexF64, nx, ny); st = zeros(ComplexF64, nx, ny)
        TDH.scd_sd!(sd, B1, B2, ops, sc); TDH.scd_st!(st, B1, B2, ops, sc)
        sd_ref = sd_reference(B1, B2, ops); st_ref = st_reference(B1, B2, ops)
        e_sd = relerr(sd, sd_ref); e_st = relerr(st, st_ref)
        @printf("  [G-F7-1a] 𝔰_d src-vs-independent-tex rel dev = %.3e   𝔰_t = %.3e  (double-entry)\n", e_sd, e_st)
        @test e_sd <= 1e-12
        @test e_st <= 1e-12

        # assembled 𝒮 == d₄⁽⁰⁾𝔰_d + d₃⁽⁰⁾𝔰_t recomputed independently from sd_ref/st_ref
        S = zeros(ComplexF64, nx, ny)
        r = scd_source!(S, B1, B2, ops)
        mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; ctau = conj.(ops.gh.tau_0)
        ωtot = 2ω; DR = zeros(ComplexF64, nx, ny); lw = zeros(ComplexF64, nx, ny)
        d4sd = zeros(ComplexF64, nx, ny); TDH.scd_thornprime!(d4sd, sd_ref, (-im * ωtot) .* sd_ref, DR, ops, 1)
        d4sd .+= (4 .* mu .+ cmu) .* sd_ref
        d3st = zeros(ComplexF64, nx, ny); TDH.scd_edth_prime!(d3st, lw, st_ref, (-im * ωtot) .* st_ref, ops, -1, -2, 4)
        d3st .+= (4 .* pi0 .- ctau) .* st_ref
        e_S = relerr(S, d4sd .+ d3st)
        @printf("  [G-F7-1a] assembled 𝒮 == d₄𝔰_d+d₃𝔰_t rel dev = %.3e\n", e_S)
        @test e_S <= 1e-12
        @test relerr(r.sd, sd) == 0.0 && relerr(r.st, st) == 0.0   # scd_source! returns the same 𝔰_d,𝔰_t

        # applicator grounding: scd_edth!(gh.raise[-1]) === the locked recon_edth! (bit-exact)
        f = B1.psi3; DTf = (-im * ω) .* f
        o1 = zeros(ComplexF64, nx, ny); rr1 = zeros(ComplexF64, nx, ny)
        o2 = zeros(ComplexF64, nx, ny); rr2 = zeros(ComplexF64, nx, ny)
        TDH.scd_edth!(o1, rr1, f, DTf, ops, -1, -1, 2)
        TDH.recon_edth!(o2, rr2, f, DTf, prob.gh_pm, -1, -1)
        de = maximum(abs, o1 .- o2)
        o3 = zeros(ComplexF64, nx, ny); lw3 = zeros(ComplexF64, nx, ny)
        o4 = zeros(ComplexF64, nx, ny); lw4 = zeros(ComplexF64, nx, ny)
        TDH.scd_edth_prime!(o3, lw3, f, DTf, ops, -1, -1, 2)
        TDH.recon_edth_prime!(o4, lw4, f, DTf, prob.gh_pm, -1, -1)
        dep = maximum(abs, o3 .- o4)
        o5 = zeros(ComplexF64, nx, ny); DR5 = zeros(ComplexF64, nx, ny)
        o6 = zeros(ComplexF64, nx, ny); DR6 = zeros(ComplexF64, nx, ny)
        TDH.scd_thornprime!(o5, f, DTf, DR5, ops, 2)               # computes ∂_R f internally (set_d1!)
        TDH.set_d1!(DR6, f, prob.pm.dr)                            # recon_thorn_prime! takes ∂_R f as arg
        TDH.recon_thorn_prime!(o6, f, DTf, DR6, prob.gh_pm, 2)
        dtp = maximum(abs, o5 .- o6)
        @printf("  [G-F7-1a] src applicators === locked recon_*: ð %.1e  ð′ %.1e  Δ %.1e (bit-exact)\n", de, dep, dtp)
        @test de == 0.0 && dep == 0.0 && dtp == 0.0
        # the pre_* coefficient fields ARE m-independent (why one gh serves all sectors)
        dm = maximum(abs, prob.gh_pm.pre_thorn_prime_DT .- prob.gh_nm.pre_thorn_prime_DT) +
             maximum(abs, prob.gh_pm.pre_edth_prime_lowered .- prob.gh_nm.pre_edth_prime_lowered)
        @printf("  [G-F7-1a] pre_* fields pm≡nm (m-independent): max diff = %.1e\n", dm)
        @test dm == 0.0
    end

    # ---- G-F7-1b: STRUCTURAL H_± grading (owner's slot proofs) -------------------------------------
    @testset "G-F7-1b H_± grading: supp(𝒮)⊂{(2,0),(1,1)}, Π_{--}𝒮=0 (owner slot proofs)" begin
        nx, ny = 32, 24
        prob = _prob(nx, ny); ops = SCDSourceOps(prob)
        ω = Complex{Float64}(QF.omega)
        B = rand_bundle(Float64, nx, ny, ω, 2; seed=7)
        Hp, Hm = split_pm(B)                                # Hp = H_+ content, Hm = H_- content
        S(b1, b2) = (SS = zeros(ComplexF64, nx, ny); scd_source!(SS, b1, b2, ops); SS)
        S_20 = S(Hp, Hp)       # (2,0): H_+ Weyl × H_+ metric
        S_11 = S(Hp, Hm)       # (1,1): H_+ Weyl × H_- metric
        S_02 = S(Hm, Hm)       # (0,2): H_- Weyl × H_- metric  →  Π_{--}
        S_m2 = S(Hm, Hp)       # H_- in the Weyl slot          →  must vanish
        Sfull = S(B, B)
        n20 = maximum(abs, S_20); n11 = maximum(abs, S_11)
        @printf("  [G-F7-1b] channel magnitudes:  (2,0)=%.4e  (1,1)=%.4e  (0,2)=%.2e  H_-→Weyl=%.2e\n",
                n20, n11, maximum(abs, S_02), maximum(abs, S_m2))
        @test n20 > 1e-6 && n11 > 1e-6                       # both admitted channels are populated …
        @test maximum(abs, S_02) == 0.0                      # … (0,2)=H_-² is IDENTICALLY zero  ⇒ Π_{--}𝒮=0
        @test maximum(abs, S_m2) == 0.0                      # H_- never enters the Weyl (source) slot
        # supp(𝒮) ⊂ {(2,0),(1,1)}: the full source is EXACTLY the sum of the two admitted channels
        e_supp = relerr(Sfull, S_20 .+ S_11)
        @printf("  [G-F7-1b] 𝒮[B,B] == 𝒮_(2,0)+𝒮_(1,1)  rel dev = %.3e  ⇒ supp⊂{(2,0),(1,1)}\n", e_supp)
        @test e_supp <= 1e-12

        # per-term bidegree table (𝔰_d, 𝔰_t) — realize the grading via the ±m structure and compare to
        # the owner's eight-term table (aggregate; organizations differ term-by-term, stated honestly).
        sc = TDH._SCDScratch(Float64, nx, ny)
        function bidegrees(assemble!)
            grid = Dict{Tuple{Int,Int},Float64}()
            for (p, (b1, b2)) in ((( 2,0),(Hp,Hp)), ((1,1),(Hp,Hm)), ((1,1),(Hm,Hp)), ((0,2),(Hm,Hm)))
                m = zeros(ComplexF64, nx, ny); assemble!(m, b1, b2, ops, sc)
                grid[p] = max(get(grid, p, 0.0), maximum(abs, m))
            end
            grid
        end
        gd = bidegrees(TDH.scd_sd!); gt = bidegrees(TDH.scd_st!)
        @printf("  [G-F7-1b] 𝔰_d bidegree mag:  (2,0)=%.3e  (1,1)=%.3e  (0,2)=%.1e\n", gd[(2,0)], gd[(1,1)], gd[(0,2)])
        @printf("  [G-F7-1b] 𝔰_t bidegree mag:  (2,0)=%.1e  (1,1)=%.3e  (0,2)=%.1e\n", gt[(2,0)], gt[(1,1)], gt[(0,2)])
        @test gd[(2,0)] > 1e-6 && gd[(1,1)] > 1e-6 && gd[(0,2)] == 0.0   # 𝔰_d: (2,0)+(1,1), no H_-²
        @test gt[(2,0)] == 0.0 && gt[(1,1)] > 1e-6 && gt[(0,2)] == 0.0   # 𝔰_t: pure (1,1)
        println("  [G-F7-1b] owner eight-term table (teukolsky-category.tex:1332-1340): supp(T₁₋₄,₇)⊂{(2,0),(1,1)},")
        println("           supp(T₅)⊂{(2,0)}, T₆=T₈=0 ⇒ supp(𝒮)⊂{(2,0),(1,1)}, Π_{--}𝒮=0. Ripley's 𝔰_d/𝔰_t is a")
        println("           DIFFERENT organization of the same 𝒮 (not a term-by-term T_i map); the AGGREGATE H_±")
        println("           support + Π_{--}=0 match exactly (basis-independent). H_+ = +m tower, H_- = mirror-conj.")
    end

    # ---- G-F7-1c: bilinearity + exchange / one-mixed counting pin ----------------------------------
    @testset "G-F7-1c bilinearity 𝒮[c₁B₁,c₂B₂]=c₁c₂𝒮[B₁,B₂] + one-mixed pin [AXIOM-F7]" begin
        nx, ny = 32, 24
        prob = _prob(nx, ny); ops = SCDSourceOps(prob)
        ω = Complex{Float64}(QF.omega)
        B1 = rand_bundle(Float64, nx, ny, ω, 2; seed=3)
        B2 = rand_bundle(Float64, nx, ny, ω, 2; seed=5)
        c1 = 0.37 - 1.4im; c2 = -2.1 + 0.6im
        S0 = zeros(ComplexF64, nx, ny); scd_source!(S0, B1, B2, ops)
        Sc = zeros(ComplexF64, nx, ny); scd_source!(Sc, scale_bundle(B1, c1), scale_bundle(B2, c2), ops)
        e_bil = relerr(Sc, (c1 * c2) .* S0)
        @printf("  [G-F7-1c] bilinearity rel dev = %.3e   (𝒮[c₁B₁,c₂B₂] = c₁c₂ 𝒮[B₁,B₂])\n", e_bil)
        @test e_bil <= 1e-12

        # additivity in each slot (a second, independent linearity witness)
        B2b = rand_bundle(Float64, nx, ny, ω, 2; seed=9)
        Ssum = zeros(ComplexF64, nx, ny); scd_source!(Ssum, B1, SCDBundle{Float64}(
            (getfield(B2, f) .+ getfield(B2b, f) for f in TDH._SCD_BUNDLE_FIELDS)..., ω, 2), ops)
        Sa = zeros(ComplexF64, nx, ny); scd_source!(Sa, B1, B2, ops)
        Sb = zeros(ComplexF64, nx, ny); scd_source!(Sb, B1, B2b, ops)
        e_add = relerr(Ssum, Sa .+ Sb)
        @printf("  [G-F7-1c] additivity in slot-2 rel dev = %.3e\n", e_add)
        @test e_add <= 1e-12

        # one-mixed counting pin: H_- lives ONLY in the metric slot ⇒ 𝒮[Hp,Hm]≠0 but 𝒮[Hm,Hp]=0, so the
        # mixed channel is single (NOT doubled) — the owner's R_{+-} counting (comparison/2410-22...md). [AXIOM-F7]
        Hp, Hm = split_pm(rand_bundle(Float64, nx, ny, ω, 2; seed=13))
        Spm = zeros(ComplexF64, nx, ny); scd_source!(Spm, Hp, Hm, ops)
        Smp = zeros(ComplexF64, nx, ny); scd_source!(Smp, Hm, Hp, ops)
        @printf("  [G-F7-1c] mixed one-sided: |𝒮[H_+,H_-]|=%.4e  |𝒮[H_-,H_+]|=%.1e (single mixed ⇒ un-doubled)\n",
                maximum(abs, Spm), maximum(abs, Smp))
        @test maximum(abs, Spm) > 1e-6 && maximum(abs, Smp) == 0.0
        println("  [AXIOM-F7] scd_source!(S,B1,B2): Weyl/source factor←B1, metric/operator factor←B2 (ordered).")
        println("             Diagonal physical source = 𝒮[B,B]; symmetric bilinear ½(𝒮[B1,B2]+𝒮[B2,B1]).")
        println("             H_- confined to the metric slot ⇒ 𝒮[H_-,·]≡0, mixed channel single (R_{+-} pin).")
    end

    # ---- G-F7-1d: QNM smoke on the F6-A (2,2) tower ------------------------------------------------
    @testset "G-F7-1d QNM smoke (F6-A (2,2) tower: finite, m=4/m=0 nonzero, Π_{--}=0, 𝒮∝e^{−2iωt})" begin
        nx, ny = 48, 24
        prob = _prob(nx, ny); ops = SCDSourceOps(prob)
        id = hphc_qnm_initialdata(prob, QF; amp=complex(0.1, 0.0), lam_from=:transport)
        r4 = scd_source_sector(id, id, 2, :pp, ops)   # m=4 sector
        r0 = scd_source_sector(id, id, 2, :pm, ops)   # m=0 sector
        @printf("  [G-F7-1d] sector m=%d: max|𝒮|=%.4e (finite=%s)   sector m=%d: max|𝒮|=%.4e (finite=%s)\n",
                r4.m_ang, maximum(abs, r4.S), all(isfinite, r4.S), r0.m_ang, maximum(abs, r0.S), all(isfinite, r0.S))
        @test all(isfinite, r4.S) && all(isfinite, r0.S)
        @test r4.m_ang == 4 && r0.m_ang == 0
        @test maximum(abs, r4.S) > 1e-6 && maximum(abs, r0.S) > 1e-6

        # Π_{--} on the physical (2,2) diagonal source (m=4 sector): split the +m bundle H_+/H_-.
        Bp = scd_bundle(id, :pm, 2, ops)
        Hp, Hm = split_pm(Bp)
        S02 = zeros(ComplexF64, nx, ny); scd_source!(S02, Hm, Hm, ops)
        S20 = zeros(ComplexF64, nx, ny); scd_source!(S20, Hp, Hp, ops)
        S11 = zeros(ComplexF64, nx, ny); scd_source!(S11, Hp, Hm, ops)
        @printf("  [G-F7-1d] m=4 physical: Π_{--}=|𝒮_(0,2)|=%.1e ; 𝒮[B,B]==𝒮_(2,0)+𝒮_(1,1) rel dev=%.3e\n",
                maximum(abs, S02), relerr(r4.S, S20 .+ S11))
        @test maximum(abs, S02) == 0.0                                  # Π_{--}𝒮=0 on real QNM data
        @test relerr(r4.S, S20 .+ S11) <= 1e-12

        # 𝒮 ∝ e^{−2iωt} EXACTLY: assemble at two snapshots (all bundle fields ×e^{−iωt}, coherent-ω).
        ω = Complex{Float64}(id.omega); tsnap = 0.83
        c = exp(-im * ω * tsnap)
        B1t = scale_bundle(Bp, c); B2t = scale_bundle(Bp, c)
        St = zeros(ComplexF64, nx, ny); scd_source!(St, B1t, B2t, ops)
        pred = (c^2) .* r4.S
        e_scale = relerr(St, pred)
        @printf("  [G-F7-1d] two-snapshot time scaling: 𝒮(t)=e^{−2iωt}𝒮(0) rel dev = %.3e  (e^{−2iωt}=%.4f%+.4fi)\n",
                e_scale, real(c^2), imag(c^2))
        @test e_scale <= 1e-10
        println("  [G-F7-1d] every source factor ∝ e^{−iωt} (DT=−iωX, D4 fact) ⇒ the bilinear 𝒮 ∝ e^{−2iωt}.")
    end

    # ---- G-F7-1e: BigFloat-256 eltype-clean -------------------------------------------------------
    @testset "G-F7-1e BigFloat-256 smoke (𝒮 eltype Complex{BigFloat}, matches Float64)" begin
        nx, ny = 32, 22
        idf = hphc_qnm_initialdata(_prob(nx, ny; T=Float64), QF; amp=complex(0.1, 0.0))
        opsf = SCDSourceOps(_prob(nx, ny; T=Float64))
        rf = scd_source_sector(idf, idf, 2, :pp, opsf)
        rb = setprecision(BigFloat, 256) do
            probb = _prob(nx, ny; T=BigFloat)
            idb = hphc_qnm_initialdata(probb, QF; amp=complex(0.1, 0.0))
            opsb = SCDSourceOps(probb)
            scd_source_sector(idb, idb, 2, :pp, opsb)
        end
        @test eltype(rb.S) === Complex{BigFloat}
        @test eltype(rb.sd) === Complex{BigFloat} && eltype(rb.st) === Complex{BigFloat}
        d = maximum(abs.(ComplexF64.(rb.S) .- rf.S)) / max(maximum(abs, rf.S), 1e-300)
        @printf("  [G-F7-1e] 𝒮 eltype Complex{BigFloat} ✓ ; max|BF−F64|/scale = %.3e\n", d)
        @test d <= 1e-8
        # a≠0 axiom guard: SCDSourceOps must REFUSE a Kerr problem (μ₀ non-constant ⇒ μ₀→0 h_ll
        # unpacking is [FUTURE]). Build the Kerr problem OUTSIDE the try (isolate the guard's throw).
        kerr = HPHCReconProblem(; nx=16, ny=24, m=2, spin=-2, bhs=0.5, T=Float64)
        threw = false
        try
            SCDSourceOps(kerr)
        catch e
            threw = e isa ArgumentError
        end
        @printf("  [G-F7-1e] a≠0 (bhs=0.5) SCDSourceOps rejected (μ₀ spread=%.3f, non-constant) = %s\n",
                maximum(abs, kerr.gh_pm.mu_0 .- (-1.0)), threw)
        @test threw
        println("  [G-F7-1e] hygiene + firewall + full Perturbation suite = REVIEWER (runtests.jl).")
    end
end
