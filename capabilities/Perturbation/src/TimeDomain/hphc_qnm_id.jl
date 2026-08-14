# hphc_qnm_id.jl — F6-A: constraint-consistent QNM initial data for the HPHC reconstruction hierarchy.
#
# PHYSICS / PROVENANCE. F3_OBSTRUCTION.md established that a BROADBAND (Gaussian) slice cannot supply
# constraint-consistent reconstruction initial data: the on-Σ IVP is under-determined AND dynamically
# non-preservable (the Δ-transport tautology `∂_T X = k_X` gives 0=0, no radial ODE). Its recommended
# unblock is PATH-1: a SINGLE-FREQUENCY / QNM slice. With ψ₄ ∝ e^{−iωt} the whole tower is ∝ e^{−iωt},
# so `∂_T → −iω` is KNOWN and every Δ-transport `∂_T X = k_X` becomes a genuine FIRST-ORDER RADIAL ODE
# (one per angular column). This file builds that tower and seeds it as reconstruction ID.
#
# THE ODE (DERIVED FROM THE CODE, not from prose — F6A_SPEC.md binding order). For a reconstruction
# field X the D4 evolution RHS is (recon_finish!, hphc_reconstruction.jl:324-333)
#     k_X = alg_X  −  pre_thorn_prime_DR·(R·∂_R X + falloff·X) / pre_thorn_prime_DT ,   alg_X = alg_diag·X + src ,
# where `src` consumes ONLY already-solved upstream fields (ð-of-k terms use k_Y = −iωY, known) and
# `alg_diag·X` is X's own algebraic self-coupling. Substituting the SINGLE-FREQUENCY fact ∂_T X = −iωX
# (which is exactly what breaks F3 §2's tautology — ∂_T X is now KNOWN) gives, with d ≡
# pre_thorn_prime_DR/pre_thorn_prime_DT,
#     −iω X = alg_diag·X + src − d·(R·∂_R X + falloff·X)
#  ⟺  d·R·∂_R X = (alg_diag + iω − d·falloff)·X + src            [the linear first-order radial ODE].
# We realize ∂_R by the SAME 4th-order one-sided-closure stencil the evolution uses (NumericsCore.set_d1!),
# so per angular column this is the banded linear system  M·X = −src  with
#     M = diag(alg_diag + iω − d·falloff) − diag(d·R)·D1 ,   D1 = set_d1!(I).
# NOTHING is hardcoded: `alg_diag` is EXTRACTED from each assembler by evaluating it at (self=1, upstream=0);
# `src` by (self=0, upstream given); `d`, `pre_*` come from Spacetime.NPGHP via HPHCReconGHP; `D1` from
# NumericsCore. {T}-generic throughout (Complex{BigFloat} flows).
#
# BC / branch [AXIOM-F6A]. The outgoing-only branch (§V rule: the homogeneous ≡ ingoing piece set to
# zero; Formalism 1403/1501/1513) is realized DISCRETELY by the FD operator's one-sided boundary
# CLOSURES (set_d1!, fd_uniform.jl:46-65) — the identical operator the hyperbolic evolution uses, no
# ghost zone / no extra BC row added. At scri (i=1, R=0) the þ′ divisor `d = (R/cl²)/(2+4MR/cl²) → 0`
# collapses row 1 to `iω·X[1] = −src[1]`, i.e. X[1] = −src[1]/(iω) — which is EXACTLY the RLGP-F §V
# outgoing scri closed form (verified field-by-field in G-F6d). So the discrete "no BC" IS the outgoing
# branch; each field's peeling `falloff` is its _RECON_SBF entry.
#
# a = 0 THIS ROUND (the shipped h5 fixture is Schwarzschild). At a=0: μ₀ = −1/cl² (const), τ₀ = π₀ = 0,
# and the GHP ð/ð′ reduce to pure angular raise/lower (pre_edth_DT = pre_edth = 0), so the coefficient
# matrix M is angular-column-INDEPENDENT (one nx×nx solve serves all ny columns). Kerr QNM ID (a≠0, a
# y-dependent M per column and a nontrivial angular_coef spread) is [FUTURE].
#
# FINDING (honest, decisive — see the test G-F6a/G-F6b). The transport tower satisfies each Δ-transport
# to machine precision (∂_T X = −iωX, ≤1e-15) and its scri values reproduce the §V forms (≤1e-15), BUT
# res_bianchi3(t=0) is NONZERO and resolution-INDEPENDENT (≈0.7, not the truncation floor): the D-direction
# Bianchi-3 CONSTRAINT and the Δ-direction TRANSPORTS are INDEPENDENT discrete identities that do not
# coincide pointwise for the peeling-stripped tower (they differ by a clean O(1): the scri λ ratio is
# exactly −2/3). The PAYOFF is dynamical: under evolution res_bianchi3(t) DECAYS as the QNM envelope
# e^{−|Imω|t} (the coherent single-mode ringdown) instead of GROWING — resolving F3's Candidate-A
# regrowth. res(t=0)=0 alone is NOT sufficient (the :bianchi3 tower below sets it to machine-zero but,
# with λ off its riem-10 transport, it regrows like Candidate-A). See G-F6b's three-way comparison.

# ---------------------------------------------------------------------------------------------------
# (1) ψ₄ QNM seed from the baked fixture  (Id.set_qnm! convention; Id.jl:120-174)
# ---------------------------------------------------------------------------------------------------

# Clenshaw evaluation of a ChebyshevT series (0-based `coef`) at x ∈ [−1,1]; {T}-generic.
@inline function _qnm_cheby(coef::AbstractVector{Complex{T}}, x::T) where {T}
    n = length(coef)
    b1 = zero(Complex{T}); b2 = zero(Complex{T})
    @inbounds for k in n:-1:2
        b1, b2 = coef[k] + 2x * b1 - b2, b1
    end
    return @inbounds coef[1] + x * b1 - b2
end

"""
    hphc_qnm_seed(prob::HPHCTeukolskyProblem, qnm; overtone_n=[0], amp=[complex(1.0)]) -> (; f0, p0, omega)

Build the peeling-stripped ψ₄ QNM seed `(f0, p0) = (ψ, ∂_Tψ)` on `prob`'s grid from the baked fixture
NamedTuple `qnm` (e.g. `QNM_A0_L2_M2_N0`), transcribing `Id.set_qnm!` (Id.jl:120-174) at the D4
`p = −iω f` reduction:

  f(R,y) = ChebyshevT(qnm.radial_coef)(x),  x = 2R/maximum(Rv) − 1  (Id.jl:138,148),
           × Σ_l (−1)^l qnm.angular_coef[l+1]·swal(spin, mv, l+lmin, y)  (Id.jl:149-151),
  then rescaled by amp/max|f| (Id.jl:152-159),  p = −iω f  (Id.jl:167).

**P4 fix (upstream key/shape mismatch, F6_prep_qnm.md §2).** The upstream reader consumes `id_overtone_n`
and `id_qnm_amp` as EQUAL-LENGTH VECTORS (`enumerate(overtone_n)`/`overtone_amp[ni]`, TeukEvolution.jl:297-308),
but ships a scalar example → `KeyError`/`MethodError`. This API takes `overtone_n::Vector` and `amp::Vector`
index-aligned by overtone. For the single-`n` fixture use `overtone_n=[qnm.overtone_n]`, `amp=[A]`; the seed
sums the requested overtones (only `qnm.overtone_n` is present in this fixture, so a length-1 vector).
`Complex{T}` element type flows from `prob` (BigFloat inherited).
"""
function hphc_qnm_seed(prob::HPHCTeukolskyProblem{T}, qnm;
                       overtone_n::AbstractVector{<:Integer}=[qnm.overtone_n],
                       amp::AbstractVector=[one(Complex{T})]) where {T}
    length(overtone_n) == length(amp) ||
        throw(ArgumentError("P4: overtone_n and amp must be equal-length vectors (got $(length(overtone_n)) vs $(length(amp)))"))
    all(n -> n == qnm.overtone_n, overtone_n) ||
        throw(ArgumentError("fixture $(qnm.overtone_n === nothing ? "" : "n=$(qnm.overtone_n)") carries only overtone n=$(qnm.overtone_n); requested $(overtone_n) [FUTURE: multi-overtone h5]"))
    nx, ny, mv = prob.nx, prob.ny, prob.mv
    Rv, Yv = prob.Rv, prob.Yv
    spin = qnm.spin; lmin = qnm.lmin
    rc = Complex{T}.(qnm.radial_coef)
    ac = Complex{T}.(qnm.angular_coef)
    ω = Complex{T}(qnm.omega)
    Rmax = maximum(Rv)
    f0 = zeros(Complex{T}, nx, ny)
    for (ov, a) in zip(overtone_n, amp)      # P4 vector loop (single n here; general shape kept)
        aux = zeros(Complex{T}, nx, ny)
        maxv = zero(T)
        @inbounds for j in 1:ny, i in 1:nx
            rad = _qnm_cheby(rc, 2 * Rv[i] / Rmax - one(T))
            ang = zero(Complex{T})
            for l in 0:(length(ac) - 1)
                ang += T((-1)^l) * ac[l + 1] * swal(spin, mv, l + lmin, Yv[j])
            end
            aux[i, j] = rad * ang
            maxv = max(maxv, abs(aux[i, j]))
        end
        f0 .+= (Complex{T}(a) / maxv) .* aux
    end
    p0 = (-im * ω) .* f0
    return (; f0, p0, omega = ω)
end

"""
    hphc_qnm_amplitude(prob::HPHCTeukolskyProblem, f0, qnm) -> A

The single complex QNM amplitude `A ≡ ₋₂𝒜_lm` calibrated ONCE from the seeded ψ₄ at scri (i=1): since
the peeling-stripped ψ₄ scri value is `A·swal(spin,mv,l,y)` (RLGP-F eq:psi-4-asym stripped, 1407),
`A = f0[1,j] / swal(spin,mv,l,y_j)` (angular-mode-independent; evaluated on the first angular node with
nonzero swal). Used by G-F6d to calibrate the §V oracle.
"""
function hphc_qnm_amplitude(prob::HPHCTeukolskyProblem{T}, f0::AbstractMatrix{Complex{T}}, qnm) where {T}
    spin, mv, l, lmin = qnm.spin, prob.mv, qnm.l, qnm.lmin
    for j in 1:prob.ny
        s = swal(spin, mv, l, prob.Yv[j])
        abs(s) > sqrt(eps(T)) && return f0[1, j] / s
    end
    error("no nonzero swal node for amplitude calibration")
end

# ---------------------------------------------------------------------------------------------------
# (2) the per-field radial linear-solve  M·X = −src   (the ∂_T→−iω transport ODE, code-discretized)
# ---------------------------------------------------------------------------------------------------

# D1 = the NumericsCore.set_d1! differentiation matrix (apply the stencil to the identity columns): the
# EXACT operator the evolution uses (interior 4th-order central + one-sided 4th-order boundary closures).
function _qnm_d1_matrix(::Type{T}, nx::Integer, dr::T) where {T}
    D1 = zeros(Complex{T}, nx, nx)
    set_d1!(D1, Matrix{Complex{T}}(Id, nx, nx), dr)
    return D1
end

# Solve M·X = −src for a single reconstruction field. `algd` = the per-node algebraic self-coefficient
# (alg_diag, EXTRACTED from the assembler), a length-nx vector (a=0 ⇒ angular-column-independent, asserted).
# `dco` = pre_thorn_prime_DR/pre_thorn_prime_DT per node. Builds the banded M and solves all ny columns.
function _qnm_radial_solve(algd::AbstractVector{Complex{T}}, src::AbstractMatrix{Complex{T}},
                           falloff::Integer, Rv::AbstractVector{T}, dco::AbstractVector{T},
                           D1::AbstractMatrix{Complex{T}}, ω::Complex{T}) where {T}
    nx = length(Rv)
    M = zeros(Complex{T}, nx, nx)
    @inbounds for i in 1:nx
        M[i, i] += algd[i] + im * ω - dco[i] * falloff
        for k in 1:nx
            M[i, k] += -dco[i] * Rv[i] * D1[i, k]
        end
    end
    return M \ (-src)
end

# Guarded extraction of the angular-column-independent (a=0) coefficient vectors from a per-mode GHP set.
# Asserts y-independence (loud error on a≠0 — [FUTURE: per-column M for Kerr]).
function _qnm_coeff_vectors(gh::HPHCReconGHP{T}) where {T}
    dco = real.(gh.pre_thorn_prime_DR ./ gh.pre_thorn_prime_DT)   # nx×ny
    yspread = maximum(abs, dco .- dco[:, 1:1])
    yspread <= sqrt(eps(T)) ||
        throw(ArgumentError("hphc_qnm_tower requires angular-column-independent coefficients (a=0); " *
                            "y-spread=$(yspread) ⇒ a≠0 is [FUTURE]. Use a Schwarzschild (bhs=0) problem."))
    return real.(dco[:, 1])   # nx real vector
end

# ---------------------------------------------------------------------------------------------------
# (3) the tower assembler — for ONE mode, in D4 DAG order (ψ₄ given ⇒ ψ₃,λ ⇒ ψ₂,h_m̄m̄ ⇒ π ⇒ h_lm̄)
# ---------------------------------------------------------------------------------------------------

"""
    hphc_qnm_tower_mode(gh, prob1, f4, ω; lam_from=:transport) -> (; psi3, lam, psi2, hmbmb, pi, hlmb)

Construct the single-mode constraint-consistent QNM reconstruction tower from the ψ₄ seed `f4` at
frequency `ω`, by the ∂_T→−iω radial linear solve for each field in DAG order (ψ₄ → {ψ₃, λ} →
{ψ₂, h_m̄m̄} → π → h_lm̄). Each field's `alg_diag` and `src` are EXTRACTED from the LIVE D4 assembler
(`recon_set_*_k!`) at (self=1, upstream=0) and (self=0, upstream) respectively; the ð-sources use the
LIVE `recon_edth!` with DT = −iω·field (a≠0-ready; DT is annihilated at a=0). `gh` is the mode's
`HPHCReconGHP`, `prob1` its `HPHCTeukolskyProblem`.

`lam_from` selects the λ branch (the F6 physics fork, honest):
  * `:transport` (default) — λ solves its riem-10 Δ-transport (dynamically preserved; res_bianchi3(t)
    then DECAYS as the QNM envelope, G-F6b).
  * `:bianchi3` — λ solved ALGEBRAICALLY from the D-direction constraint res_bianchi3 = 0 (makes
    res_bianchi3(t=0)=0 to machine precision, but λ is then OFF its transport, so under evolution it
    regrows like F3 Candidate-A — the negative control in G-F6b).
"""
function hphc_qnm_tower_mode(gh::HPHCReconGHP{T}, prob1::HPHCTeukolskyProblem{T},
                             f4::AbstractMatrix{Complex{T}}, ω::Complex{T};
                             lam_from::Symbol=:transport) where {T}
    nx, ny = prob1.nx, prob1.ny
    Rv, dr = prob1.Rv, prob1.dr
    mu0, tau0, pi0 = gh.mu_0, gh.tau_0, gh.pi_0
    dco = _qnm_coeff_vectors(gh)
    D1 = _qnm_d1_matrix(T, nx, dr)
    Z = zeros(Complex{T}, nx, ny); O = ones(Complex{T}, nx, ny)
    tmp = zeros(Complex{T}, nx, ny); raised = zeros(Complex{T}, nx, ny)
    sbf = _RECON_SBF
    # ð of an already-solved upstream field g of GHP type (spin,boost): DT = −iω g (a≠0-ready).
    edth(g, spin, boost) = (o = zeros(Complex{T}, nx, ny);
                            recon_edth!(o, raised, g, (-im * ω) .* g, gh, spin, boost); o)

    # --- ψ₃ : alg_diag from recon_set_psi3_k!(self=1,ψ4=0,ðψ4=0); src at (self=0, ψ4, ðψ4) ---
    p4e = edth(f4, sbf.psi4_f[1], sbf.psi4_f[2])
    a3 = zeros(Complex{T}, nx, ny); recon_set_psi3_k!(a3, O, Z, Z, mu0, tau0, Rv)
    s3 = zeros(Complex{T}, nx, ny); recon_set_psi3_k!(s3, Z, f4, p4e, mu0, tau0, Rv)
    psi3 = _qnm_radial_solve(a3[:, 1], s3, sbf.psi3[3], Rv, dco, D1, ω)

    # --- λ ---
    al = zeros(Complex{T}, nx, ny); recon_set_lam_k!(al, O, Z, mu0, Rv)
    sl = zeros(Complex{T}, nx, ny); recon_set_lam_k!(sl, Z, f4, mu0, Rv)
    lam_tr = _qnm_radial_solve(al[:, 1], sl, sbf.lam[3], Rv, dco, D1, ω)

    lam = lam_tr
    if lam_from === :bianchi3
        # λ from the D-constraint res_bianchi3 = ð′ψ₃ + 4π₀ψ₃ − þψ₄ + ρ₀ψ₄ − 3ψ₂₀λ = 0.
        psi3k = zeros(Complex{T}, nx, ny)
        recon_set_psi3_k!(psi3k, psi3, f4, p4e, mu0, tau0, Rv); recon_finish!(psi3k, psi3, tmp, gh, sbf.psi3[3], dr)
        p3ep = zeros(Complex{T}, nx, ny); low = zeros(Complex{T}, nx, ny)
        recon_edth_prime!(p3ep, low, psi3, psi3k, gh, sbf.psi3[1], sbf.psi3[2])
        p4dr = zeros(Complex{T}, nx, ny); set_d1!(p4dr, f4, dr)
        p4th = zeros(Complex{T}, nx, ny)
        recon_thorn!(p4th, f4, (-im * ω) .* f4, p4dr, gh, sbf.psi4_f[1], sbf.psi4_f[2], sbf.psi4_f[3])
        lam = (p3ep .+ 4 .* pi0 .* psi3 .- p4th .+ gh.rho_0 .* f4) ./ (3 .* gh.psi2_0)
    elseif lam_from !== :transport
        throw(ArgumentError("lam_from must be :transport or :bianchi3 (got $(lam_from))"))
    end

    # --- ψ₂ (from ψ₃) : needs ðψ₃ ---
    p3e = edth(psi3, sbf.psi3[1], sbf.psi3[2])
    a2 = zeros(Complex{T}, nx, ny); recon_set_psi2_k!(a2, O, Z, Z, mu0, tau0, Rv)
    s2 = zeros(Complex{T}, nx, ny); recon_set_psi2_k!(s2, Z, psi3, p3e, mu0, tau0, Rv)
    psi2 = _qnm_radial_solve(a2[:, 1], s2, sbf.psi2[3], Rv, dco, D1, ω)

    # --- h_m̄m̄ (from λ) ---
    ah = zeros(Complex{T}, nx, ny); recon_set_hmbmb_k!(ah, O, Z, mu0, Rv)
    sh = zeros(Complex{T}, nx, ny); recon_set_hmbmb_k!(sh, Z, lam, mu0, Rv)
    hmbmb = _qnm_radial_solve(ah[:, 1], sh, sbf.hmbmb[3], Rv, dco, D1, ω)

    # --- π (from λ, h_m̄m̄, ψ₃) : no self-coupling ⇒ alg_diag = 0 ---
    api = zeros(Complex{T}, nx, ny)   # recon_set_pi_k! has no π argument ⇒ alg_diag ≡ 0
    spi = zeros(Complex{T}, nx, ny); recon_set_pi_k!(spi, lam, hmbmb, psi3, mu0, tau0, pi0, Rv)
    pil = _qnm_radial_solve(api[:, 1], spi, sbf.pi[3], Rv, dco, D1, ω)

    # --- h_lm̄ (from h_m̄m̄, π) ---
    ahl = zeros(Complex{T}, nx, ny); recon_set_hlmb_k!(ahl, O, Z, Z, mu0, tau0, Rv)
    shl = zeros(Complex{T}, nx, ny); recon_set_hlmb_k!(shl, Z, hmbmb, pil, mu0, tau0, Rv)
    hlmb = _qnm_radial_solve(ahl[:, 1], shl, sbf.hlmb[3], Rv, dco, D1, ω)

    return (; psi3, lam, psi2, hmbmb, pi = pil, hlmb)
end

"""
    hphc_qnm_initialdata(prob::HPHCReconProblem, qnm; amp=complex(0.1), lam_from=:transport)
        -> (; pm, nm, omega, amp, A_pm, A_nm)

Assemble the full ±m constraint-consistent QNM reconstruction ID: ψ₄(f,p) seed + the transport tower
(ψ₃,λ,ψ₂,h_m̄m̄,π,h_lm̄) for BOTH modes, then the ±m-coupled μh_ll (its cross-mode source consumes the
partner mode's ð h_lm̄ etc. — recon_set_muhll_k!). Returns per-mode NamedTuples over `_RECON_FIELDS`
(the 9-field slot order) ready to pack into the 18-field state. `amp` is the single QNM amplitude used
for both modes (Id.set_qnm! normalizes |ψ₄|max→amp); `A_pm/A_nm` are the calibrated ₋₂𝒜.
"""
function hphc_qnm_initialdata(prob::HPHCReconProblem{T}, qnm;
                              amp=complex(0.1, 0.0), lam_from::Symbol=:transport) where {T}
    ampv = [Complex{T}(amp)]
    ov = [qnm.overtone_n]
    seed_p = hphc_qnm_seed(prob.pm, qnm; overtone_n=ov, amp=ampv)
    seed_m = hphc_qnm_seed(prob.nm, qnm; overtone_n=ov, amp=ampv)
    ω = seed_p.omega
    f4p, p4p = seed_p.f0, seed_p.p0
    f4m, p4m = seed_m.f0, seed_m.p0
    tp = hphc_qnm_tower_mode(prob.gh_pm, prob.pm, f4p, ω; lam_from=lam_from)
    tm = hphc_qnm_tower_mode(prob.gh_nm, prob.nm, f4m, ω; lam_from=lam_from)
    # μh_ll : ±m coupled. Build both modes' ð(h_lm̄,π,h_m̄m̄), then recon_set_muhll_k! + finish.
    # CRITICAL: the PARTNER mode's ð uses the PARTNER's GHP (its raise matrix differs by m), exactly as
    # the evolution does (_recon_muhll! reads sc_p from the partner's own _recon_chain_upto_hlmb!).
    nx = prob.pm.nx
    dco = _qnm_coeff_vectors(prob.gh_pm); D1 = _qnm_d1_matrix(T, nx, prob.pm.dr)
    dco_m = _qnm_coeff_vectors(prob.gh_nm)
    muhll_p = _qnm_muhll(prob.gh_pm, prob.gh_nm, prob.pm, tp, tm, ω, prob.pm.Rv, dco, D1)
    muhll_m = _qnm_muhll(prob.gh_nm, prob.gh_pm, prob.nm, tm, tp, ω, prob.nm.Rv, dco_m, D1)
    A_pm = hphc_qnm_amplitude(prob.pm, f4p, qnm)
    A_nm = hphc_qnm_amplitude(prob.nm, f4m, qnm)
    pm = (; psi4_f=f4p, psi4_p=p4p, psi3=tp.psi3, lam=tp.lam, psi2=tp.psi2,
          hmbmb=tp.hmbmb, pi=tp.pi, hlmb=tp.hlmb, muhll=muhll_p)
    nm = (; psi4_f=f4m, psi4_p=p4m, psi3=tm.psi3, lam=tm.lam, psi2=tm.psi2,
          hmbmb=tm.hmbmb, pi=tm.pi, hlmb=tm.hlmb, muhll=muhll_m)
    return (; pm, nm, omega=ω, amp=Complex{T}(amp), A_pm, A_nm)
end

# μh_ll radial solve for one mode using its own tower `me` (with own GHP `gh`) and the partner's tower
# `pa` (with the partner's GHP `gh_p` — its raise matrix differs by m) — the ±m coupling.
function _qnm_muhll(gh::HPHCReconGHP{T}, gh_p::HPHCReconGHP{T}, prob1::HPHCTeukolskyProblem{T}, me, pa,
                    ω::Complex{T}, Rv::AbstractVector{T}, dco::AbstractVector{T},
                    D1::AbstractMatrix{Complex{T}}) where {T}
    nx, ny = prob1.nx, prob1.ny
    mu0, tau0, pi0 = gh.mu_0, gh.tau_0, gh.pi_0
    Z = zeros(Complex{T}, nx, ny); O = ones(Complex{T}, nx, ny)
    raised = zeros(Complex{T}, nx, ny)
    sbf = _RECON_SBF
    edth(g, spin, boost, gg) = (o = zeros(Complex{T}, nx, ny);
                                recon_edth!(o, raised, g, (-im * ω) .* g, gg, spin, boost); o)
    hlmb_e  = edth(me.hlmb,  sbf.hlmb[1], sbf.hlmb[2], gh)     # own ð h_lm̄  (own GHP)
    pi_e    = edth(me.pi,    sbf.pi[1],   sbf.pi[2],   gh)     # own ð π
    hmbmb_e_nm = edth(pa.hmbmb, sbf.hmbmb[1], sbf.hmbmb[2], gh_p)   # partner ð h_m̄m̄  (partner GHP)
    hlmb_e_nm  = edth(pa.hlmb,  sbf.hlmb[1],  sbf.hlmb[2],  gh_p)   # partner ð h_lm̄
    # alg_diag (self=1, all sources 0)
    am = zeros(Complex{T}, nx, ny)
    recon_set_muhll_k!(am, O, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, mu0, tau0, pi0, Rv)
    # src (self=0, sources given)
    sm = zeros(Complex{T}, nx, ny)
    recon_set_muhll_k!(sm, Z, hlmb_e, me.hlmb, pi_e, me.pi, me.psi2,
                       pa.pi, hmbmb_e_nm, pa.hmbmb, hlmb_e_nm, pa.hlmb, mu0, tau0, pi0, Rv)
    return _qnm_radial_solve(am[:, 1], sm, sbf.muhll[3], Rv, dco, D1, ω)
end

# ---------------------------------------------------------------------------------------------------
# (4) evolution driver — pack the QNM ID into the 18-field state and step (ONE rk4! + P3 filter)
# ---------------------------------------------------------------------------------------------------

"""
    solve_hphc_recon_qnm(prob::HPHCReconProblem, qnm; nt, ts, amp=complex(0.1),
                         lam_from=:transport, seed_tower=true) -> (; times, pm, nm, res3)

Evolve the HPHC reconstruction problem seeded with the constraint-consistent QNM initial data (via
[`hphc_qnm_initialdata`](@ref)) for `nt` RK4 steps, recording every `ts` steps. Identical stepping to
[`solve_hphc_recon`](@ref) (one `rk4!` over the packed 18-field ±m state; P3 KO+killtop post-filter on
{ψ₄f,ψ₄p,λ,ψ₃}), but (i) the ψ₄ seed is the QNM Chebyshev profile (not a Gaussian), (ii) the transport
tower is injected (unless `seed_tower=false`, the zero-tower control), and (iii) it additionally returns
`res3` — the +m `res_bianchi3` interior-max time series (the G-F6a/G-F6b diagnostic).

`lam_from=:transport` (dynamically-preserved λ ⇒ res3 DECAYS as e^{−|Imω|t}) or `:bianchi3`
(res3(t=0)=0 but regrows). `seed_tower=false` seeds only ψ₄ (tower=0 ⇒ res3 regrows — the F3 baseline).
"""
function solve_hphc_recon_qnm(prob::HPHCReconProblem{T}, qnm; nt::Integer, ts::Integer,
                              amp=complex(0.1, 0.0), lam_from::Symbol=:transport,
                              seed_tower::Bool=true) where {T}
    nx, ny = prob.pm.nx, prob.pm.ny
    dt = prob.pm.dt
    u = zeros(Complex{T}, nx * ny * 18)
    U = reshape(u, nx, ny, 18)
    id = hphc_qnm_initialdata(prob, qnm; amp=amp, lam_from=lam_from)
    @views U[:, :, 1] .= id.pm.psi4_f;  @views U[:, :, 2] .= id.pm.psi4_p
    @views U[:, :, 10] .= id.nm.psi4_f; @views U[:, :, 11] .= id.nm.psi4_p
    if seed_tower
        for (k, f) in enumerate(_RECON_FIELDS)
            k <= 2 && continue
            @views U[:, :, k]     .= getfield(id.pm, f)
            @views U[:, :, 9 + k] .= getfield(id.nm, f)
        end
    end
    sc_pm = _ReconScratch(T, nx, ny); sc_nm = _ReconScratch(T, nx, ny)
    epsKO = T(1) / 2; tmpf = zeros(Complex{T}, nx, ny)

    function rhs!(du, uu, ::Real)
        UU = reshape(uu, nx, ny, 18); dUU = reshape(du, nx, ny, 18)
        fl_pm = _mode_views(UU, 0); kk_pm = _mode_views(dUU, 0)
        fl_nm = _mode_views(UU, 9); kk_nm = _mode_views(dUU, 9)
        _recon_chain_upto_hlmb!(kk_pm, fl_pm, prob.pm, prob.gh_pm, sc_pm)
        _recon_chain_upto_hlmb!(kk_nm, fl_nm, prob.nm, prob.gh_nm, sc_nm)
        _recon_muhll!(kk_pm, fl_pm, sc_pm, fl_nm, sc_nm, prob.pm, prob.gh_pm)
        _recon_muhll!(kk_nm, fl_nm, sc_nm, fl_pm, sc_pm, prob.nm, prob.gh_nm)
        return du
    end
    function postfilter!(uu)
        UU = reshape(uu, nx, ny, 18)
        for (off, kt) in ((0, prob.kt_pm), (9, prob.kt_nm))
            for (slot, s) in ((1, -2), (2, -2), (4, -2), (3, -1))
                v = @view UU[:, :, off + slot]
                ko_filter!(v, tmpf, epsKO); tmpf .= v; angular_matrix_mult!(v, tmpf, kt[s])
            end
        end
        return nothing
    end

    snaps = Dict(f => Matrix{Complex{T}}[] for f in _RECON_FIELDS)
    snaps_nm = Dict(f => Matrix{Complex{T}}[] for f in _RECON_FIELDS)
    times = T[]; res3 = T[]
    function cb(uu, t, step)
        step >= 1 && postfilter!(uu)
        if step % ts == 0
            UU = reshape(uu, nx, ny, 18)
            for (k, f) in enumerate(_RECON_FIELDS)
                push!(snaps[f], copy(@view UU[:, :, k]))
                push!(snaps_nm[f], copy(@view UU[:, :, 9 + k]))
            end
            snap = NamedTuple{_RECON_FIELDS}(Tuple(snaps[f][end] for f in _RECON_FIELDS))
            res = hphc_recon_residual3(prob, snap)
            push!(times, T(t)); push!(res3, maximum(abs, @view res[4:nx-3, :]))
        end
        return nothing
    end
    rk4!(rhs!, u, (zero(T), nt * dt), dt; callback=cb)
    pm = NamedTuple{_RECON_FIELDS}(Tuple(snaps[f] for f in _RECON_FIELDS))
    nm = NamedTuple{_RECON_FIELDS}(Tuple(snaps_nm[f] for f in _RECON_FIELDS))
    return (; times, pm, nm, res3)
end
