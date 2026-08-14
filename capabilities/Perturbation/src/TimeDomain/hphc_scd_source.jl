# hphc_scd_source.jl — the second-order Teukolsky source (𝔰_d, 𝔰_t, 𝒮) bilinear assembler.
#
# 𝒮 = d₄⁽⁰⁾𝔰_d + d₃⁽⁰⁾𝔰_t, d₄⁽⁰⁾ = Δ+4μ+μ̄, d₃⁽⁰⁾ = ð′+4π−τ̄, in outgoing radiation gauge with the
# RLGP background tetrad (arXiv:2010.00162 eq:schematic_source + eq:def_mathfrak_s, tex:924-981;
# the full 𝔰_d/𝔰_t operator expansions are double-entried in test_td_hphc_scd.jl). Unsuperscripted
# μ,μ̄,π,π̄,τ,τ̄ are BACKGROUND NP scalars (Spacetime.NPGHP); every ⁽¹⁾/h_· field is source input.
#
# BILINEAR STRUCTURE (owner tex:1271-1341): 𝒮 is built as an ordered bilinear over two sector
# bundles — the Weyl factor {Ψ₂,Ψ₃,Ψ₄} always from B1 (the unreflected H₊ slot), the metric/
# operator factor from B2. The H₋ fields (h_lm, h_mm, π̄) live ONLY in the metric slot, so
# Π₋₋𝒮 = 0 is manifest and the mixed channel is single-counted (owner tex:1332-1340). The
# diagonal source is 𝒮[B,B]; the general symmetric bilinear is ½(𝒮[B1,B2]+𝒮[B2,B1]).
#
# MODE MIXING: bundles at azimuthal a1, a2 land the product in the (a1+a2) sector. The H₋ metric
# fields are the ±m conjugate partners (h_lm = conj h_lm̄[−m], …; hphc_reconstruction.jl:287-312).
#
# GHP application uses the same coefficient machinery as the reconstruction (Spacetime.NPGHP);
# at a=0 the pre_* fields are m-independent and ð/ð′ reduce to pure raise/lower; for a QNM bundle
# the DT argument is −iωX exactly. New raise/lower matrices at the H₋ positive spins are gated
# equal to the locked reconstruction matrices where they overlap.
#
# a=0 this round (μ₀ = −1/cl² constant, τ₀=π₀=0 ⇒ h_ll = −cl²·μh_ll cleanly); Kerr a≠0 is
# [FUTURE]. {T}-generic. This file assembles 𝒮 only — no Ψ₄⁽²⁾ solve here.

# ---------------------------------------------------------------------------------------------------
# (0) the GHP-type / peeling table of the source's composite quantities (DERIVED, not guessed)
# ---------------------------------------------------------------------------------------------------
# GHP (spin s, boost b) and peeling `falloff` of every quantity a Δ/ð/ð′ is applied to. Derived from
# the leg content (l:b+1, n:b−1, m:s+1, m̄:s−1; Δ=þ′ adds (0,−1), ð adds (+1,0), ð′ adds (−1,0)) and
# the _RECON_SBF peeling of the tower (h_ll = μh_ll/μ₀, μ₀ boost −1 ⇒ h_ll boost +2 falloff 3; the
# conjugate-partner H_− fields inherit the falloff of their +m partner). Cross-checked: 𝔰_d type
# (−2,−1) fall 1, 𝔰_t (−1,−2) fall 1, and 𝒮 = d₄𝔰_d = d₃𝔰_t ⇒ type (−2,−2) fall 1 == the scd_f field
# (TeukEvolution.jl:208-216, spin/boost/falloff = −2,−2,1). Values used below:
#   Ψ₄ (−2,−2,fall1)  Ψ₃ (−1,−1,fall2)  Ψ₂ (0,0)   λ (−2,−1)   π (−1,0)   π̄ (+1,0)
#   h_ll (0,2,fall3)  h_lm̄ (−1,+1,fall2)  h_m̄m̄ (−2,0)   h_lm (+1,+1,fall2)   h_mm (+2,0)
#   𝔰_d (−2,−1,fall1)   𝔰_t (−1,−2,fall1)

# ---------------------------------------------------------------------------------------------------
# (1) source operator machinery — coefficients SOURCED from HPHCReconGHP, ladders built per (spin,m)
# ---------------------------------------------------------------------------------------------------

"""
    SCDSourceOps{T}

Operator machinery for the second-order source assembler: the m-INDEPENDENT `HPHCReconGHP` coefficient
fields (`pre_*`, `mu_0`, `tau_0`, `pi_0`), the radial grid `(Rv,dr,cl)`, and a lazily-built cache of
spin raise/lower collocation matrices keyed by `(dir,spin,m_ang)` (dir ∈ {:raise,:lower}). The GHP
applicators [`scd_edth!`](@ref)/[`scd_edth_prime!`](@ref)/[`scd_thornprime!`](@ref) reproduce the
locked `recon_*` forms but accept an arbitrary-spin ladder (needed for the H_− positive spins).
Build with [`SCDSourceOps(::HPHCReconProblem)`](@ref).
"""
struct SCDSourceOps{T<:AbstractFloat}
    gh::HPHCReconGHP{T}
    Rv::Vector{T}
    dr::T
    cl::T
    ny::Int
    ladder::Dict{Tuple{Symbol,Int,Int},Matrix{T}}
end

"""
    SCDSourceOps(prob::HPHCReconProblem) -> SCDSourceOps

Assemble the source operator machinery on the reconstruction problem. Uses `prob.gh_pm` as the
(m-independent) coefficient source and asserts the a=0 axiom `mu_0 = −1/cl²` CONSTANT (
loud error at a≠0 — the μ₀→0 near-scri unpacking of `h_ll` is [FUTURE]).
"""
function SCDSourceOps(prob::HPHCReconProblem{T}) where {T}
    gh = prob.gh_pm
    cl = prob.pm.cl
    tgt = -one(T) / cl^2
    spread = maximum(abs, gh.mu_0 .- tgt)
    spread <= sqrt(eps(T)) ||
        throw(ArgumentError("SCDSourceOps: [AXIOM-F6A-a0] requires mu_0 = −1/cl² CONSTANT (a=0); " *
                            "got spread=$(spread) ⇒ a≠0 h_ll=μh_ll/μ₀ unpacking is [FUTURE]. Use bhs=0."))
    return SCDSourceOps{T}(gh, prob.pm.Rv, prob.pm.dr, T(cl), prob.pm.ny,
                           Dict{Tuple{Symbol,Int,Int},Matrix{T}}())
end

"Spin raise (`:raise`, ð) or lower (`:lower`, ð′) collocation matrix at `(spin, m_ang)`, cached."
function scd_ladder(ops::SCDSourceOps{T}, dir::Symbol, spin::Integer, m_ang::Integer) where {T}
    get!(ops.ladder, (dir, Int(spin), Int(m_ang))) do
        dir === :raise ? swal_raising_matrix(T, ops.ny, spin, m_ang) :
        dir === :lower ? swal_lowering_matrix(T, ops.ny, spin, m_ang) :
        throw(ArgumentError("dir must be :raise or :lower"))
    end
end

"""
    scd_edth!(out, raised, f, DT, ops, spin, boost, m_ang) -> out

GHP `ð` on a type-`(spin,boost)` field `f` at azimuthal `m_ang`, with `DT = ∂_T f`. Same form as the
locked `recon_edth!` (`out = pre_edth_DT·DT + pre_edth_raised·(ð_swal f) + p·pre_edth·f`,
`p = ghp_weights(spin,boost)[1]`) but with the ladder `scd_ladder(ops,:raise,spin,m_ang)`, so any spin
is admissible. Gated === `recon_edth!` where the ladders coincide.
"""
function scd_edth!(out, raised, f, DT, ops::SCDSourceOps, spin::Integer, boost::Integer, m_ang::Integer)
    angular_matrix_mult!(raised, f, scd_ladder(ops, :raise, spin, m_ang))
    p, _ = NPGHP.ghp_weights(spin, boost)
    gh = ops.gh
    @inbounds for j in axes(out, 2), i in axes(out, 1)
        out[i, j] = gh.pre_edth_DT[i, j] * DT[i, j] + gh.pre_edth_raised[i, j] * raised[i, j] +
                    p * gh.pre_edth[i, j] * f[i, j]
    end
    return out
end

"""
    scd_edth_prime!(out, lowered, f, DT, ops, spin, boost, m_ang) -> out

GHP `ð′` on a type-`(spin,boost)` field at azimuthal `m_ang` (`q = ghp_weights(spin,boost)[2]`). Same
form as the locked `recon_edth_prime!` with the ladder `scd_ladder(ops,:lower,spin,m_ang)`.
"""
function scd_edth_prime!(out, lowered, f, DT, ops::SCDSourceOps, spin::Integer, boost::Integer, m_ang::Integer)
    angular_matrix_mult!(lowered, f, scd_ladder(ops, :lower, spin, m_ang))
    _, q = NPGHP.ghp_weights(spin, boost)
    gh = ops.gh
    @inbounds for j in axes(out, 2), i in axes(out, 1)
        out[i, j] = gh.pre_edth_prime_DT[i, j] * DT[i, j] + gh.pre_edth_prime_lowered[i, j] * lowered[i, j] +
                    q * gh.pre_edth_prime[i, j] * f[i, j]
    end
    return out
end

"""
    scd_thornprime!(out, f, DT, DR, ops, falloff) -> out

GHP `Δ` (=þ′) on a peeling-`falloff` field `f`, `DT = ∂_T f`. Verbatim `recon_thorn_prime!` form
(`out = pre_thorn_prime_DT·DT + pre_thorn_prime_DR·(R·∂_R f + falloff·f)`); m-independent. Computes
`∂_R f` into `DR` via `NumericsCore.set_d1!` (the same one-sided-closure stencil the evolution uses).
"""
function scd_thornprime!(out, f, DT, DR, ops::SCDSourceOps, falloff::Integer)
    set_d1!(DR, f, ops.dr)
    gh = ops.gh
    R = ops.Rv
    @inbounds for j in axes(out, 2), i in axes(out, 1)
        out[i, j] = gh.pre_thorn_prime_DT[i, j] * DT[i, j] +
                    gh.pre_thorn_prime_DR[i, j] * (R[i] * DR[i, j] + falloff * f[i, j])
    end
    return out
end

# ---------------------------------------------------------------------------------------------------
# (2) the sector bundle — all source ingredients at one azimuthal number m_ang and QNM frequency ω
# ---------------------------------------------------------------------------------------------------

"""
    SCDBundle{T}

One quadratic-source factor: the eleven first-order ingredient fields at a definite azimuthal number
`m_ang` and QNM frequency `omega` (so `∂_T X = −iω X` uniformly). Split by the H_±
grading:
  * H_+ Weyl / source slot   : `psi4, psi3, psi2`
  * H_+ metric / operator slot: `hll, hlmb, hmbmb, lam, pi`
  * H_− metric slot (conjugate-partner): `hlm, hmm, pibar`
Built with [`scd_bundle`](@ref) (from an `hphc_qnm_initialdata` result) or by hand (the formal
H_± monomial towers of the structural gate).
"""
struct SCDBundle{T<:AbstractFloat}
    psi4::Matrix{Complex{T}}
    psi3::Matrix{Complex{T}}
    psi2::Matrix{Complex{T}}
    hll::Matrix{Complex{T}}
    hlmb::Matrix{Complex{T}}
    hmbmb::Matrix{Complex{T}}
    lam::Matrix{Complex{T}}
    pi::Matrix{Complex{T}}
    hlm::Matrix{Complex{T}}
    hmm::Matrix{Complex{T}}
    pibar::Matrix{Complex{T}}
    omega::Complex{T}
    m_ang::Int
end

const _SCD_BUNDLE_FIELDS = (:psi4, :psi3, :psi2, :hll, :hlmb, :hmbmb, :lam, :pi, :hlm, :hmm, :pibar)

"Zero bundle of grid size `(nx,ny)` at `(omega, m_ang)` — the base for hand-built formal towers."
function scd_zero_bundle(::Type{T}, nx::Integer, ny::Integer, omega, m_ang::Integer) where {T<:AbstractFloat}
    Z() = zeros(Complex{T}, nx, ny)
    return SCDBundle{T}(Z(), Z(), Z(), Z(), Z(), Z(), Z(), Z(), Z(), Z(), Z(), Complex{T}(omega), Int(m_ang))
end

"""
    scd_bundle(id, block, m_pos, ops) -> SCDBundle

Build a sector bundle from an `hphc_qnm_initialdata` result `id` (NamedTuples `id.pm`,`id.nm` over
`_RECON_FIELDS` + `id.omega`) with positive azimuthal number `m_pos` (the `+m` of the tower).
`block=:pm` gives the `+m` sector (Weyl/H_+ from `id.pm`, H_− = `conj(id.nm)`); `block=:nm` gives the
`−m` sector (Weyl/H_+ from `id.nm`, H_− = `conj(id.pm)`) — the ±m conjugate-partner mechanism
(hphc_reconstruction.jl:287-312). `h_ll = μh_ll/mu_0` (a=0: ×(−cl²)).
"""
function scd_bundle(id, block::Symbol, m_pos::Integer, ops::SCDSourceOps{T}) where {T}
    self, part = block === :pm ? (id.pm, id.nm) :
                 block === :nm ? (id.nm, id.pm) :
                 throw(ArgumentError("block must be :pm or :nm"))
    m_ang = block === :pm ? Int(m_pos) : -Int(m_pos)
    ω = Complex{T}(id.omega)
    mu = ops.gh.mu_0
    hll = self.muhll ./ mu                               # h_ll = μh_ll / μ₀   (a=0: μ₀=−1/cl² const)
    return SCDBundle{T}(copy(self.psi4_f), copy(self.psi3), copy(self.psi2),
                        hll, copy(self.hlmb), copy(self.hmbmb), copy(self.lam), copy(self.pi),
                        conj.(part.hlmb), conj.(part.hmbmb), conj.(part.pi), ω, m_ang)
end

# ---------------------------------------------------------------------------------------------------
# (3) 𝔰_d, 𝔰_t assemblers — the eq:def_mathfrak_s bilinears (Weyl←B1, metric←B2)
# ---------------------------------------------------------------------------------------------------

# per-assembly scratch (radial/angular deriv outputs + operator outputs), allocated once.
struct _SCDScratch{T}
    DR::Matrix{Complex{T}}; raised::Matrix{Complex{T}}; lowered::Matrix{Complex{T}}
    op::Matrix{Complex{T}}
end
_SCDScratch(::Type{T}, nx, ny) where {T} = _SCDScratch{T}(ntuple(_ -> zeros(Complex{T}, nx, ny), 4)...)

"""
    scd_sd!(sd, B1, B2, ops, sc) -> sd

Assemble 𝔰_d (arXiv:2010.00162 eq:def_mathfrak_s) into `sd`: Weyl factor {Ψ₄,Ψ₃,Ψ₂} from `B1`, metric
factor {h_ll,h_lm̄,h_m̄m̄,h_lm,λ,π} from `B2`. `DT = −iω X` per factor; Δ/ð/ð′ via `ops`. `sd` has GHP
type (−2,−1), falloff 1.
"""
function scd_sd!(sd, B1::SCDBundle{T}, B2::SCDBundle{T}, ops::SCDSourceOps{T}, sc::_SCDScratch{T}) where {T}
    a1, a2 = B1.m_ang, B2.m_ang
    ω1, ω2 = B1.omega, B2.omega
    mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; tau0 = ops.gh.tau_0
    cpi = conj.(pi0); ctau = conj.(tau0)
    op = sc.op
    fill!(sd, zero(Complex{T}))

    # T_a : ½ h_ll (Δ+μ)Ψ₄
    scd_thornprime!(op, B1.psi4, (-im * ω1) .* B1.psi4, sc.DR, ops, 1)                 # ΔΨ₄ (fall 1)
    @inbounds @. sd += (1 // 2) * B2.hll * (op + mu * B1.psi4)
    # T_b : Ψ₄ · ½ (ð+π̄+2τ) h_lm̄
    scd_edth!(op, sc.raised, B2.hlmb, (-im * ω2) .* B2.hlmb, ops, -1, 1, a2)           # ð h_lm̄
    @inbounds @. sd += B1.psi4 * ((1 // 2) * (op + (cpi + 2 * tau0) * B2.hlmb))
    # T_c : Ψ₄ · (Δ−μ+μ̄) h_ll
    scd_thornprime!(op, B2.hll, (-im * ω2) .* B2.hll, sc.DR, ops, 3)                   # Δ h_ll (fall 3)
    @inbounds @. sd += B1.psi4 * (op + (-mu + cmu) * B2.hll)
    # T_d : Ψ₄ · (−½)(ð′−5π−4τ̄) h_lm       [H_−]
    scd_edth_prime!(op, sc.lowered, B2.hlm, (-im * ω2) .* B2.hlm, ops, 1, 1, a2)       # ð′ h_lm
    @inbounds @. sd += B1.psi4 * (-(1 // 2) * (op + (-5 * pi0 - 4 * ctau) * B2.hlm))
    # T_e : −½ Ψ₃ · (ð+π̄+τ) h_m̄m̄
    scd_edth!(op, sc.raised, B2.hmbmb, (-im * ω2) .* B2.hmbmb, ops, -2, 0, a2)         # ð h_m̄m̄
    @inbounds @. sd += -(1 // 2) * B1.psi3 * (op + (cpi + tau0) * B2.hmbmb)
    # T_f : −½ Ψ₃ · (Δ−2μ+μ̄) h_lm̄
    scd_thornprime!(op, B2.hlmb, (-im * ω2) .* B2.hlmb, sc.DR, ops, 2)                 # Δ h_lm̄ (fall 2)
    @inbounds @. sd += -(1 // 2) * B1.psi3 * (op + (-2 * mu + cmu) * B2.hlmb)
    # T_g : − h_lm̄ (Δ Ψ₃)
    scd_thornprime!(op, B1.psi3, (-im * ω1) .* B1.psi3, sc.DR, ops, 2)                 # Δ Ψ₃ (fall 2)
    @inbounds @. sd += -B2.hlmb * op
    # T_h : + ½ h_m̄m̄ (ð Ψ₃)
    scd_edth!(op, sc.raised, B1.psi3, (-im * ω1) .* B1.psi3, ops, -1, -1, a1)          # ð Ψ₃
    @inbounds @. sd += (1 // 2) * B2.hmbmb * op
    # T_i : + 4 π⁽¹⁾ Ψ₃
    @inbounds @. sd += 4 * B2.pi * B1.psi3
    # T_j : − 3 λ⁽¹⁾ Ψ₂
    @inbounds @. sd += -3 * B2.lam * B1.psi2
    return sd
end

"""
    scd_st!(st, B1, B2, ops, sc) -> st

Assemble 𝔰_t (arXiv:2010.00162 eq:def_mathfrak_s): Weyl Ψ₄ from `B1`, H_− metric {h_lm,h_mm,π̄} from
`B2`. Every term is (1,1) in the H_± grading. `st` has GHP type (−1,−2), falloff 1.
"""
function scd_st!(st, B1::SCDBundle{T}, B2::SCDBundle{T}, ops::SCDSourceOps{T}, sc::_SCDScratch{T}) where {T}
    a1, a2 = B1.m_ang, B2.m_ang
    ω1, ω2 = B1.omega, B2.omega
    mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; tau0 = ops.gh.tau_0; ctau = conj.(tau0)
    op = sc.op
    fill!(st, zero(Complex{T}))

    # U_a : − h_lm (Δ+μ+2μ̄) Ψ₄
    scd_thornprime!(op, B1.psi4, (-im * ω1) .* B1.psi4, sc.DR, ops, 1)                 # ΔΨ₄ (fall 1)
    @inbounds @. st += -B2.hlm * (op + (mu + 2 * cmu) * B1.psi4)
    # U_b : + ½ h_mm (ð′ Ψ₄)
    scd_edth_prime!(op, sc.lowered, B1.psi4, (-im * ω1) .* B1.psi4, ops, -2, -2, a1)   # ð′ Ψ₄
    @inbounds @. st += (1 // 2) * B2.hmm * op
    # U_c : + Ψ₄ π̄⁽¹⁾
    @inbounds @. st += B1.psi4 * B2.pibar
    # U_d : + Ψ₄ (− Δ h_lm)
    scd_thornprime!(op, B2.hlm, (-im * ω2) .* B2.hlm, sc.DR, ops, 2)                   # Δ h_lm (fall 2)
    @inbounds @. st += -B1.psi4 * op
    # U_e : + Ψ₄ (ð′−½π−½τ̄) h_mm
    scd_edth_prime!(op, sc.lowered, B2.hmm, (-im * ω2) .* B2.hmm, ops, 2, 0, a2)       # ð′ h_mm
    @inbounds @. st += B1.psi4 * (op + (-(1 // 2) * pi0 - (1 // 2) * ctau) * B2.hmm)
    return st
end

# ---------------------------------------------------------------------------------------------------
# (4) 𝒮 = d₄⁽⁰⁾ 𝔰_d + d₃⁽⁰⁾ 𝔰_t   (the outer background operators)
# ---------------------------------------------------------------------------------------------------

"""
    scd_source!(S, B1, B2, ops; sc) -> (; S, sd, st)

Assemble the full second-order Teukolsky source 𝒮 = d₄⁽⁰⁾𝔰_d + d₃⁽⁰⁾𝔰_t (arXiv:2010.00162
eq:schematic_source) into `S`, from the two sector bundles `B1` (Weyl slot) and `B2` (metric slot).
Returns the assembled source and the intermediate `sd`,`st`. `d₄⁽⁰⁾ = Δ+4μ+μ̄`, `d₃⁽⁰⁾ = ð′+4π−τ̄`
(background operators, GHP-degree 0 in H_±). The composite `∂_T = −i(ω1+ω2)`. The result sits in the
`(a1+a2)` azimuthal sector and has GHP type (−2,−2), falloff 1 (= the second-order Ψ₄ field).
"""
function scd_source!(S, B1::SCDBundle{T}, B2::SCDBundle{T}, ops::SCDSourceOps{T};
                     sc::_SCDScratch{T}=_SCDScratch(T, size(S, 1), size(S, 2))) where {T}
    nx, ny = size(S)
    mtot = B1.m_ang + B2.m_ang
    ωtot = B1.omega + B2.omega
    mu = ops.gh.mu_0; cmu = conj.(mu); pi0 = ops.gh.pi_0; tau0 = ops.gh.tau_0; ctau = conj.(tau0)

    sd = zeros(Complex{T}, nx, ny); st = zeros(Complex{T}, nx, ny)
    scd_sd!(sd, B1, B2, ops, sc)
    scd_st!(st, B1, B2, ops, sc)

    op = sc.op
    fill!(S, zero(Complex{T}))
    # d₄⁽⁰⁾ 𝔰_d = Δ𝔰_d + (4μ+μ̄)𝔰_d      (𝔰_d: fall 1)
    scd_thornprime!(op, sd, (-im * ωtot) .* sd, sc.DR, ops, 1)
    @inbounds @. S += op + (4 * mu + cmu) * sd
    # d₃⁽⁰⁾ 𝔰_t = ð′𝔰_t + (4π−τ̄)𝔰_t      (𝔰_t: type (−1,−2), ð′ at the target sector m_ang=mtot)
    scd_edth_prime!(op, sc.lowered, st, (-im * ωtot) .* st, ops, -1, -2, mtot)
    @inbounds @. S += op + (4 * pi0 - ctau) * st
    return (; S, sd, st)
end

"""
    scd_source_sector(id1, id2, m_pos, sector, ops) -> (; S, sd, st, m_ang, omega)

Convenience driver: assemble 𝒮 in a named `sector` of the QNM tower pair `id1`,`id2` (each an
`hphc_qnm_initialdata` result at positive azimuthal `m_pos`). `sector = :pp` builds the (+m,+m)→(+2m)
source (`B1,B2` = the `+m` bundles); `sector = :pm` the (+m,−m)→(0) source (`B1` = `+m`, `B2` = `−m`).
For the (2,2) tower: `:pp` ⇒ the m=4 sector, `:pm` ⇒ the m=0 sector.
"""
function scd_source_sector(id1, id2, m_pos::Integer, sector::Symbol, ops::SCDSourceOps{T}) where {T}
    B1, B2 = sector === :pp ? (scd_bundle(id1, :pm, m_pos, ops), scd_bundle(id2, :pm, m_pos, ops)) :
             sector === :pm ? (scd_bundle(id1, :pm, m_pos, ops), scd_bundle(id2, :nm, m_pos, ops)) :
             throw(ArgumentError("sector must be :pp (→+2m) or :pm (→0)"))
    nx = ops.Rv |> length; ny = ops.ny
    S = zeros(Complex{T}, nx, ny)
    r = scd_source!(S, B1, B2, ops)
    return (; r.S, r.sd, r.st, m_ang = B1.m_ang + B2.m_ang, omega = B1.omega + B2.omega)
end
