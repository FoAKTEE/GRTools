# QNM/qqnm_source_tables.jl — the BAKED second-order QQNM source-coefficient tables a_ij^(ν)/b_ij^(μ):
# the owner's certified tables as {T}-generic rational evaluators in (r, Γ, Γ̄), frozen as a
# GENERATED artifact.  Runtime is pure Julia; the offline generator lives under
# test/fixtures/gen_qqnm_source_tables.jl and re-derives these literals from the .tex + checks sha256.
#
# ─────────────────────────────── PROVENANCE LEDGER (sha256, 2026-07-06) ────────────────────────────
# Source: ref-code/QQNM-coeff/teukolsky-category/coeff/  (owner Ma–Yang tables; F7-3-certified basis)
#   spin2-plus-plus-coefficients.generated.tex        ea1668805fb4b0e5d6596bb86356b76f703012fc0f1b5745140d155f73e62301
#   spin2-plus-plus-coefficients-full.generated.tex   f26e3bea15a42b17a2353a5b769e4d12a4f0fa4ef37b663a6d691eb3aec880dc
#   spin2-plus-minus-sixteen-basis.generated.tex      ccffcda6d124029a39504fcc6da465455cafbe87eb1419324a3e587f9ae8c727
#   spin2-plus-minus-four-basis.generated.tex         bf562422665ba7473dabad991a8c02c04b45ab064d8486d681561d22c1d19c98
#   spin2-minus-minus-zero-proof.generated.tex        15356b9bbdffcd34280b5db844991767018d0d25c67e5203c011c0a56f6d424b
#   spin2-minus-minus-org-cck-slot-proof.generated.tex 3830a58e73a6b91fdaa9249498733ebe805b2b50661bf7a8b8d96892a2bcaacd
#
# ─────────────────────────────── CHANNELS ─────────────────────────────────────────────────────────
#   ++  A_ν = Σ_{0≤i≤j≤6} a_ij^(ν) X_i X_j,  X_i = D₂†ⁱℛ  (Paper I:591 eq:final_source_term;
#       owner spin2-plus-plus).  ν∈{1,2,3,4} ↔ the four angular structures
#       [L₂†S]², [S·L₁†L₂†S], [3ia sinθ·S·L₂†S], [a sinθ·S]².  16-basis (X_0..X_6).
#   +−  B_μ^{+-} = Σ_{0≤i,j≤6} b_{ij,+-}^(μ) X_i^+ X_j^-  (owner:1503 eq:app-pm-Bmu-compact),
#       μ=1..16 ↔ the mixed angular atoms 𝒵_μ (owner:1516 eq:app-pm-sixteen-angular-basis).
#       NONZERO radial support (i,j)∈{(4,0),(4,1),(4,2),(5,0),(5,1),(6,0)} (owner:1550).
#   −−  STRUCTURAL ZERO: Π_{--}S₄⁽²⁾=0 ⇒ Q_src^{--}(r)=0, R_{--}=0 (owner:1585 eq:mm-source-zero;
#       spin2-minus-minus-zero-proof).  No table — [`qqnm_mm_is_zero`](@ref) asserts + cites.
#
# ─────────────────────────────── CONVENTIONS ──────────────────────────────────────────────────────
# Paper Γ = r + ia cosθ (GrTool ζ̄), Γ̄ = r − ia cosθ (GrTool ζ); at a=0 Γ=Γ̄=r.  Under the on-shell
# identity Γ + Γ̄ = 2r the machine r-forms below equal the Appendix-C Γ-forms (test-gated).
# ii = im (the +− tables carry explicit i from the mirror-channel factor).
#
# [TYPO, Paper-I:591]. The printed A_1 (D₂†R)² coefficient 3(Γ̄−15Γ)/(2Γ̄⁵) is wrong; we bake the
# owner's corrected a_11^(1) = 3(2Γ̄−15Γ)/(2Γ̄⁵), which reproduces the explicit Schwarzschild Q(r)
# of Paper-I:1427 to machine precision (the printed value drifts 0.98–1.02).

# ═══════════════════════════════ ++  A_ν  (16-basis, 0≤i≤j≤6) ══════════════════════════════════════
"""
    qqnm_pp_aij(ν, ::Type{C}, Γ, Γb, r) -> Dict{Tuple{Int,Int},C}

The nonzero scalar coefficients `a_ij^(ν)` of the ++ second-order source `A_ν = Σ a_ij^(ν) X_i X_j`
(`X_i = D₂†ⁱℛ`), `ν∈{1,2,3,4}`, as {T}-generic rationals in (Γ, Γ̄, r).  Keys are `(i,j)` with
`0≤i≤j≤6`.  Owner spin2-plus-plus table (Paper I:591); see the module PROVENANCE + [TYPO-FOUND].
"""
function qqnm_pp_aij(ν::Integer, ::Type{C}, Γ, Γb, r) where {C}
    Γ = C(Γ); Γb = C(Γb); r = C(r)
    d = Dict{Tuple{Int,Int},C}()
    if ν == 1
        d[(0,0)] = 45*(Γb-2r)/(2*Γb^7);  d[(0,1)] = (90r-48*Γb)/Γb^6;   d[(0,2)] = (45*Γb-84r)/(2*Γb^5)
        d[(0,3)] = -6*(Γb-2r)/Γb^4;      d[(0,4)] = (7*Γb-18r)/(8*Γb^3); d[(0,5)] = r/(4*Γb^2)
        d[(1,1)] = (51*Γb-90r)/(2*Γb^5); d[(1,2)] = -6*(4*Γb-7r)/Γb^4;   d[(1,3)] = (13*Γb-24r)/(2*Γb^3)
        d[(1,4)] = (9r-4*Γb)/(4*Γb^2);   d[(1,5)] = -Γ/(8*Γb);           d[(2,2)] = (45*Γb-78r)/(8*Γb^3)
        d[(2,3)] = (11r-6*Γb)/(2*Γb^2);  d[(2,4)] = -Γ/(2*Γb);           d[(3,3)] = -3*Γ/(8*Γb)
    elseif ν == 2
        d[(0,1)] = -3/Γb^5;              d[(0,2)] = 3r/Γb^5;             d[(0,3)] = 3*(Γb-2r)/(2*Γb^4)
        d[(0,4)] = (3r-2*Γb)/(2*Γb^3);   d[(0,5)] = (3*Γb-4r)/(8*Γb^2);  d[(0,6)] = Γ/(16*Γb)
        d[(1,1)] = 3/Γb^4;               d[(1,2)] = -3*(Γb+2r)/(2*Γb^4); d[(1,3)] = -(Γb-3r)/Γb^3
        d[(1,4)] = (7*Γb-12r)/(8*Γb^2);  d[(1,5)] = Γ/(4*Γb);            d[(2,2)] = 3r/(2*Γb^3)
        d[(2,3)] = (3*Γb-8r)/(4*Γb^2);   d[(2,4)] = 7*Γ/(16*Γb);         d[(3,3)] = Γ/(4*Γb)
    elseif ν == 3
        d[(0,1)] = 5*(5*Γb-6r)/Γb^7;     d[(0,2)] = -5*(3*Γb-4r)/Γb^6;   d[(0,3)] = (7*Γb-12r)/(2*Γb^5)
        d[(0,4)] = -(Γb-6r)/(6*Γb^4);    d[(0,5)] = -(3*Γb+2r)/(24*Γb^3);d[(0,6)] = 1/(24*Γb)
        d[(1,1)] = (30r-26*Γb)/Γb^6;     d[(1,2)] = (57*Γb-68r)/(2*Γb^5);d[(1,3)] = (30r-23*Γb)/(3*Γb^4)
        d[(1,4)] = (25*Γb-42r)/(24*Γb^3);d[(1,5)] = r/(6*Γb^2);          d[(2,2)] = -3*(5*Γb-6r)/(2*Γb^4)
        d[(2,3)] = (45*Γb-58r)/(12*Γb^3);d[(2,4)] = (16r-9*Γb)/(24*Γb^2);d[(3,3)] = (3r-2*Γb)/(6*Γb^2)
    elseif ν == 4
        d[(1,1)] = 15*(3*Γ-4*Γb)/(2*Γb^7); d[(1,2)] = (72*Γb-60r)/Γb^6;  d[(1,3)] = -9*(5*Γb-4r)/(2*Γb^5)
        d[(1,4)] = (4*Γb-3r)/Γb^4;         d[(1,5)] = (Γ-2*Γb)/(8*Γb^3); d[(2,2)] = -9*(5*Γb-4r)/(2*Γb^5)
        d[(2,3)] = (12*Γb-9r)/Γb^4;        d[(2,4)] = (2r-3*Γb)/(2*Γb^3);d[(3,3)] = (6r-9*Γb)/(8*Γb^3)
    else
        throw(ArgumentError("qqnm_pp_aij: ν must be 1..4 (got $ν)"))
    end
    d
end

# Appendix-C Γ-forms for the F7-3b named-entry gate (a DISTINCT string, equal to the r-form only under
# Γ+Γ̄=2r).  A representative spanning A_1..A_4 incl. the [TYPO-FOUND] a_11^(1) and the flagship a_06^(2).
"""
    qqnm_pp_aij_gammaform(name) -> (Γ,Γb,r)->value

Hand-written Γ-form (Paper I:591 / Appendix-C) of a named ++ entry, keyed by `"a_ij^(ν)"`.  Distinct
string from [`qqnm_pp_aij`](@ref)'s r-form; the two agree only under the on-shell `Γ+Γ̄=2r` (gate
G-Q3-T1).  `a_11^(1)` uses the OWNER-CORRECTED `3(2Γ̄−15Γ)/(2Γ̄⁵)` (see [TYPO-FOUND]).
"""
function qqnm_pp_aij_gammaform(name::AbstractString)
    tbl = Dict{String,Function}(
        "a_00^(1)" => (Γ,Γb,r)-> -45*Γ/(2*Γb^7),
        "a_01^(1)" => (Γ,Γb,r)-> 3*(15*Γ-Γb)/Γb^6,
        "a_11^(1)" => (Γ,Γb,r)-> 3*(2*Γb-15*Γ)/(2*Γb^5),      # owner-corrected (Paper printed Γ̄−15Γ)
        "a_02^(1)" => (Γ,Γb,r)-> (3*Γb-42*Γ)/(2*Γb^5),
        "a_06^(2)" => (Γ,Γb,r)-> Γ/(16*Γb),                    # flagship 6th-order
        "a_03^(2)" => (Γ,Γb,r)-> -3*Γ/(2*Γb^4),
        "a_33^(2)" => (Γ,Γb,r)-> Γ/(4*Γb),
        "a_04^(3)" => (Γ,Γb,r)-> -(Γb-6r)/(6*Γb^4),
        "a_11^(4)" => (Γ,Γb,r)-> 15*(3*Γ-4*Γb)/(2*Γb^7),
    )
    haskey(tbl, name) || throw(ArgumentError("qqnm_pp_aij_gammaform: unknown entry $name"))
    tbl[name]
end

# ═══════════════════════════════ +−  B_μ  (16-basis, mixed) ════════════════════════════════════════
# Nonzero radial support (owner:1550): (i,j)∈{(4,0),(4,1),(4,2),(5,0),(5,1),(6,0)}.
const QQNM_PM_SUPPORT = ((4,0),(4,1),(4,2),(5,0),(5,1),(6,0))

"""
    qqnm_pm_bij(μ, ::Type{C}, Γ, Γb, r) -> Dict{Tuple{Int,Int},C}

The nonzero mixed coefficients `b_{ij,+-}^(μ)` of the +− source `B_μ^{+-} = Σ b_{ij,+-}^(μ) X_i^+ X_j^-`
(owner:1503, eq:app-pm-Bmu-compact), `μ=1..16` for the mixed angular atom `𝒵_μ`.  {T}/Complex-generic
rationals in (Γ, Γ̄, r) — the explicit `im` carries the mirror-channel factor 𝔠_-.  `(i,j)` keys lie
in [`QQNM_PM_SUPPORT`](@ref).
"""
function qqnm_pm_bij(μ::Integer, ::Type{C}, Γ, Γb, r) where {C}
    Γ = C(Γ); Γb = C(Γb); r = C(r); ii = C(im)
    d = Dict{Tuple{Int,Int},C}()
    if μ == 1
        d[(4,1)] = -1/(4*Γ)
    elseif μ == 2
        d[(4,0)] = (4r-3*Γb)/(8*Γ^3);  d[(4,1)] = 3r/(8*Γ^2);  d[(5,0)] = r/(4*Γ^2)
    elseif μ == 3
        d[(4,1)] = -1/(8*Γ)
    elseif μ == 4
        d[(4,0)] = (Γ*Γb+3*Γ*r-2*Γb*r)/(8*Γ^3*Γb);  d[(4,1)] = r/(8*Γ^2);  d[(5,0)] = r/(8*Γ^2)
    elseif μ == 5
        d[(4,1)] = (1//1)*(13/Γ - 9/Γb)/16
    elseif μ == 6
        d[(4,0)] = -(-4*Γb^2+Γb*r+6*r^2)/(8*Γ^3*Γb);  d[(4,1)] = (7*Γ-4*Γb)/(16*Γ^2)
        d[(4,2)] = Γb/(8*Γ);  d[(5,0)] = -3*(Γb-r)/(8*Γ^2);  d[(5,1)] = 3*Γb/(16*Γ);  d[(6,0)] = Γb/(16*Γ)
    elseif μ == 7
        d[(4,1)] = ii*(Γ-Γb)/(8*Γ)
    elseif μ == 8
        d[(4,1)] = ii*(3*Γ-2*Γb)/(8*Γ^2*Γb)
    elseif μ == 9
        d[(4,1)] = ii*(Γb/Γ - 1)/8
    elseif μ == 10
        d[(4,1)] = ii*(Γb-3r)/(2*Γ^2*Γb)
    elseif μ == 11
        d[(4,0)] = -ii*(2*Γb^2-5*Γb*r+3*r^2)/(8*Γ^2*Γb);  d[(4,1)] = ii*r*(Γ-Γb)/(8*Γ^2)
        d[(5,0)] = -ii*r*(Γb-r)/(8*Γ^2)
    elseif μ == 12
        d[(4,1)] = -ii*(7*Γb^2-21*Γb*r+6*r^2)/(8*Γ^3*Γb);  d[(4,2)] = -ii/(4*Γ)
        d[(5,0)] = ii*(3*Γ+2r)/(8*Γ^3);  d[(5,1)] = -3*ii/(8*Γ);  d[(6,0)] = -ii/(8*Γ)
    elseif μ == 13
        d[(4,1)] = ii*(Γb-r)/(8*Γ)
    elseif μ == 14
        d[(4,1)] = ii*(Γb-3r)/(4*Γ^2*Γb)
    elseif μ == 15
        d[(4,1)] = -(Γb-r)^2/(4*Γ)
    elseif μ == 16
        d[(4,1)] = (-6*Γ^2+3*Γ*Γb+4*Γb^2)/(8*Γ^3*Γb^2)
    else
        throw(ArgumentError("qqnm_pm_bij: μ must be 1..16 (got $μ)"))
    end
    d
end

# The mixed angular atoms 𝒵_μ (owner:1516). Symbolic metadata for the four-generator reduction
# (owner derive-angular-four-basis.m): each 𝒵_μ is a word in {S_A, Y_0, Y_1, Y_2, ℒ_0, s_a^k}.
# Stored as a NamedTuple per μ for the +− angular assembler (Deliverable B).  s_a = a sinθ, ℒ_0 the
# local angular derivative; Y_0,Y_1,Y_2 the reflected H_- CCK angular slots.
const QQNM_PM_ATOMS = (
    (μ=1,  word="L0(S_A L0 Y0)",  sapow=0),  (μ=2,  word="L0(S_A Y1)",    sapow=0),
    (μ=3,  word="L0(Y0 L0 S_A)",  sapow=0),  (μ=4,  word="S_A L0 Y1",     sapow=0),
    (μ=5,  word="S_A Y0",         sapow=0),  (μ=6,  word="S_A Y2",         sapow=0),
    (μ=7,  word="L0(S_A Y0)",     sapow=-1), (μ=8,  word="L0(S_A Y0)",     sapow=1),
    (μ=9,  word="S_A L0 Y0",      sapow=-1), (μ=10, word="S_A L0 Y0",      sapow=1),
    (μ=11, word="S_A Y1",         sapow=-1), (μ=12, word="S_A Y1",         sapow=1),
    (μ=13, word="Y0 L0 S_A",      sapow=-1), (μ=14, word="Y0 L0 S_A",      sapow=1),
    (μ=15, word="S_A Y0",         sapow=-2), (μ=16, word="S_A Y0",         sapow=2),
)

# ═══════════════════════════════ −−  structural zero ══════════════════════════════════════════════
"""
    qqnm_mm_is_zero() -> Bool

The −− channel is the EXACT structural zero: `Π_{--}S₄⁽²⁾=0` (supp_{H±}(S₄⁽²⁾)⊂{(2,0),(1,1)}, no
`H_-²` monomial), hence `Q_src^{--}(r)=0` and `R_{--}=0` before any angular projection or normal
ordering (owner:1571 eq:mm-projector-zero, :1585 eq:mm-source-zero; spin2-minus-minus-zero-proof).
No coefficient table exists for this channel.  Returns `true` (the certified assertion).
"""
qqnm_mm_is_zero() = true
