# completion.jl — the completion sector + the metric assembly / gauge-match layer (LY:872–968).
#
# For the monopole (mass-shift) shell the only nonzero completion field is γ^{(1)} (=0 here, real
# trace); the dipole/angular-momentum completion (α,β,τ,κ,σ,ε,ρ,Ψ0,Ψ1) is [FUTURE] (all zero for the
# spherically symmetric monopole, LY:333). This file provides the pieces the end-to-end gate needs:
#   (§3.6) the traceful-ORG assembly h_ab → h_μν (LY:890–897),
#   (ANCHOR-0) the exact linearized Frauendiener metric 𝔥_μν (LY:670–685),
#   (ANCHOR-14) the ξ_t gauge generator + its defining ODE RHS (LY:933–967).

# ================================================================================================
# (§3.6) traceful-ORG assembly  h_ab^{(1)} → h_μν^{(1)}   (LY:890–897)
# ================================================================================================
# The traceful-ORG condition h_μν n^ν = 0 + sphericity give (LY:893–895)
#   h_tt = f h_tr = f² h_rr = (f²/4) h_ll ,     h_φφ = sin²θ h_θθ = r² sin²θ h_mm̄ .
# With h_ll = numer/f² (HllField) the f² cancels in h_tt, leaving Laurents.

"h_tt^{(1)}(r) = (f²/4) h_ll = numer_hll(r)/4  (LY:893)."
recon_htt(h::HllField, r) = (r < h.R ? leval(h.numer_in, r) : leval(h.numer_out, r)) / 4
"h_rr^{(1)}(r) = h_ll/4  (LY:893)."
recon_hrr(h::HllField, r) = hll_eval(h, r) / 4
"h_tr^{(1)}(r) = (f/4) h_ll  (LY:893)."
recon_htr(h::HllField{T}, r) where {T} = (1 - 2h.M / r) / 4 * hll_eval(h, r)
"h_φφ^{(1)}(r,θ)/sin²θ = r² h_mm̄(r)  (LY:895)."
recon_hφφ_over_sin2(h_mm̄::ShellProfile, r) = r^2 * peval(h_mm̄, r)

# ================================================================================================
# (ANCHOR-0) the EXACT linearized Frauendiener metric 𝔥_μν^{(1)}  (LY:670–685) — the match target.
# ================================================================================================
"""
    exact_htt(M, R, r) -> value

`𝔥_tt^{(1)}` of LY:672–679: `2M/r` outside, `(2M/Rf(R))(1−3M/r+MR/r²)` inside — the linearization of
the exact Frauendiener thin-shell metric (independent of the transport chain; the G-E6 target).
"""
function exact_htt(M, R, r)
    T = float(promote_type(typeof(M), typeof(R), typeof(r)))
    Mt, Rt, rt = T(M), T(R), T(r)
    fR = 1 - 2Mt / Rt
    rt < Rt ? 2Mt / (Rt * fR) * (1 - 3Mt / rt + Mt * Rt / rt^2) : 2Mt / rt
end
"`𝔥_rr^{(1)} = 𝔥_tt/f²` (LY:673)."
function exact_hrr(M, R, r)
    fr = 1 - 2 * float(promote_type(typeof(M), typeof(R), typeof(r)))(M) / float(r)
    exact_htt(M, R, r) / fr^2
end
"`𝔥_φφ^{(1)}/sin²θ` of LY:680–683: `(2M/Rf(R))r(r−R)` inside, `0` outside."
function exact_hφφ_over_sin2(M, R, r)
    T = float(promote_type(typeof(M), typeof(R), typeof(r)))
    Mt, Rt, rt = T(M), T(R), T(r)
    fR = 1 - 2Mt / Rt
    rt < Rt ? 2Mt / (Rt * fR) * rt * (rt - Rt) : zero(T)
end

# ================================================================================================
# (ANCHOR-14) the ξ_t gauge generator (LY:941–967) + its defining ODE RHS (LY:937/954).
# ================================================================================================
"""
    gauge_ode_rhs(M, R, r) -> value

The RHS of the ξ_t gauge ODE `(∂_r − f′/f)ξ_t^{(1)} = ·` (LY:937 outside, LY:954 inside):
`2M/(rf)` outside, `(2M/Rf(R)f)(1−3M/r+MR/r²)` inside. By construction this EQUALS the reconstructed
`h_tr^{(1)}` — so the single residual `h_tr` is exactly removable by ξ_t, sending `h_μν → 𝔥_μν`
(LY:968). The G-E6 gauge gate asserts `recon_htr == gauge_ode_rhs`.
"""
function gauge_ode_rhs(M, R, r)
    T = float(promote_type(typeof(M), typeof(R), typeof(r)))
    Mt, Rt, rt = T(M), T(R), T(r)
    fr = 1 - 2Mt / rt
    fR = 1 - 2Mt / Rt
    rt < Rt ? 2Mt / (Rt * fR * fr) * (1 - 3Mt / rt + Mt * Rt / rt^2) : 2Mt / (rt * fr)
end

"""
    gauge_xi_t(M, R, r; Ctp=0, Ctm=0) -> value

The closed-form gauge generator `ξ_t^{(1)}` (LY:941–946 outside, LY:959–967 inside):

    ξ_t⁺ = f[C⁺ + 2M log(r/2M−1)] − 4M²/r,
    ξ_t⁻ = f[C⁻ + (2M/Rf(R))(r + M log(r/2M−1))] − 2M²/r,

with the residual-gauge constants `C_t^±` free (**[AXIOM-E]**, LY:947/968 default 0). Provided so the
gate can independently confirm `(∂_r − f′/f)ξ_t = gauge_ode_rhs` (the generator integrates the RHS).
"""
function gauge_xi_t(M, R, r; Ctp = 0, Ctm = 0)
    T = float(promote_type(typeof(M), typeof(R), typeof(r)))
    Mt, Rt, rt = T(M), T(R), T(r)
    fr = 1 - 2Mt / rt
    fR = 1 - 2Mt / Rt
    lg = log(rt / (2Mt) - 1)
    rt < Rt ? fr * (T(Ctm) + 2Mt / (Rt * fR) * (rt + Mt * lg)) - 2Mt^2 / rt :
              fr * (T(Ctp) + 2Mt * lg) - 4Mt^2 / rt
end
