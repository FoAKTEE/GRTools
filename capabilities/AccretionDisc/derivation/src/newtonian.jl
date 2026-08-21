# newtonian.jl — the PAPER-FACING gate of route R3: the Newtonian limit of the derived relativistic
# entries must reproduce Eq. (50) of arXiv:2608.17003 (main.tex M:623-643) at leading order.
#
# This is NOT an oracle comparison: (50) is transcribed here from the paper's own TeX and rebuilt in
# the series ring, and it is the only place in stages 1-3 where route R3 is confronted with a printed
# result.  Back-substitution (stage2.jl) certifies that the derived system SOLVES the perturbed
# equations; this certifies that it is the SAME system the paper writes down — i.e. that the
# master-variable conventions (48),(49), the sign of ω_d (51) and the Fourier convention (45a) are the
# paper's, not merely self-consistent.
#
# Scaling (the standard 1/c counting, s² ~ 1/c²; the recipe route R2 also used):
#   ν → s²ν₁,  ψ → log r + s²ψ₁,  μ → s²μ₁,  ω̃ → s³w₁,  v → s v₁,  ω_d → s ω_{d1},
#   p → s²p₁,  ε → Σ,  Γ₁ = O(1),  ĥ_tt → s²ĥ_tt1  (Φ̂ = −ĥ_tt/2, Eq. (52)),  ĥ_tφ → s³,
#   ĥ_φφ → s²r²,  ĥ_rr, ĥ_zz → s²,  y0 = δÛ → s²Y0,  y1 = −iδû^r → s Y1,
#   ĥ_tr = ĥ_rφ = 0 (a circular equatorial orbit, cf. `MasterInputs`).
# e^ν, e^ψ, e^μ and W = 1/√(1−v²) are built by EXACT series exponentiation (`sexp`, `sinvsqrt1m`), so
# nothing is truncated by hand.

@variables rN nu1 nu1p psi1 psi1p mu1 mu1p w1 w1p v1 v1p p1 p1p Sig Sigp wd1
@variables htt1 htt1p htp1 htp1p hpp1 hpp1r hrr1 hrr1r hzz1 hzz1r Y0 Y1

"The Newtonian coefficients that a gate point must bind (all O(1) rationals)."
const NEWTONIAN_SYMBOLS = [rN, nu1, nu1p, psi1, psi1p, mu1, mu1p, w1, w1p, v1, v1p, p1, p1p,
                           Sig, Sigp, wd1, htt1, htt1p, htp1, htp1p, hpp1, hpp1r,
                           hrr1, hrr1r, hzz1, hzz1r, Y0, Y1]

"""
    newtonian_point(rng) -> Dict

Random exact-rational values for every Newtonian coefficient, for `Γ₁` and for `m`.
"""
function newtonian_point(rng::AbstractRNG)
    pt = Dict{Num,RatB}(x => _rand_signed(rng) for x in NEWTONIAN_SYMBOLS)
    pt[rN] = _rand_rat(rng) + 2                 # r > 0, bounded away from 0
    pt[Sig] = _rand_rat(rng)
    pt[p1] = _rand_rat(rng)
    pt[v1] = _rand_rat(rng)
    pt[Gamma1] = _rand_rat(rng)
    pt[m] = RatB(rand(rng, 1:9))
    pt
end

"""
    newtonian_env(pt) -> Dict

Bind every symbol of the derivation alphabet to its Laurent series in `s` at the point `pt`.
"""
function newtonian_env(pt)
    g(x) = pt[x]
    env = Dict{Any,LSeries}()
    put!(sym, val) = (env[Symbolics.value(sym)] = val)
    put!(enu, sexp(series_mono(g(nu1), 2)))
    put!(nu_r, series_mono(g(nu1p), 2))
    put!(epsi, series_const(g(rN)) * sexp(series_mono(g(psi1), 2)))
    put!(psi_r, series_const(inv(g(rN))) + series_mono(g(psi1p), 2))
    put!(emu, sexp(series_mono(g(mu1), 2)))
    put!(mu_r, series_mono(g(mu1p), 2))
    put!(omt, series_mono(g(w1), 3))
    put!(omt_r, series_mono(g(w1p), 3))
    put!(v, series_mono(g(v1), 1))
    put!(v_r, series_mono(g(v1p), 1))
    put!(W, sinvsqrt1m(series_mono(g(v1)^2, 2)))
    put!(p, series_mono(g(p1), 2))
    put!(p_r, series_mono(g(p1p), 2))
    put!(eps, series_const(g(Sig)))
    put!(eps_r, series_const(g(Sigp)))
    put!(Gamma1, series_const(pt[Gamma1]))
    put!(m, series_const(pt[m]))
    put!(omega_d, series_mono(g(wd1), 1))
    put!(y0, series_mono(g(Y0), 2))
    put!(y1, series_mono(g(Y1), 1))
    put!(htt, series_mono(g(htt1), 2))
    put!(htt_r, series_mono(g(htt1p), 2))
    put!(htphi, series_mono(g(htp1), 3))
    put!(htphi_r, series_mono(g(htp1p), 3))
    put!(hphiphi, series_mono(g(rN)^2 * g(hpp1), 2))
    put!(hphiphi_r, series_mono(g(rN)^2 * g(hpp1r), 2))
    put!(hrr, series_mono(g(hrr1), 2))
    put!(hrr_r, series_mono(g(hrr1r), 2))
    put!(hzz, series_mono(g(hzz1), 2))
    put!(hzz_r, series_mono(g(hzz1r), 2))
    for x in (htr, htr_r, hrphi, hrphi_r)
        put!(x, series_zero())
    end
    env
end

"""
    ref50(pt, env) -> NamedTuple

Eq. (50) of the paper, transcribed from `ref-paper/arXiv-2608.17003/src/main.tex` M:623-643 and
evaluated in the series ring:

    δu_φ  = mδÛ/(rω_d) + mΦ̂/(rω_d) − iδû^rκ²/(2ω_dΩ)            (50a)
    C₀₀   = 2mΩ/(rω_d) + p′/(Γ₁p) − Σ′/Σ                          (50b)
    C₀₁   = −ω_d + [κ² + p′(−Σp′/(Γ₁p) + Σ′)/Σ²]/ω_d              (50c)
    C₁₀   = −m²/(r²ω_d) + ω_dΣ/(Γ₁p)                              (50d)
    C₁₁   = −1/r − mκ²/(2rω_dΩ) − p′/(Γ₁p)                        (50e)
    S₀    = 2mΦ̂Ω/(rω_d) − Φ̂′                                      (50f)
    S₁    = −m²Φ̂/(r²ω_d)                                          (50g)

with Σ = ε, Φ̂ = −ĥ_tt/2 (52), κ² = 4Ω² + 2rΩΩ′ the epicyclic frequency (M:645) built from the SAME
Ω = ω̃ + ve^{ν−ψ} (10) as the relativistic side, and δÛ = y0, −iδû^r = y1 so that the last term of
(50a) is +y1κ²/(2ω_dΩ).
"""
function ref50(pt, env)
    S(x) = series_const(pt[x])
    Ω = eval_series(Ω_expr, env)
    Ωp = eval_series(Ω_r_expr, env)
    r = S(rN)
    wd = series_mono(pt[wd1], 1)
    pv = series_mono(pt[p1], 2)
    pp = series_mono(pt[p1p], 2)
    Σ = S(Sig)
    Σp = S(Sigp)
    Γ = series_const(pt[Gamma1])
    mm = series_const(pt[m])
    Φ = series_mono(-pt[htt1] // 2, 2)
    Φp = series_mono(-pt[htt1p] // 2, 2)
    κ2 = 4 * Ω^2 + 2 * r * Ω * Ωp
    (duphi_low = mm * series_mono(pt[Y0], 2) / (r * wd) + mm * Φ / (r * wd) +
                 series_mono(pt[Y1], 1) * κ2 / (2 * wd * Ω),
     C00 = 2 * mm * Ω / (r * wd) + pp / (Γ * pv) - Σp / Σ,
     C01 = -wd + (κ2 + pp * (-Σ * pp / (Γ * pv) + Σp) / Σ^2) / wd,
     C10 = -(mm^2) / (r^2 * wd) + wd * Σ / (Γ * pv),
     C11 = -inv(r) - mm * κ2 / (2 * r * wd * Ω) - pp / (Γ * pv),
     S0 = 2 * mm * Φ * Ω / (r * wd) - Φp,
     S1 = -(mm^2) * Φ / (r^2 * wd))
end

"Leading s-order of each entry, forced by y0 ~ s², y1 ~ s and dy/dr ~ y/r."
const NEWTONIAN_ORDER = (C00 = 0, C01 = 1, C10 = -1, C11 = 0, S0 = 2, S1 = 1, duphi_low = 1)

_lo(x::LSeries) = lorder(x)
function _pair_order(a, b)
    # order of (complex R3 entry) − (real reference): compare real and imaginary parts
    oa = lorder(a)
    ob = lorder(b)
    oa === nothing && return ob
    ob === nothing && return oa
    min(oa, ob)
end

"""
    newtonian_certificates(; npts=6, seed=2608, verbose=true)

Per entry: (a) the scaled R3 expression starts EXACTLY at the expected order (nothing more singular),
and (b) R3 − (50) starts at least one order higher, i.e. the leading Newtonian coefficients agree.
Also re-derives ERRATUM E-M1 independently: the covariant δû_φ matches r × (50a), NOT (50a).
"""
function newtonian_certificates(; npts::Int = 6, seed::Int = 2608, verbose::Bool = true)
    rng = MersenneTwister(seed)
    r3 = (C00 = MASTER_ONSHELL[1], C01 = MASTER_ONSHELL[2], C10 = MASTER_ONSHELL[3],
          C11 = MASTER_ONSHELL[4], S0 = MASTER_ONSHELL[5], S1 = MASTER_ONSHELL[6],
          duphi_low = ELIMINATED_ONSHELL.duphi_low)
    ks = keys(r3)
    okorder = Dict(k => true for k in ks)
    okmatch = Dict(k => true for k in ks)
    gap = Dict(k => 99 for k in ks)
    em1_r, em1_bare = true, true
    for _ in 1:npts
        pt = newtonian_point(rng)
        env = newtonian_env(pt)
        ref = ref50(pt, env)
        for k in ks
            a = eval_series(r3[k], env)
            are, aim = a isa LSeries ? (a, series_zero()) : a
            b = ref[k]
            expr_order(r3[k], env) == NEWTONIAN_ORDER[k] || (okorder[k] = false)
            d = _pair_order(are - b, aim)
            (d !== nothing && d < NEWTONIAN_ORDER[k] + 1) && (okmatch[k] = false)
            d !== nothing && (gap[k] = min(gap[k], d - NEWTONIAN_ORDER[k]))
            if k === :duphi_low
                dr = _pair_order(are - series_const(pt[rN]) * b, aim)
                (dr !== nothing && dr < 2) && (em1_r = false)
                (d === nothing || d >= 2) && (em1_bare = false)
            end
        end
    end
    C = Certificate[]
    for k in ks
        push!(C, certify_bool("Newtonian order of $(k) is s^$(NEWTONIAN_ORDER[k])", "(50)",
                              okorder[k], "no term more singular than the leading one";
                              verbose, modality = :exact_series, npoints = npts))
        k === :duphi_low && continue
        push!(C, certify_bool("Newtonian limit  $(k) ≡ (50)", "(50)", okmatch[k],
                              "R3 − (50) starts $(gap[k]) order(s) above leading";
                              verbose, modality = :exact_series, npoints = npts))
    end
    push!(C, certify_bool("E-M1: δû_φ (covariant) ≡ r × (50a), NOT (50a)", "(50a)",
                          em1_r && em1_bare,
                          em1_r && em1_bare ?
                          "matches r×(50a); does NOT match (50a) as printed — (50a) is the orthonormal component" :
                          "unexpected (matches_r=$(em1_r), differs_from_printed=$(em1_bare))";
                          verbose, modality = :exact_series, npoints = npts))
    C
end
