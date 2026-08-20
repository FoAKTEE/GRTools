"""
    PostNewtonian.MST

Symbolic Mano-Suzuki-Takasugi coefficient engine: a truncated-`eps` Laurent algebra (`ES{T}`), the
renormalized-angular-momentum `nu`-solve from the pole-cleared characteristic equation, the MST
coefficients `a_n` as `eps`-series, and the lift of those `eps`-series to PN [`PNSeries`](@ref) (the
standard `eta`-form and the dense `gamma`-FreqRep form). `eps = 2 omega`, so an `eps`-power IS a PN
order. `lambda` (the SWSH eigenvalue) enters as an external `gamma=a*omega` Taylor series, auto-pulled
from `Harmonics.SpinWeightedSpheroidalEigenvalueSeries` for `a != 0`.

Leaf after `Series`; `RadialMode`/`SourceMode` build on the `eps`-series accessor `_mst_eps_series`.
The imaginary unit is carried as the Symbolics symbol `imagunit_sym()` (`𝕚`) so a coefficient stays a
single `Num` (a `Complex{Num}` cannot live in `PNSeries`' `Vector{Num}`); substitute `=> im` to
evaluate numerically.
"""
module MST

using ..Series: PNSeries
using Symbolics: Num
import Symbolics
using Harmonics: SpinWeightedSpheroidalEigenvalueSeries

export ES, MSTCoefficients, imagunit_sym

# ----------------------------------------------------------------------------
# Truncated Laurent eps-series ES{T}: `lead` = lowest stored eps-power, `c[i]` =
# coeff of eps^(lead+i-1). Always truncated so the O[] term is at eps^(ord+1), i.e.
# stored powers run lead .. ord. T is Complex{Rational{BigInt}} (exact, a==0/perfect-
# square kappa) or ComplexF64 (numeric, irrational kappa). Laurent (negative lead) is
# needed for the intermediate continued-fraction ratios at the resonant indices
# n≈-l, -l-1; the physical a_n / characteristic eq come out as proper power series.
# ----------------------------------------------------------------------------

"""
    ES{T}

Truncated Laurent `ε`-series, the internal coefficient carrier of the MST engine: `c[i]` is the
coefficient of `ε^(lead + i - 1)`, always truncated so the `O[]` term sits at `ε^(ord+1)` (stored
powers run `lead .. ord`). `T` is `Complex{Rational{BigInt}}` (exact arithmetic: `a = 0` or
perfect-square `κ`) or `ComplexF64` (numeric, irrational `κ`). The Laurent part (negative `lead`)
is needed for the intermediate continued-fraction ratios at the resonant indices `n ≈ -ℓ, -ℓ-1`;
the physical `aₙ` and the characteristic equation come out as proper power series.
"""
struct ES{T}
    lead::Int
    c::Vector{T}        # c[1] is coeff of eps^lead; length(c) == ord - lead + 1 (>=0)
    ord::Int
end

# normalize: trim leading zeros (raising lead) and trailing beyond ord; an all-zero
# series collapses to lead=ord+1, empty c (pure O[]).
function _es_norm(lead::Int, c::Vector{T}, ord::Int; tol = 0) where {T}
    # drop entries past ord
    hi = ord - lead + 1
    if length(c) > hi
        c = c[1:max(hi, 0)]
    end
    # trim leading (numerically/exactly) zero coeffs
    i = 1
    while i <= length(c) && (iszero(c[i]) || (tol > 0 && abs(c[i]) <= tol))
        i += 1
    end
    if i > length(c)
        return ES{T}(ord + 1, T[], ord)
    end
    ES{T}(lead + (i - 1), c[i:end], ord)
end

_es_iszero(a::ES) = isempty(a.c)

# scalar constant (eps^0)
_es_const(x::T, ord::Int) where {T} = _es_norm(0, T[x], ord)
# series from explicit coeffs c0,c1,... (powers 0,1,...)
_es_from(coeffs::AbstractVector{T}, ord::Int) where {T} = _es_norm(0, collect(T, coeffs), ord)

# add
function _es_add(a::ES{T}, b::ES{T}, ord::Int) where {T}
    (_es_iszero(a)) && return b
    (_es_iszero(b)) && return a
    lead = min(a.lead, b.lead)
    hi = ord
    n = hi - lead + 1
    n <= 0 && return ES{T}(ord + 1, T[], ord)
    out = zeros(T, n)
    @inbounds for (i, x) in enumerate(a.c)
        e = a.lead + i - 1; e > ord && break
        out[e - lead + 1] += x
    end
    @inbounds for (i, x) in enumerate(b.c)
        e = b.lead + i - 1; e > ord && break
        out[e - lead + 1] += x
    end
    _es_norm(lead, out, ord)
end

# negate / scalar-scale
_es_scale(a::ES{T}, s) where {T} = _es_iszero(a) ? a : _es_norm(a.lead, [s * x for x in a.c], a.ord)

# multiply (Cauchy product, truncated at ord)
function _es_mul(a::ES{T}, b::ES{T}, ord::Int) where {T}
    (_es_iszero(a) || _es_iszero(b)) && return ES{T}(ord + 1, T[], ord)
    lead = a.lead + b.lead
    n = ord - lead + 1
    n <= 0 && return ES{T}(ord + 1, T[], ord)
    out = zeros(T, n)
    @inbounds for (i, x) in enumerate(a.c)
        iszero(x) && continue
        ei = a.lead + i - 1
        for (j, y) in enumerate(b.c)
            e = ei + b.lead + j - 1
            e > ord && break
            out[e - lead + 1] += x * y
        end
    end
    _es_norm(lead, out, ord)
end

# reciprocal: 1/(eps^lead u) = eps^-lead u^-1 ; u(0)!=0. Result truncated at ord.
function _es_inv(a::ES{T}, ord::Int) where {T}
    _es_iszero(a) && throw(ArgumentError("_es_inv: division by an all-zero (pure O[]) series"))
    u0 = a.c[1]
    rlead = -a.lead                      # leading power of the inverse
    n = ord - rlead + 1                  # number of inverse coeffs to compute
    n <= 0 && return ES{T}(ord + 1, T[], ord)
    inv0 = inv(u0)
    out = zeros(T, n)
    out[1] = inv0
    @inbounds for k in 1:(n - 1)
        acc = zero(T)
        for j in 1:k
            (j + 1 <= length(a.c)) || continue
            acc += a.c[j + 1] * out[k - j + 1]
        end
        out[k + 1] = -acc * inv0
    end
    _es_norm(rlead, out, ord)
end

_es_div(a::ES{T}, b::ES{T}, ord::Int) where {T} = _es_mul(a, _es_inv(b, ord), ord)

# leading eps-power (== ord+1 if all-zero)
_es_lead(a::ES; tol = 0) = _es_iszero(a) ? (a.ord + 1) : a.lead

# coeff of eps^k (0 if outside stored range)
function _es_coeff(a::ES{T}, k::Int) where {T}
    (_es_iszero(a) || k < a.lead || k - a.lead + 1 > length(a.c)) && return zero(T)
    a.c[k - a.lead + 1]
end

# ----------------------------------------------------------------------------
# MST recurrence coefficients alpha_n, beta_n, gamma_n as eps-series (MST.m 217-227).
# nu is itself an eps-series (l + corrections); q=a numeric; kappa=sqrt(1-a^2),
# tau=(eps - m a)/kappa. lambda is an eps-series (passed in). s,m,n integers.
# We build each coefficient by composing eps-series products/quotients.
# ----------------------------------------------------------------------------

struct _MSTParams{T}
    s::Int
    m::Int
    a::T            # numeric spin (Rational or Float64 wrapped in T's real type)
    kappa::T        # sqrt(1-a^2)
    tau::ES{T}      # eps-series of (eps - m a)/kappa
    lambda::ES{T}   # eps-series of the SWSH eigenvalue
    ord::Int
end

# (off + nu) as an ES given integer/scalar offset `off`
_es_shift(nu::ES{T}, off) where {T} = _es_add(nu, _es_const(T(off), nu.ord), nu.ord)
# (k*nu + off) as ES  (for (2n+2nu+c)-type factors: k=2, off=2n+c)
_es_scale2(nu::ES{T}, k, off) where {T} = _es_add(_es_scale(nu, T(k)), _es_const(T(off), nu.ord), nu.ord)
# eps^p monomial with coefficient `coef`
_es_mono(coef::T, p::Int, ord::Int) where {T} = _es_norm(p, T[coef], ord)

# POLE-CLEARED coefficients: alpha_n,beta_n,gamma_n each have 1/(n+nu) or 1/(n+nu+1)
# poles that hit the resonance n≈-l when nu=l (exactly). The CF ratios only use
# alpha/beta/gamma inside ratios where multiplying num & denom by (n+nu)(n+nu+1)
# clears ALL three poles simultaneously (alpha's 1/(n+nu+1), gamma's 1/(n+nu),
# beta's 1/((n+nu)(n+nu+1))). So we build the cleared
#   ahat_n = alpha_n (n+nu)(n+nu+1),  bhat_n = beta_n (n+nu)(n+nu+1),
#   ghat_n = gamma_n (n+nu)(n+nu+1)
# which are pole-free polynomials/series in eps (no division by a vanishing series).
# The ratio R_n = -ghat_n/(bhat_n + ahat_n R_{n+1}) is identical to the unscaled one.

# (n+nu)(n+nu+1)  (the common clearing factor)
_nnp1(nu::ES{T}, n::Int, ord::Int) where {T} = _es_mul(_es_shift(nu, n), _es_shift(nu, n + 1), ord)

# ahat_n = (I eps kappa (n+nu+1+s+I eps)(n+nu+1+s-I eps)(n+nu+1+I tau)) (n+nu) / (2n+2nu+3)
function _ahat(p::_MSTParams{T}, nu::ES{T}, n::Int) where {T}
    ord = p.ord; I = one(T) * im
    nps = _es_shift(nu, n + 1 + p.s)
    f1 = _es_add(nps, _es_mono(I, 1, ord), ord)
    f2 = _es_add(nps, _es_mono(-I, 1, ord), ord)
    f3 = _es_add(_es_shift(nu, n + 1), _es_scale(p.tau, I), ord)
    num = _es_mul(_es_mul(f1, f2, ord), f3, ord)
    num = _es_mul(_es_mono(I * p.kappa, 1, ord), num, ord)
    num = _es_mul(num, _es_shift(nu, n), ord)                  # × (n+nu) [the (n+nu+1) cancels]
    _es_div(num, _es_scale2(nu, 2, 2n + 3), ord)              # ÷ (2n+2nu+3)
end

# bhat_n = beta_n (n+nu)(n+nu+1)
#        = (-lam - s(s+1) + (n+nu)(n+nu+1) + eps^2 + eps(eps-m a)) (n+nu)(n+nu+1)
#          + eps(eps-m a)(s^2+eps^2)
function _bhat(p::_MSTParams{T}, nu::ES{T}, n::Int) where {T}
    ord = p.ord; s = p.s
    nn = _nnp1(nu, n, ord)
    epseps = _es_mono(one(T), 2, ord)
    eps = _es_mono(one(T), 1, ord)
    epsema = _es_mul(eps, _es_from(T[T(-p.m) * p.a, one(T)], ord), ord)
    s2eps2 = _es_from(T[T(s)^2, zero(T), one(T)], ord)
    poly = _es_add(_es_scale(p.lambda, T(-1)), _es_const(T(-s * (s + 1)), ord), ord)
    poly = _es_add(poly, nn, ord)
    poly = _es_add(poly, _es_add(epseps, epsema, ord), ord)
    _es_add(_es_mul(poly, nn, ord), _es_mul(epsema, s2eps2, ord), ord)
end

# ghat_n = -( I eps kappa (n+nu-s+I eps)(n+nu-s-I eps)(n+nu-I tau) ) (n+nu+1) / (2n+2nu-1)
function _ghat(p::_MSTParams{T}, nu::ES{T}, n::Int) where {T}
    ord = p.ord; I = one(T) * im
    nms = _es_shift(nu, n - p.s)
    f1 = _es_add(nms, _es_mono(I, 1, ord), ord)
    f2 = _es_add(nms, _es_mono(-I, 1, ord), ord)
    f3 = _es_add(_es_shift(nu, n), _es_scale(p.tau, -I), ord)
    num = _es_mul(_es_mul(f1, f2, ord), f3, ord)
    num = _es_mul(_es_mono(I * p.kappa, 1, ord), num, ord)
    num = _es_mul(num, _es_shift(nu, n + 1), ord)             # × (n+nu+1) [the (n+nu) cancels]
    _es_scale(_es_div(num, _es_scale2(nu, 2, 2n - 1), ord), T(-1))
end

# ----------------------------------------------------------------------------
# Continued-fraction ratios (minimal solution), pole-cleared.
#   R_n = f_n/f_{n-1} = -ghat_n / (bhat_n + ahat_n R_{n+1});  seed R_{Ntail+1}=0.
#   L_n = f_n/f_{n+1} = -ahat_n / (bhat_n + ghat_n L_{n-1});  seed L_{-Ntail-1}=0.
# Evaluated as Laurent eps-series. `Ntail` controls CF depth; the recurrence
# converges as a power series at fixed eps-order.
# ----------------------------------------------------------------------------

function _ratio_R(p::_MSTParams{T}, nu::ES{T}, n::Int, Ntail::Int) where {T}
    ord = p.ord
    R = ES{T}(ord + 1, T[], ord)        # seed 0
    for k in Ntail:-1:n
        den = _es_add(_bhat(p, nu, k), _es_mul(_ahat(p, nu, k), R, ord), ord)
        R = _es_div(_es_scale(_ghat(p, nu, k), T(-1)), den, ord)
    end
    R
end

function _ratio_L(p::_MSTParams{T}, nu::ES{T}, n::Int, Ntail::Int) where {T}
    ord = p.ord
    L = ES{T}(ord + 1, T[], ord)        # seed 0
    for k in (-Ntail):n
        den = _es_add(_bhat(p, nu, k), _es_mul(_ghat(p, nu, k), L, ord), ord)
        L = _es_div(_es_scale(_ahat(p, nu, k), T(-1)), den, ord)
    end
    L
end

# ----------------------------------------------------------------------------
# Solve nu = l + sum_{k>=2} nuC[k] eps^k order-by-order from the characteristic
# equation char(nu) = bhat_0 + ahat_0 R_1 + ghat_0 L_{-1} = 0 (pole-cleared).
# At each eps-order p, char_p is AFFINE in the new unknown nuC[p]; we solve it by
# two evaluations (see _solve_nu). Corrections start at eps^2 (PN.wl: nuMST =
# l + sum dnu eps^(i+2)).
# ----------------------------------------------------------------------------

# set the coeff of eps^p in an ES (p within 0..ord), returning a new ES.
function _es_setcoeff(nu::ES{T}, pwr::Int, val::T) where {T}
    ord = nu.ord
    lead = min(_es_iszero(nu) ? pwr : nu.lead, pwr)
    n = ord - lead + 1
    out = zeros(T, n)
    if !_es_iszero(nu)
        @inbounds for (i, x) in enumerate(nu.c)
            e = nu.lead + i - 1; e > ord && break
            out[e - lead + 1] = x
        end
    end
    (pwr >= lead && pwr <= ord) && (out[pwr - lead + 1] = val)
    _es_norm(lead, out, ord)
end

function _solve_nu(p::_MSTParams{T}, l::Int, Ntail::Int) where {T}
    ord = p.ord
    nu = _es_const(T(l), ord)            # nu = l + (corrections, initially 0)
    # corrections start at eps^2 (PN.wl: nuMST = l + sum dnu eps^(i+2)).
    # char_p is AFFINE in nuC[porder]: char_p(x) = B + A x. We must NEVER evaluate at
    # nuC[porder]=0, because then (n+nu) at the resonance n=-l collapses to identically
    # zero and the CF tail degenerates (all-zero denominator). Instead probe at TWO
    # nonzero values x=1,2 (which keep (n+nu) ~ x eps^porder nonzero through the
    # resonance) and solve the affine relation:  A=char_p(2)-char_p(1), B=2 char_p(1)-char_p(2),
    # nuC = -B/A = (2 c1 - c2)/(c2 - c1).
    for porder in 2:ord
        nu1 = _es_setcoeff(nu, porder, _es_coeff(nu, porder) + one(T))
        nu2 = _es_setcoeff(nu, porder, _es_coeff(nu, porder) + T(2))
        c1 = _es_coeff(_char(p, nu1, Ntail), porder)
        c2 = _es_coeff(_char(p, nu2, Ntail), porder)
        A = c2 - c1
        iszero(A) && continue            # no constraint here (e.g. odd orders, a=0)
        B = T(2) * c1 - c2
        nu = _es_setcoeff(nu, porder, _es_coeff(nu, porder) - B / A)
    end
    nu
end

# Characteristic equation (pole-cleared, ×ν(ν+1) at n=0):
#   bhat_0 + ahat_0 R_1 + ghat_0 L_{-1} = 0,  with R_1 = f_1/f_0, L_{-1} = f_{-1}/f_0.
function _char(p::_MSTParams{T}, nu::ES{T}, Ntail::Int) where {T}
    ord = p.ord
    b0 = _bhat(p, nu, 0)
    a0R = _es_mul(_ahat(p, nu, 0), _ratio_R(p, nu, 1, Ntail), ord)
    g0L = _es_mul(_ghat(p, nu, 0), _ratio_L(p, nu, -1, Ntail), ord)
    _es_add(_es_add(b0, a0R, ord), g0L, ord)
end

# ----------------------------------------------------------------------------
# Assemble a_n eps-series from the ratios:  a_0 = 1; a_n = a_{n-1} R_n (n>0);
# a_n = a_{n+1} L_n (n<0).  Returns Dict(n => ES), n in [-nmax, nmax].
# ----------------------------------------------------------------------------

function _solve_a(p::_MSTParams{T}, nu::ES{T}, nmax::Int, Ntail::Int) where {T}
    ord = p.ord
    a = Dict{Int,ES{T}}()
    a[0] = _es_const(one(T), ord)
    for n in 1:nmax
        a[n] = _es_mul(a[n - 1], _ratio_R(p, nu, n, Ntail), ord)
    end
    for n in -1:-1:-nmax
        a[n] = _es_mul(a[n + 1], _ratio_L(p, nu, n, Ntail), ord)
    end
    a
end

# ----------------------------------------------------------------------------
# eps-series (coeff of eps^k) -> PNSeries in :eta with coeffs  c_k * omega^k at
# eta-exponent 3k (eps -> 2 omega eta^3, i.e. eps^k -> (2 omega)^k eta^(3k)).
# `orderEta` is the requested PN order (O[eta]^(orderEta+1)). `ωsym` the omega Num.
# Leading/trailing eta-orders track the eps leading order * 3.
# ----------------------------------------------------------------------------

function _eps_to_pnseries(es::ES{T}, orderEta::Int, ωsym::Num; tol = 0) where {T}
    nmax_eta = orderEta + 1                    # O[eta]^(orderEta+1)
    lead_eps = _es_lead(es; tol = tol)         # lowest nonzero eps-power (negative leads can't occur for physical a_n/nu)
    if lead_eps < 0 || lead_eps > es.ord || 3 * lead_eps >= nmax_eta
        # pure O[] (within requested eta-window): empty coeffs at nmin=nmax=nmax_eta
        return PNSeries(:η, Num[], nmax_eta, nmax_eta, 1)
    end
    nmin_eta = 3 * lead_eps
    coeffs = Num[]
    for k in lead_eps:es.ord
        ek = 3 * k
        ek >= nmax_eta && break
        while nmin_eta + length(coeffs) < ek    # pad zero eta-slots (step 3)
            push!(coeffs, Num(0))
        end
        c = _es_coeff(es, k)
        cval = (tol > 0 && abs(c) <= tol) ? Num(0) : _to_num(c) * (Num(2)^k) * ωsym^k
        push!(coeffs, cval)
    end
    PNSeries(:η, coeffs, nmin_eta, nmax_eta, 1)
end

# eps-series -> FreqRep PNSeries in :gamma: gamma^k coeff = (eps^k coeff)*(2 omega)^k,
# consecutive (step-1, den=1), nmin = eps leading order, O[gamma]^order. This is MMA's
# MSTCoefficientsInternalFreq form (KerrMSTSeries[..,order-1] carries O[eps]^order, then
# ChangeSeriesParameter[#,2 omega] // PowerCounting[#,gamma]); so the O[] is at gamma^order.
function _eps_to_freqseries(es::ES{T}, order::Int, ωsym::Num; tol = 0) where {T}
    nmaxg = order                              # O[gamma]^order
    lead = _es_lead(es; tol = tol)
    if lead < 0 || lead > es.ord || lead >= nmaxg
        return PNSeries(:γ, Num[], nmaxg, nmaxg, 1)
    end
    coeffs = Num[]
    for k in lead:min(es.ord, nmaxg - 1)
        c = _es_coeff(es, k)
        cval = (tol > 0 && abs(c) <= tol) ? Num(0) : _to_num(c) * (Num(2)^k) * ωsym^k
        push!(coeffs, cval)
    end
    PNSeries(:γ, coeffs, lead, nmaxg, 1)
end

# The imaginary unit as a Symbolics symbol (keeps coefficients a single `Num`; a
# `Complex{Num}` cannot be stored in PNSeries' `Vector{Num}`). Substitute `imagunit_sym() => im`
# to evaluate numerically. MMA's InputForm uses `I`; here it is the symbol `𝕚`.
const _IUNIT_SYM = Symbolics.variable(:𝕚)
"The Symbolics symbol used for the imaginary unit in MSTCoefficients output (substitute `=> im` to evaluate)."
imagunit_sym() = _IUNIT_SYM

# Complex{Rational}/ComplexF64 -> Symbolics Num, using `𝕚` for the imaginary part.
function _to_num(c)
    re = real(c); ie = imag(c)
    Num(re) + Num(ie) * _IUNIT_SYM
end

# ----------------------------------------------------------------------------
# Raw eps-series accessor: returns (nu::ES, aeps::Dict{Int,ES}, kappa::T, av::T)
# computing the renormalized angular momentum nu and the MST coefficients a_n as
# exact (or float) eps-series to internal eps-order `epsOrd`. Used by the radial-PN
# builder (slice 3), which needs the numeric eps-coefficients directly. `nkeys`
# bounds the |n| range of a_n returned. The eps-series have eps == 2 omega.
# ----------------------------------------------------------------------------
function _mst_eps_series(s::Integer, l::Integer, m::Integer, a::Real, epsOrd::Integer;
        lambda_series = nothing, nkeys::Integer = epsOrd + 2 * Int(l) + 4, Ntail = :auto)
    s in (0, -1, -2) || throw(ArgumentError("_mst_eps_series supports s in {0,-1,-2}; got s=$s"))
    l >= abs(s) || throw(ArgumentError("need l >= |s|; got l=$l, s=$s"))
    # internal eps-order: resonance at n≈-l raises CF-denominator leading order by ~2l.
    ordInt = Int(epsOrd) + 2 * Int(l) + 4
    Ntailv = Ntail === :auto ? ordInt + Int(l) + 6 : Int(Ntail)
    T, av, kappa = _spin_type(a)
    tau = _es_div(_es_from(T[T(-m) * T(av), one(T)], ordInt), _es_const(kappa, ordInt), ordInt)
    if lambda_series === nothing && !iszero(av)        # a≠0: pull λ(γ=aω) from the SWSH unit
        lambda_series = SpinWeightedSpheroidalEigenvalueSeries(Int(s), Int(l), Int(m); order=ordInt)
    end
    lamcoeffs = lambda_series === nothing ? _default_lambda(s, l, m, ordInt) : lambda_series
    lam = _lambda_eps(lamcoeffs, av, T, ordInt)
    p = _MSTParams{T}(Int(s), Int(m), T(av), kappa, tau, lam, ordInt)
    nu = _solve_nu(p, Int(l), Ntailv)
    aeps = _solve_a(p, nu, Int(nkeys), Ntailv)
    (nu, aeps, kappa, av, ordInt)
end

# ----------------------------------------------------------------------------
# Public API: MSTCoefficients (PN.wl 895-899, standard eta-form).
# ----------------------------------------------------------------------------

"""
    MSTCoefficients(s, l, m, a, ωvar, (pvar, order); FreqRep=false, lambda_series=nothing, Ntail=:auto)
        -> Dict{Symbol,Any}

PN-expanded MST coefficients (Teukolsky`PN`MSTCoefficients). Returns a Dict with
`:νMST => PNSeries` and `:aMST => Dict(n => PNSeries)`. Spec: PN.wl 871-904; MST
recurrence MST.m 217-227.

Two output forms (selected by `FreqRep`):
  * `FreqRep=false` (default) — STANDARD η-form `MSTCoefficients[s,l,m,a,ω,{η,order}]`
    (PN.wl 895). Each entry is a series in `pvar` to real order `order` (`O[pvar]^(order+1)`)
    with the PN scaling `ε -> 2 ω η^3`, so only η^{3k} slots are populated. Keys n over
    `-(epsOrd)-2 .. epsOrd`, epsOrd = ceil(order/3)+1.
  * `FreqRep=true` — FREQUENCY-rep `MSTCoefficients[...,{γ,order},"FreqRep"->True]`
    (PN.wl 900, MSTCoefficientsInternalFreq). The dense ε-series renamed to `pvar=γ`:
    `γ^k` coeff = (ε^k coeff)·(2 ω)^k, consecutive (den=1), `O[γ]^(order+1)`. Keys n over
    `-(order-1)-2 .. order-1`.

Coefficients are `Num` polynomials in `ωvar` and the imaginary-unit symbol `imagunit_sym()`
(substitute `=> im` to evaluate). Spin `a` numeric (Float64 or Rational). `lambda_series`
= SWSH eigenvalue Taylor coeffs `[c0,c1,...]` in γ=a*ω; for `a != 0` it is AUTO-supplied from
`SpinWeightedSpheroidalEigenvalueSeries(s,l,m)` when not given (pass it to override).
"""
function MSTCoefficients(s::Integer, l::Integer, m::Integer, a::Real,
        ωvar::Symbol, porder::Tuple{Symbol,<:Integer};
        FreqRep::Bool = false, lambda_series = nothing, Ntail = :auto)
    pvar, order = porder
    order >= 1 || throw(ArgumentError("order must be >= 1; got $order"))
    s in (0, -1, -2) || throw(ArgumentError("MSTCoefficients port supports s in {0,-1,-2}; got s=$s"))
    l >= abs(s) || throw(ArgumentError("need l >= |s|; got l=$l, s=$s"))

    # OUTPUT eps-order: standard eta-form needs eps up to ceil(order/3) (eps^k -> eta^3k,
    # replsMST: KerrMSTSeries[..,ceil(order/3+1)]); FreqRep needs eps up to `order` (gamma=eps).
    epsOrd = FreqRep ? order : max(2, cld(order, 3) + 1)
    # INTERNAL eps-order: the resonance at n≈-l pushes the leading order of CF
    # denominators UP by ~2l (they cancel through low orders; only eps^(2l+..) survives).
    # Compute at extended order so the leading term is captured, then truncate output.
    ordInt = epsOrd + 2 * Int(l) + 4
    Ntailv = Ntail === :auto ? ordInt + Int(l) + 6 : Int(Ntail)

    # choose exact or float arithmetic for the spin
    T, av, kappa = _spin_type(a)

    # tau = (eps - m a)/kappa  as eps-series
    tau = _es_div(_es_from(T[T(-m) * T(av), one(T)], ordInt), _es_const(kappa, ordInt), ordInt)

    # lambda series in gamma_swsh = a*omega = a*eps/2: substitute gamma -> a eps/2.
    # a≠0: auto-pull the SWSH eigenvalue Taylor series (the λ-input is owned by the SWSH unit).
    if lambda_series === nothing && !iszero(av)
        lambda_series = SpinWeightedSpheroidalEigenvalueSeries(Int(s), Int(l), Int(m); order=ordInt)
    end
    lamcoeffs = lambda_series === nothing ? _default_lambda(s, l, m, ordInt) : lambda_series
    lam = _lambda_eps(lamcoeffs, av, T, ordInt)

    p = _MSTParams{T}(Int(s), Int(m), T(av), kappa, tau, lam, ordInt)

    nu = _solve_nu(p, Int(l), Ntailv)
    # a-key range matches MMA: standard -(epsOrd)-2..epsOrd; FreqRep -(order-1)-2..order-1
    keylo, keyhi = FreqRep ? (-(order - 1) - 2, order - 1) : (-(epsOrd) - 2, epsOrd)
    aeps = _solve_a(p, nu, max(-keylo, keyhi), Ntailv)

    ωsym = Symbolics.variable(ωvar)
    # numeric tolerance for "is this eps-coeff zero" when using Float arithmetic
    tol = T <: Complex{<:AbstractFloat} ? 1e-9 : 0
    emit(es) = FreqRep ? _eps_to_freqseries(es, order, ωsym; tol = tol) :
                         _eps_to_pnseries(es, order, ωsym; tol = tol)
    νps = emit(nu)
    aMSTd = Dict{Int,PNSeries}()
    for n in keylo:keyhi
        es = get(aeps, n, _es_const(zero(T), ordInt))
        aMSTd[n] = emit(es)
    end

    # rename :η coeffs to the requested ηvar if different
    # rename the internal param (:η standard / :γ FreqRep) to the requested `pvar`
    internalvar = FreqRep ? :γ : :η
    if pvar != internalvar
        νps = _rename(νps, pvar)
        for (k, v) in aMSTd
            aMSTd[k] = _rename(v, pvar)
        end
    end
    Dict{Symbol,Any}(:νMST => νps, :aMST => aMSTd)
end

_rename(s::PNSeries, v::Symbol) = PNSeries(v, s.x0, copy(s.coeffs), s.nmin, s.nmax, s.den)

# spin -> (element type, value, kappa). Exact Complex{Rational{BigInt}} is used ONLY
# when kappa = sqrt(1-a^2) is itself rational (a==0 -> kappa=1; or 1-a^2 a perfect
# rational square). Otherwise kappa is irrational (e.g. sqrt(99)/10 at a=1/10) and we
# fall back to ComplexF64 (numeric path; verified by numeric probe to ~1e-10).
function _spin_type(a::Real)
    if a isa Integer || a isa Rational
        ar = Rational{BigInt}(a)
        k2 = 1 - ar^2
        krat = _try_sqrt_rational(k2)
        if krat !== nothing
            T = Complex{Rational{BigInt}}
            return T, ar, T(krat)
        end
    end
    return ComplexF64, Float64(a), ComplexF64(sqrt(1 - Float64(a)^2))
end

# exact rational sqrt of a nonneg rational, or nothing if irrational.
function _try_sqrt_rational(x::Rational{BigInt})
    x < 0 && return nothing
    n = numerator(x); d = denominator(x)
    sn = isqrt(n); sd = isqrt(d)
    (sn * sn == n && sd * sd == d) ? Rational{BigInt}(sn, sd) : nothing
end

# lambda eps-series from gamma-Taylor coeffs `lamc` (γ = a*omega = a*eps/2):
#   lambda(eps) = sum_i lamc[i] (a eps/2)^i, truncated to eps^ord.
function _lambda_eps(lamc, av, ::Type{T}, ord::Int) where {T}
    v = zeros(T, ord + 1)
    half_a = T(av) / 2
    pw = one(T)                 # (a/2)^i tracked, eps-power i
    for i in 0:min(length(lamc) - 1, ord)
        v[i + 1] += T(_rat_or_float(lamc[i + 1], T)) * pw
        pw *= half_a
    end
    _es_from(v, ord)
end

_rat_or_float(x, ::Type{Complex{Rational{BigInt}}}) = Rational{BigInt}(x)
_rat_or_float(x, ::Type{ComplexF64}) = ComplexF64(x)

# default lambda series in γ=a*ω: the leading (spherical) value l(l+1)-s(s+1).
# This is EXACT and complete for a==0 (γ=0, higher coeffs unused). For a!=0 the
# higher γ-coefficients are needed; the caller MUST pass `lambda_series` (the SWSH
# eigenvalue Taylor coeffs), per the spec (λ is owned by the SWSH unit). We do NOT
# pull in the SWSH package here to keep slice-2 dependency-light; a!=0 without a
# supplied series throws below.
_default_lambda(s::Integer, l::Integer, m::Integer, ::Int) =
    Rational{BigInt}[l * (l + 1) - s * (s + 1)]

end # module MST
