"""
    PostNewtonian.RadialMode

η-graded symbolic radial Coulomb/MST solution and the asymptotic mode amplitudes. The minimal/maximal
Teukolsky radial solutions are built as [`PNSeries`](@ref) in `η` whose coefficients are `Num`
polynomials in `r`, `ω` and the imaginary unit [`imagunit_sym()`](@ref); `ν` (=νMST) and `a_n` (=aMST)
come from [`MST`](@ref)'s `_mst_eps_series` as `ε`-series. Adds the transcendental series ops
(`exp`/`log1p`/power-by-series/Pochhammer/Gamma-ratio/Hurwitz-ζ/sin π) the Coulomb construction needs.

Public: [`TeukolskyRadialPN`](@ref) (the In/Up radial pair) and the amplitude builders
`_A_amp`/`_K_nu`/`_K_ratio`/`_B_trans`/`_C_trans`/`_B_inc` consumed by `SourceMode`.
"""
module RadialMode

using ..Series: PNSeries
import ..Series
import ..MST
using ..MST: imagunit_sym
using Symbolics: Num, Differential
import Symbolics
import SpecialFunctions

export TeukolskyRadialPN

# ---------------------------------------------------------------------------
# eta-graded series with Num coefficients (reuses the slice-1 PNSeries type, den=1).
# A `PNSeries(:η, coeffs, nmin, nmax, 1)` with `coeffs[i]` = Num at eta-power nmin+i-1.
# We add the transcendental series ops (exp, log1p, power-by-series, pochhammer,
# gamma-ratio) that the Coulomb construction needs but slice 1 doesn't provide.
# ---------------------------------------------------------------------------

# the imaginary-unit symbol shared with slice 2 (substitute => im to evaluate)
_IU() = imagunit_sym()

# zero / constant eta-series helpers (den=1)
_zero_eta(nmax::Int) = PNSeries(:η, Num[], nmax, nmax, 1)
function _const_eta(c, nmax::Int)
    cN = Num(c)
    iszero(cN) && return _zero_eta(nmax)
    PNSeries(:η, Num[cN], 0, nmax, 1)
end
# monomial c * eta^p
function _mono_eta(c, p::Int, nmax::Int)
    cN = Num(c)
    (iszero(cN) || p >= nmax) && return _zero_eta(nmax)
    PNSeries(:η, Num[cN], p, nmax, 1)
end

# truncate a PNSeries to O[eta]^(nmaxnew) (drop higher slots), keeping nmin.
# If the whole series sits at/above nmaxnew (nmin >= nmaxnew), collapse to a pure
# O[eta]^nmaxnew (empty coeffs at nmin=nmax=nmaxnew).
function _trunc_eta(s::PNSeries, nmaxnew::Int)
    nmaxnew >= s.nmax && return s
    s.nmin >= nmaxnew && return _zero_eta(nmaxnew)
    keep = max(nmaxnew - s.nmin, 0)
    PNSeries(:η, s.x0, s.coeffs[1:min(keep, length(s.coeffs))], s.nmin, nmaxnew, s.den)
end

# expand each coefficient to a flat sum-of-products (prevents deep unexpanded trees from
# `inv`/`*` blowing up downstream symbolic work). Applied at hot points, not everywhere.
_expand_coeffs(s::PNSeries) =
    PNSeries(s.var, s.x0, [Symbolics.expand(c) for c in s.coeffs], s.nmin, s.nmax, s.den)

# leading eta-power (nmax if all-zero)
function _lead_eta(s::PNSeries)
    i = findfirst(c -> !iszero(c), s.coeffs)
    i === nothing ? s.nmax : s.nmin + i - 1
end

# coeff of eta^p
function _coeff_eta(s::PNSeries, p::Int)
    (p < s.nmin || p - s.nmin + 1 > length(s.coeffs)) && return Num(0)
    s.coeffs[p - s.nmin + 1]
end

# exp of a series with NO constant term (lead order >= 1): exp(u) = sum u^k/k!,
# truncated at O[eta]^nmax. u must have _lead_eta(u) >= 1.
function _exp_series(u::PNSeries, nmax::Int)
    lead = _lead_eta(u)
    lead >= 1 || throw(ArgumentError("_exp_series needs a series with no constant term (lead>=1); got lead=$lead"))
    res = _const_eta(1, nmax)         # 1 + u + u^2/2 + ...
    term = _const_eta(1, nmax)
    k = 1
    while k * lead < nmax
        term = _trunc_eta(term * u, nmax) * (1 // k)
        res = res + term
        k += 1
    end
    res
end

# log(1 + u) for a series u with lead>=1: sum_{k>=1} (-1)^(k+1) u^k/k.
function _log1p_series(u::PNSeries, nmax::Int)
    lead = _lead_eta(u)
    lead >= 1 || throw(ArgumentError("_log1p_series needs lead>=1; got $lead"))
    res = _zero_eta(nmax)
    term = _const_eta(1, nmax)        # u^k accumulator start (will be *u)
    k = 1
    while k * lead < nmax
        term = _trunc_eta(term * u, nmax)
        res = res + term * (Num((-1)^(k + 1)) // k)
        k += 1
    end
    res
end

# base^expo as a series, where base = b0*(1+w) (b0 the Num leading coeff at some eta
# power p0) and expo is EITHER a Num (constant) or a PNSeries. Returns a series.
#   base^expo = b0^expo * exp(expo * log(1+w)) * eta^(p0*expo-ish)
# We require base's leading eta-power p0 to satisfy p0*expo integer-graded; handled by
# the callers which factor out the pure-power part explicitly. This helper expects
# `base` already normalized to lead 0 with leading coeff `b0` (a Num), i.e. base = b0(1+w).
function _powseries_normalized(b0::Num, w::PNSeries, expo, nmax::Int)
    # b0^expo
    pref = _b0_pow(b0, expo)
    inner = _log1p_series(w, nmax)            # log(1+w), lead>=1
    expinner = expo isa PNSeries ? _trunc_eta(expo * inner, nmax) : inner * Num(expo)
    eseries = _exp_series(expinner, nmax)
    _scale_series(eseries, pref)
end

# b0^expo where expo may be Num or PNSeries(eta). If expo is a PNSeries we split
# expo = e0 + edelta (e0 the eta^0 Num part, edelta lead>=1) and use
# b0^expo = b0^e0 * exp(edelta * log b0).  b0^e0 with e0 a Num: literal power.
function _b0_pow(b0::Num, expo)
    if expo isa PNSeries
        e0 = _coeff_eta(expo, 0)
        return b0^e0     # the edelta*log(b0) part is folded by the caller's exp
    else
        return b0^Num(expo)
    end
end

# scale every coeff of a series by a Num
_scale_series(s::PNSeries, c) = PNSeries(:η, s.x0, [Num(c) * cc for cc in s.coeffs], s.nmin, s.nmax, s.den)

# ---------------------------------------------------------------------------
# (small) building blocks for the Coulomb construction (a==0 specialization)
# ---------------------------------------------------------------------------

# Pochhammer (x)_n = x(x+1)...(x+n-1) for integer n>=0 (n<0 via reciprocal of the
# shifted product). `x` is a PNSeries; returns a PNSeries. n is a (small) integer.
function _pochhammer_series(x::PNSeries, n::Int, nmax::Int)
    if n == 0
        return _const_eta(1, nmax)
    elseif n > 0
        res = _const_eta(1, nmax)
        for k in 0:(n - 1)
            res = _trunc_eta(res * (x + Num(k)), nmax)
        end
        return res
    else
        # (x)_n = 1 / (x-1)(x-2)...(x-|n|)  ... actually (x)_n for n<0 =
        # 1/((x-1)(x-2)...(x+n)) = prod_{k=1}^{-n} 1/(x-k)
        res = _const_eta(1, nmax)
        for k in 1:(-n)
            res = _trunc_eta(res * inv(x - Num(k)), nmax)
        end
        return res
    end
end

# Gamma[k + x]/Gamma[k] as an eta-series, k a POSITIVE integer, x a series with lead>=1.
#   = exp( sum_{j>=1} c_j x^j ),  c_1 = psi(k),  c_j = (-1)^j Zeta(j,k)/j (j>=2),
# where psi(k) = -EulerGamma + H_{k-1}, Zeta(j,k) = Hurwitz zeta = Zeta(j) - sum_{i=1}^{k-1} i^-j.
# Constants kept EXACT as Num (EulerGamma, pi via Zeta(2)=pi^2/6, Zeta(j) for j>=3).
function _gamma_ratio_series(k::Int, x::PNSeries, nmax::Int)
    k >= 1 || throw(ArgumentError("_gamma_ratio_series: k must be a positive integer; got $k"))
    lead = _lead_eta(x)
    lead >= 1 || throw(ArgumentError("_gamma_ratio_series needs x with lead>=1"))
    # accumulate S = sum_j c_j x^j (a series with lead>=1), then exp(S).
    S = _zero_eta(nmax)
    xpow = _const_eta(1, nmax)
    j = 1
    while j * lead < nmax
        xpow = _trunc_eta(xpow * x, nmax)
        cj = j == 1 ? _digamma_int(k) : Num((-1)^j // j) * _hurwitz_zeta(j, k)
        S = S + _scale_series(xpow, cj)
        j += 1
    end
    _exp_series(S, nmax)
end

# Gamma[K + x] for ANY integer K and a series x (lead>=1), as an eta-Laurent series.
# Uses Gamma[K+x] = Gamma[1+x] * (1+x)_{K-1} for K>=1, and the pole forms for K<=0:
#   Gamma[x]     = Gamma[1+x]/x                              (K=0)
#   Gamma[K+x]   = Gamma[1+x]/(x (x-1)...(x+K))   (K<=-1)
# where Gamma[1+x] = _gamma_ratio_series(1,x) (= exp(sum c_j x^j), pole-free).
# The result may be Laurent (a 1/x^p pole) when K<=0 -- this captures the leading-order
# enhancement of the R_C^{-nu-1} Coulomb solution (the dominant Up terms).
function _gamma_full_series(K::Int, x::PNSeries, nmax::Int)
    g1x = _gamma_ratio_series(1, x, nmax)        # Gamma[1+x]/Gamma[1] = Gamma[1+x]
    if K >= 1
        # * prod_{i=1}^{K-1} (i + x)
        acc = g1x
        for i in 1:(K - 1)
            acc = _trunc_eta(acc * (x + Num(i)), nmax)
        end
        return acc
    else
        # divide by x * (x-1) * ... * (x+K)   [i.e. factors (x - t) for t=0..-K]
        denom = _const_eta(1, nmax)
        for t in 0:(-K)
            denom = _trunc_eta(denom * (x - Num(t)), nmax)   # (x - t), lead from x at t=0 -> Laurent
        end
        return _trunc_eta(g1x * inv(denom), nmax)
    end
end

# psi(k) = -EulerGamma + sum_{i=1}^{k-1} 1/i  (digamma at positive integer k)
function _digamma_int(k::Int)
    acc = -_EULERGAMMA
    for i in 1:(k - 1)
        acc += Num(1 // i)
    end
    acc
end

# Hurwitz zeta Zeta(j,k) = sum_{n>=0} 1/(n+k)^j = Zeta(j) - sum_{i=1}^{k-1} 1/i^j.
# Zeta(2) = pi^2/6 (kept symbolic via _PISYM); Zeta(j>=3) kept as a Num constant.
function _hurwitz_zeta(j::Int, k::Int)
    zj = j == 2 ? Num(1 // 6) * _PISYM^2 : _zeta_const(j)
    acc = zj
    for i in 1:(k - 1)
        acc -= Num(1 // (big(i)^j))
    end
    acc
end

# symbolic constants matching MMA's InputForm: EulerGamma, Pi, Zeta[n].
const _EULERGAMMA = Symbolics.variable(:EulerGamma)
const _PISYM = Symbolics.variable(:Pi)
_zeta_const(n::Int) = Symbolics.variable(Symbol("Zeta_$n"))

# ---------------------------------------------------------------------------
# eps-series (slice-2 ES{T}, numeric Complex coeffs) -> eta-graded PNSeries with
# Num coeffs.  eps^k -> coeff_k * 2^k  placed at eta-slot 3k (eps ~ eta^3).
# PERFORMANCE: omega is NOT carried symbolically through the build (it would bloat
# every product). Each eps^k contributes a NUMERIC coefficient at eta-slot 3k; the
# physical omega-power is fixed by the eta-grade and is reattached at the very end
# (`_reattach_omega`): a monomial r^p at eta-slot g carries omega^((g+2p)/3).
# The imaginary unit maps to the `imagunit_sym()` symbol.
# ---------------------------------------------------------------------------
function _es_to_etaseries(es, nmax::Int; tol = 0)
    lead = MST._es_lead(es; tol = tol)
    (lead < 0 || 3 * lead >= nmax) && return _zero_eta(nmax)
    coeffs = Num[]
    nmin = 3 * lead
    for k in lead:es.ord
        g = 3 * k
        g >= nmax && break
        while nmin + length(coeffs) < g
            push!(coeffs, Num(0))
        end
        c = MST._es_coeff(es, k)
        cval = (tol > 0 && abs(c) <= tol) ? Num(0) : _cnum(c) * (Num(2)^k)
        push!(coeffs, cval)
    end
    isempty(coeffs) && return _zero_eta(nmax)
    PNSeries(:η, coeffs, nmin, nmax, 1)
end

# Reattach omega to a built (omega-free) radial PNSeries. For each eta-slot g, the
# coefficient is expanded in r; the r^p term gets multiplied by omega^((g+2p)/3).
# Also handles 1/r^p (negative r-power) terms (Up). `ωsym` the omega Num.
function _reattach_omega(s::PNSeries, ωsym::Num, rsym::Num)
    newcoeffs = Vector{Num}(undef, length(s.coeffs))
    for (i, c) in enumerate(s.coeffs)
        g = s.nmin + i - 1
        newcoeffs[i] = _attach_omega_coeff(c, g, ωsym, rsym)
    end
    PNSeries(:η, s.x0, newcoeffs, s.nmin, s.nmax, s.den)
end

# multiply each r^p monomial of coefficient `c` (a polynomial or Laurent poly in r) by
# omega^((g+2p)/3).  Symbolics `coeff` cannot read powers off a fraction (1/r terms), so
# we CLEAR denominators: pick P = max negative r-power, form the polynomial pe = c*r^P
# (all r-powers now >= 0), read coeff of r^k => physical r-power (k-P), reattach omega.
function _attach_omega_coeff(c, g::Int, ωsym::Num, rsym::Num)
    iszero(c) && return Num(0)
    # combine into a single fraction num(r)/r^P, then read powers off the polynomial num.
    cf = Symbolics.simplify_fractions(Num(c))
    P = _max_inv_r_power(cf, rsym)                    # clear 1/r^P
    pe = Symbolics.expand(Symbolics.simplify_fractions(cf * rsym^P))  # polynomial in r
    kmax = Series._degree_in(pe, rsym)
    acc = Num(0)
    for k in 0:kmax
        ck = Series._coeff_at(pe, rsym, k)
        iszero(ck) && continue
        p = k - P                                    # physical r-power
        q = g + 2 * p
        (q % 3 == 0) || error("non-integer omega power: g=$g p=$p (q3=$q)")
        acc += ck * rsym^p * ωsym^(div(q, 3))
    end
    acc
end

# Largest P such that c has a 1/r^P term (0 if c is already a polynomial in r). Found by
# probing: multiply by r^P and test whether the result is fraction-free in r (i.e. the
# r->1/u substitution leaves no u in any denominator). Bounded probe.
function _max_inv_r_power(c::Num, rsym::Num)
    for P in 0:64
        cleared = Symbolics.expand(Symbolics.simplify_fractions(c * rsym^P))
        _is_r_polynomial(cleared, rsym) && return P
    end
    error("_max_inv_r_power: r-denominator deeper than 64; coefficient = $c")
end

# true if `e` is a polynomial in r (no r in any denominator). Detected by a NUMERIC finiteness
# probe: set every other symbolic variable to a generic value (the imaginary unit -> im) and r -> 0;
# a 1/r^k denominator blows up to Inf, a genuine r-polynomial stays finite. (`Symbolics.degree`
# under Symbolics 7 silently returns 0 on a `1/r` fraction instead of throwing, so the old degree
# probe mis-classified Laurent terms as polynomials and dropped the inverse-r power.)
function _is_r_polynomial(e, rsym::Num)
    en = Num(e)
    iu = _IU()
    sub = Dict{Num,Any}()
    for v in Symbolics.get_variables(Symbolics.value(en))
        vn = Num(v)
        sub[vn] = isequal(vn, rsym) ? 0.0 : (isequal(vn, iu) ? im : 2.0)
    end
    val = try
        Symbolics.value(Symbolics.substitute(en, sub))
    catch
        return false
    end
    return (val isa Number) && isfinite(ComplexF64(val))
end

# Complex numeric -> Num using the imagunit_sym symbol for the imaginary part.
function _cnum(c)
    re = real(c); ie = imag(c)
    Num(re) + Num(ie) * _IU()
end

# ---------------------------------------------------------------------------
# The Coulomb-wave radial solution R_C^{nu} (and, with `flip`, R_C^{-nu-1}).
# Spin shifts (PN.wl c["In"] l.2333, cD l.2342/2352, z[a,r] l.2330):
#   kappa = Sqrt[1-a^2]; eps = 2 omega; tau = (eps - m a)/kappa; epsp = (eps+tau)/2;
#   z = (kappa-1+r) omega;
#   R_C = c_In(z) * sum_{n,j} cD[n,j] z^(n+j)            (PN.wl 2333-2425)
#   c_In = E^(-I z) 2^nu z^nu (z-eps kappa)^(-s) (1-eps kappa/z)^(-I epsp)
#   cD["nu"][n,j] = (-1)^n (2I)^(n+j) Gamma[n+nu+1-s+I eps] PH[nu+1+s-I eps,n]
#                   PH[n+nu+1-s+I eps,j] aMST[n]
#                   / (Gamma[2n+2nu+2] PH[nu+1-s+I eps,n] PH[2n+2nu+2,j] j!)
#   cD["-nu-1"][n,j] = (-1)^n (2I)^(n+j) PH[s-I eps-nu,n] PH[n-s+I eps-nu,j]
#                      Gamma[n-s+I eps-nu] aMST[-n]
#                      / (j! PH[2n-2nu,j] PH[-s+I eps-nu,n] Gamma[2n-2nu])
# The cD coefficients depend on spin ONLY through nu(=nuMST) and aMST (the a-dependent
# slice-2 series); they carry NO bare kappa/tau/epsp. Only c_In and z do.
# At a==0: kappa=1, tau=eps, epsp=eps, z = r omega, (z-eps), (1-eps/z) — recovered
# byte-for-byte (kappa-1=0; epspser collapses to eps; epskappa=eps).
# `nu` here is the renormalized angular-momentum eta-series; for the -nu-1 solution
# the caller passes nu -> -nu-1 (a series) and uses cD["-nu-1"].
# Args: s integer; l (leading z-power / nMin bound); aMSTser n->etaseries; nuser the
# nu eta-series; epsser = eps eta-series (2 omega); kap = kappa (a Num PN-constant);
# epspser = epsp eta-series; ωsym; nmax eta-cutoff.
# `flip=false` => R_C^{nu}; `flip=true` => R_C^{-nu-1}.
# ---------------------------------------------------------------------------
# z^k for integer k, z a (multi-term) eta-series leading at eta^1 (coeff r). Repeated
# multiply (k>=0) or inv of that (k<0). At a==0, z is the single monomial r so this
# returns r^k exactly (the old `z^nu0` path).
function _z_pow_int(z::PNSeries, rsym::Num, k::Int, nmax::Int)
    k == 0 && return _const_eta(1, nmax)
    acc = _const_eta(1, nmax)
    for _ in 1:abs(k)
        acc = _trunc_eta(acc * z, nmax)
    end
    k > 0 ? acc : inv(acc)
end

function _coulomb_R(s::Int, l::Int, aMSTser::Dict{Int,PNSeries}, nuser::PNSeries,
        epsser::PNSeries, kap::Num, epspser::PNSeries, rsym::Num, nmax::Int; flip::Bool)

    iu = _IU()
    eps = epsser                                   # 2 (omega reattached later), eta^3
    # z = (kappa-1+r) omega: r*omega at eta^1 + (kappa-1)*omega at eta^3 (a const at slot 3).
    z = _mono_eta(rsym, 1, nmax) + _mono_eta(kap - Num(1), 3, nmax)
    logz = Num(log(rsym))                          # grade-0 marker (log z; omega-free,
                                                   # exact only where logs don't yet appear)

    # nu used inside the kernel: nu (flip=false) or -nu-1 (flip=true)
    nu = flip ? (_const_eta(-1, nmax) - nuser) : nuser
    nu0 = Int(round(Float64(Symbolics.value(_coeff_eta(nu, 0)))))   # integer leading nu (=l or -l-1)
    delta = nu - _const_eta(nu0, nmax)             # nu - nu0, lead>=? (lead 6 for nuser)

    # eps*kappa = 2 kappa at eta^3 (omega reattached later). At a==0 (kappa=1): 2.
    epskap = _scale_series(eps, kap)               # 2 kappa, eta^3

    # --- prefactor c_In(z) ---
    # E^(-I z)
    pf_exp = _exp_series(_scale_series(z, -iu), nmax)
    # 2^nu = 2^nu0 * exp(delta * log 2).  Use a rational power so nu0<0 (the -nu-1
    # solution has nu0 = -l-1) does not hit Julia's integer-pow DomainError.
    pf_2 = _scale_series(_exp_series(_scale_series(delta, Num(log(2))), nmax),
                         Num(Rational{BigInt}(2)^nu0))
    # z^nu = z^nu0 * exp(delta * log z).  z is no longer a monomial for a!=0, so z^nu0 is a
    # series (still leading eta^nu0): power-by-series via the leading r*omega term.
    znu0 = _z_pow_int(z, rsym, nu0, nmax)
    pf_z = _trunc_eta(znu0 * _exp_series(_scale_series(delta, logz), nmax), nmax)
    # (z - eps kappa)^(-s):  z - eps kappa = r*omega + (kappa-1)*omega - 2 kappa*omega
    #   = r*omega + (-kappa-1)*omega.  r*omega at eta^1; the const (-kappa-1) at eta^3.
    zme = _mono_eta(rsym, 1, nmax) + _mono_eta(-kap - Num(1), 3, nmax)
    pf_zme = (-s) >= 0 ? zme^(-s) : inv(zme^(s))
    # (1 - eps kappa/z)^(-I epsp).  -eps kappa/z is a series (lead eta^2 with coeff
    # -2 kappa/r at a==0; the z-tail makes it multi-term for a!=0).
    w = _scale_series(_trunc_eta(epskap * inv(z), nmax), -Num(1))   # -eps kappa / z
    expo = _scale_series(epspser, -iu)             # -I epsp
    pf_pow = _powseries_normalized(Num(1), w, expo, nmax)   # (1+w)^expo, b0=1

    cIn = _expand_coeffs(_trunc_eta(_trunc_eta(_trunc_eta(pf_exp * pf_2, nmax) * pf_z, nmax) *
                     _trunc_eta(pf_zme * pf_pow, nmax), nmax))

    # --- double sum over (n, j) ---
    # leading eta-order of R_C: -s+l (In) at flip=false. We include all (n,j) whose
    # term leading eta-order < nmax. z^(n+j) ~ eta^(n+j). The prefactor adds eta^(l-s)
    # (z^nu0 ~ eta^nu0, (z-eps)^-s ~ eta^-s). So term grade ~ nu0 - s + (n+j) + grade(cD-eps-part).
    # cIn contributes leading eta-order (l - s) (= nu0 - s, with nu0=l for flip=false,
    # nu0=-l-1 for flip=true). The double sum is multiplied by cIn, so we only need
    # `sumser` to eta-order  nmax - leadcIn. Each term's eta-order: aMST[n] leads at
    # 3|n| (eps^|n|-suppressed), z^(n+j) adds (n+j). So contributions decay quickly.
    leadcIn = nu0 - s
    sumcut = nmax - leadcIn               # sumser only needed to this eta-order
    sumcut <= 0 && return _zero_eta(nmax)
    sumser = _zero_eta(nmax)
    # n-range: positive n suppressed ~4n; negative n ~2|n|. Bound generously but finite.
    nlo = -(sumcut + l + 3)
    nhi = sumcut + 2
    for n in nlo:nhi
        an = flip ? get(aMSTser, -n, nothing) : get(aMSTser, n, nothing)
        an === nothing && continue
        _lead_eta(an) >= sumcut && continue
        cDn = _cD_term(s, l, n, nu, eps, an, iu, sumcut, flip)
        cDn === nothing && continue
        (Pn, jbuilder) = cDn
        _lead_eta(Pn) >= sumcut && continue        # Pn already past cutoff (incl. z^n via Pn? no)
        zn = _z_pow_int(z, rsym, n, sumcut)        # z^n (multi-term-safe; truncated)
        pnzn = _trunc_eta(Pn * zn, sumcut)
        _lead_eta(pnzn) >= sumcut && continue
        # j sum: term_j = jbuilder(j) * z^j ; accumulate while leading order < sumcut.
        jacc = _zero_eta(sumcut)
        consecutive_skips = 0
        for j in 0:(sumcut + 2)
            tj = _trunc_eta(jbuilder(j) * _z_pow_int(z, rsym, j, sumcut), sumcut)
            if _lead_eta(tj) >= sumcut
                consecutive_skips += 1
                consecutive_skips > (2 * l + 2) && break
                continue
            end
            consecutive_skips = 0
            jacc = jacc + tj
        end
        jacc = _expand_coeffs(jacc)
        sumser = sumser + _expand_coeffs(_trunc_eta(pnzn * jacc, sumcut))
    end

    _expand_coeffs(_trunc_eta(cIn * _expand_coeffs(sumser), nmax))
end

# cD[n,j] split into (Pn, jbuilder): Pn = the n-dependent, j-independent factor
# (a PNSeries); jbuilder(j) = the j-dependent factor (a PNSeries), so that
# cD[n,j] = Pn * jbuilder(j).  Returns `nothing` if a Gamma pole makes the n-term
# absent (non-positive integer Gamma argument with no compensating numerator).
function _cD_term(s::Int, l::Int, n::Int, nu::PNSeries, eps::PNSeries,
        an::PNSeries, iu::Num, nmax::Int, flip::Bool)
    Ieps = _scale_series(eps, iu)                  # I eps, eta^3
    nu0 = Int(round(Float64(Symbolics.value(_coeff_eta(nu, 0)))))
    if !flip
        # numerator Gamma[n+nu+1-s+I eps]; denom Gamma[2n+2nu+2]
        K1 = n + nu0 + 1 - s; x1 = (nu - _const_eta(nu0, nmax)) + Ieps     # Gamma[K1+x1]
        K2 = 2n + 2 * nu0 + 2; x2 = _scale_series(nu - _const_eta(nu0, nmax), 2)  # Gamma[K2+x2]
        gnum = _gamma_full_series(K1, x1, nmax)    # Gamma[K1+x1] (may be Laurent if K1<=0)
        gden = _gamma_full_series(K2, x2, nmax)    # Gamma[K2+x2]
        gfac = Num(1)
        # PH[nu+1+s-I eps, n] / PH[nu+1-s+I eps, n]
        phA = _poch_seriesfull(nu + _const_eta(1 + s, nmax) - Ieps, n, nmax)
        phB = _poch_seriesfull(nu + _const_eta(1 - s, nmax) + Ieps, n, nmax)
        Pn = _trunc_eta(_trunc_eta(_scale_series(gnum, gfac) * inv(gden), nmax) *
                        _trunc_eta(phA * inv(phB), nmax), nmax)
        Pn = _expand_coeffs(_scale_series(_trunc_eta(Pn * an, nmax), Num((-1)^n)))
        # jbuilder: (2I)^(n+j) PH[n+nu+1-s+I eps, j] / (PH[2n+2nu+2, j] j!)
        baseA = nu + _const_eta(n + 1 - s, nmax) + Ieps   # for PH[.,j]
        baseB = nu + _const_eta(2n + 2, nmax)              # nu*2? NO: 2n+2nu+2 -> see below
        baseB = (_scale_series(nu, 2) + _const_eta(2n + 2, nmax))
        twoI = Num(2) * iu
        jbuilder = function (j)
            num = _poch_seriesfull(baseA, j, nmax)
            den = _poch_seriesfull(baseB, j, nmax)
            _scale_series(_trunc_eta(num * inv(den), nmax), (twoI)^(n) * (twoI)^(j) / Num(factorial(big(j))))
        end
        return (Pn, jbuilder)
    else
        # flip: nu is ALREADY -nu-1 here (caller passed). Use cD["-nu-1"] verbatim:
        #   num: PH[s-I eps-nuOrig, n] = PH[s-I eps+(nu+1), n]? careful. The source uses
        #   nuMST (original). Here `nu` = -nuOrig-1, so nuOrig = -nu-1.
        nuO = _const_eta(-1, nmax) - nu             # original nuMST series
        nuO0 = Int(round(Float64(Symbolics.value(_coeff_eta(nuO, 0)))))
        # Gamma[n - s + I eps - nuOrig]; Gamma[2n - 2 nuOrig]
        K1 = n - s - nuO0; x1 = Ieps - (nuO - _const_eta(nuO0, nmax))
        K2 = 2n - 2 * nuO0; x2 = _scale_series(-(nuO - _const_eta(nuO0, nmax)), 2)
        gnum = _gamma_full_series(K1, x1, nmax)    # may be Laurent (K1<=0): the Up enhancement
        gden = _gamma_full_series(K2, x2, nmax)
        gfac = Num(1)
        # PH[s - I eps - nuOrig, n] / PH[-s + I eps - nuOrig, n]
        phA = _poch_seriesfull(_const_eta(s, nmax) - Ieps - nuO, n, nmax)
        phB = _poch_seriesfull(_const_eta(-s, nmax) + Ieps - nuO, n, nmax)
        Pn = _trunc_eta(_trunc_eta(_scale_series(gnum, gfac) * inv(gden), nmax) *
                        _trunc_eta(phA * inv(phB), nmax), nmax)
        Pn = _expand_coeffs(_scale_series(_trunc_eta(Pn * an, nmax), Num((-1)^n)))
        baseA = _const_eta(n - s, nmax) + Ieps - nuO        # for PH[n-s+I eps-nuOrig, j]
        baseB = _const_eta(2n, nmax) - _scale_series(nuO, 2) # PH[2n-2nuOrig, j]
        twoI = Num(2) * iu
        jbuilder = function (j)
            num = _poch_seriesfull(baseA, j, nmax)
            den = _poch_seriesfull(baseB, j, nmax)
            _scale_series(_trunc_eta(num * inv(den), nmax), (twoI)^(n) * (twoI)^(j) / Num(factorial(big(j))))
        end
        return (Pn, jbuilder)
    end
end

# Pochhammer (base)_n where `base` is a PNSeries (full, may have integer + series part),
# integer n (>=0 typical; n<0 handled). Returns PNSeries.
function _poch_seriesfull(base::PNSeries, n::Int, nmax::Int)
    _pochhammer_series(base, n, nmax)
end

# ---------------------------------------------------------------------------
# In / Up assembly (a==0, Default normalization). PN.wl RadialAssociationBoth 2893-2930.
#   R^In = R_C^{nu}            (+ K R_C^{-nu-1}, absent for order < gap = 4l+2)
#   R^Up = C1 R_C^{-nu-1} + C2 R_C^{nu}
#   coeffUp = -I E^(-pi eps - I pi s) Sin[pi(nu+s-I eps)] / Sin[2 pi nu]
#   C1 = coeffUp;  C2 = coeffUp * I E^(-I pi nu) Sin[pi(nu-s+I eps)]/Sin[pi(nu+s-I eps)]
# all SeriesTake'd to `order` (eta-cutoff nmax = order + leadingOrder).
# ---------------------------------------------------------------------------

# build the raw eps/nu/aMST eta-series inputs + spin shifts (kappa, epsp) for (s,l,m,a,order).
# a==0: kappa=1 (exact), epsp = eps (single eta^3 term, value 2). a!=0: kappa=Sqrt[1-a^2]
# (ComplexF64-irrational via slice-2's _spin_type), tau=(eps-m a)/kappa, epsp=(eps+tau)/2
# carries a grade-0 constant -m a/(2 kappa) (eta^0) plus the grade-3 part (1 + 1/kappa)
# (the eps/2 + 1/(2 kappa) coefficient). nu/aMST already a-dependent (slice-2 auto-lambda).
function _radial_inputs(s::Int, l::Int, m::Int, a::Real, order::Int, nmaxIn::Int)
    # eps-order budget: each eps ~ eta^3, so eta^nmax needs eps up to ceil(nmax/3).
    # The Coulomb sum's prefactor/coefficients combine several eps-graded pieces but the
    # NET eta-grade of any eps^k contribution is 3k, so eps-order ceil(nmax/3)+2 suffices.
    epsOrd = cld(nmaxIn, 3) + 2
    # a_n keys needed: the z^(n+j) sum contributes at eta-grade ~ n (n>0) or via the
    # prefactor's eta^(l-s) plus n; only |n| up to ~ ceil(nmaxIn/3)+l matters (a_n itself
    # is eps^|n|-suppressed). Keep a modest range.
    nkeys = cld(nmaxIn, 3) + l + 3
    nu_es, aeps, kappa, av, ordInt =
        MST._mst_eps_series(s, l, m, a, epsOrd; nkeys = nkeys)
    tol = 0
    nuser = _es_to_etaseries(nu_es, nmaxIn; tol = tol)           # omega-free
    epsser = _mono_eta(Num(2), 3, nmaxIn)                        # eps = 2 (omega reattached later), eta^3
    aMSTser = Dict{Int,PNSeries}()
    for (n, es) in aeps
        aMSTser[n] = _es_to_etaseries(es, nmaxIn; tol = tol)
    end
    # kappa as a Num PN-constant; epsp = (eps + tau)/2 as a 2-slot eta-series.
    #   tau = (eps - m a)/kappa  ->  grade-0: -m a/kappa,  grade-3: 2/kappa (from eps=2).
    #   epsp = (eps + tau)/2     ->  grade-0: -m a/(2 kappa),  grade-3: (2 + 2/kappa)/2 = 1 + 1/kappa.
    kapnum = _cnum(kappa)
    invkap = _cnum(inv(kappa))
    epsp0 = _cnum(-m * av * inv(kappa) / 2)        # grade-0 constant (0 at a==0)
    epsp3 = Num(1) + invkap                         # grade-3 coeff (= 2 at a==0)
    epspser = (iszero(epsp0) ? _zero_eta(nmaxIn) : _const_eta(epsp0, nmaxIn)) +
              _mono_eta(epsp3, 3, nmaxIn)
    (nuser, aMSTser, epsser, kapnum, epspser)
end

# working eta-cutoff for a given boundary condition. Up multiplies factors with NEGATIVE
# leading eta-orders (C1,C2 ~ eta^(-2l+1) from 1/Sin[2 pi nu], and R_C^{-nu-1} ~ the Up
# lead eta^(-s-l-1)); the Cauchy product truncates at the SMALLER factor nmax, so every
# factor must reach the target nmax PLUS the depth of those negative leads. We add both the
# C-coefficient depth (~2l) and the Up-lead depth (|-s-l-1|) plus slack.
function _radial_work(s::Int, l::Int, order::Int, bc::Symbol)
    if bc === :Up
        lead = -s - l - 1
        nmax = lead + order
        # 4l+2 covers the C1/C2 (~eta^(-2l+1)) and R_C^{-nu-1} negative-lead depths for the
        # s=-2 dominant mode (Up lead eta^-1; oracle-verified to ~1e-16). The deeper Up lead
        # at s=0 (eta^-3) needs more headroom: each extra unit of -lead-1 costs ~3 orders of
        # the inner Coulomb series (verified: s=0 l=2 order-4 needs work>=15, oracle-exact).
        return nmax + 4 * l + 2 + 3 * max(0, -lead - 1)
    else
        return (-s + l) + order + 4
    end
end

# the In/Up SeriesData, as a PNSeries in :eta with Num coeffs in (r, omega). `a` is the
# (numeric) Kerr spin; a==0 recovers the Schwarzschild build byte-for-byte.
# `inputs` (a (nuser, aMSTser, epsser, kappa, epspser) tuple at >= the needed work) may be
# passed to share one slice-2 eps-series computation between the In and Up builds.
function _radial_solution(s::Int, l::Int, m::Int, order::Int, bc::Symbol;
        a::Real = 0, rsym::Num = Symbolics.variable(:r), inputs = nothing)
    s in (0, -1, -2) || throw(ArgumentError("radial PN port supports s in {0,-1,-2}; got s=$s"))
    l >= abs(s) || throw(ArgumentError("need l >= |s|; got l=$l, s=$s"))
    l >= abs(m) || throw(ArgumentError("need l >= |m|; got l=$l, m=$m"))

    leadIn = -s + l                       # In leading eta-order
    leadUp = -s - l - 1                   # Up leading eta-order
    lead = bc === :In ? leadIn : leadUp
    nmax = lead + order                   # O[eta]^(lead+order): keep `order` terms
    work = _radial_work(s, l, order, bc)

    nuser, aMSTser, epsser, kap, epspser =
        inputs === nothing ? _radial_inputs(s, l, m, a, order, work) : inputs

    if bc === :In
        # InGap (PN.wl 2558-2561): 4l+2 when m*a==0 (Schwarzschild, or m==0); else 4l-1.
        gap = iszero(m * a) ? 4 * l + 2 : 4 * l - 1
        Rnu = _coulomb_R(s, l, aMSTser, nuser, epsser, kap, epspser, rsym, work; flip = false)
        if order >= gap
            # second term K * R_C^{-nu-1}; K = KAmplitude["Ratio"]. DEFERRED for now.
            throw(ArgumentError("In second term (order >= gap=$gap) not yet ported (l=$l order=$order a=$a)"))
        end
        R = Rnu
    elseif bc === :Up
        Rnu  = _coulomb_R(s, l, aMSTser, nuser, epsser, kap, epspser, rsym, work; flip = false)
        Rmnu = _coulomb_R(s, l, aMSTser, nuser, epsser, kap, epspser, rsym, work; flip = true)
        iu = _IU()
        # C1, C2 as eta-series (functions of nu, eps). Sin[pi x] expanded as series.
        C1, C2 = _up_coeffs(s, nuser, epsser, iu, work)
        R = _trunc_eta(C1 * Rmnu, work) + _trunc_eta(C2 * Rnu, work)
    else
        throw(ArgumentError("bc must be :In or :Up; got $bc"))
    end

    # SeriesTake to `order` terms from the true leading order, then reattach omega
    # (the omega-power of each r^p monomial at eta-slot g is (g+2p)/3).
    R = _seriestake_from_lead(R, lead, order)
    ωsym = Symbolics.variable(:ω)
    R = _reattach_omega(R, ωsym, rsym)
    R
end

# Sin[pi*(integer K + x_series)] = (-1)^K * sin(pi*x_series) expanded for series x
# (lead>=1): sin(pi x) = pi x - (pi x)^3/6 + ... ; the leading (-1)^K pi x sets the order.
# Returns a PNSeries.
function _sin_pi_series(K::Int, x::PNSeries, nmax::Int)
    # sin(pi x) for x lead>=1
    lead = _lead_eta(x)
    lead >= 1 || throw(ArgumentError("_sin_pi_series needs x lead>=1; got $lead"))
    px = _scale_series(x, _PISYM)
    res = _zero_eta(nmax)
    term = px
    k = 0
    # sum_{k>=0} (-1)^k px^(2k+1)/(2k+1)!
    while (2k + 1) * lead < nmax
        p = 2k + 1
        term = px^p
        res = res + _scale_series(term, Num((-1)^k // factorial(big(p))))
        k += 1
    end
    _scale_series(res, Num((-1)^K))
end

# C1, C2 for the Up combination (PN.wl 2636-2638), a==0.
#   coeffUp = -I E^(-pi eps - I pi s) Sin[pi(nu+s-I eps)] / Sin[2 pi nu]
#   nu = l + delta (delta lead 6). nu + s - I eps = (l+s) + (delta - I eps).
#   2 nu = 2l + 2 delta.  Sin[2 pi nu] = (-1)^(2l) sin(2 pi delta) ~ lead 6.
function _up_coeffs(s::Int, nuser::PNSeries, epsser::PNSeries, iu::Num, nmax::Int)
    nu0 = Int(round(Float64(Symbolics.value(_coeff_eta(nuser, 0)))))   # = l
    delta = nuser - _const_eta(nu0, nmax)
    Ieps = _scale_series(epsser, iu)
    # Sin[pi(nu+s-I eps)] = Sin[pi((nu0+s) + (delta - I eps))]
    sinA = _sin_pi_series(nu0 + s, delta - Ieps, nmax)
    # Sin[2 pi nu] = Sin[pi(2 nu0 + 2 delta)] = (-1)^(2 nu0) sin(2 pi delta)
    sin2 = _sin_pi_series(2 * nu0, _scale_series(delta, 2), nmax)
    # E^(-pi eps): eps lead 3 -> series
    Eeps = _exp_series(_scale_series(epsser, -_PISYM), nmax)
    # E^(-I pi s) = (-1)^s for integer s
    Eis = Num((-1)^s)
    coeffUp = _scale_series(_trunc_eta(_trunc_eta(Eeps * sinA, nmax) * inv(sin2), nmax),
                            -iu * Eis)
    C1 = coeffUp
    # C2 = coeffUp * I E^(-I pi nu) Sin[pi(nu-s+I eps)]/Sin[pi(nu+s-I eps)]
    # E^(-I pi nu) = E^(-I pi nu0) * exp(-I pi delta)
    Eipinu = _scale_series(_exp_series(_scale_series(delta, -iu * _PISYM), nmax),
                           Num((-1)^nu0))
    sinC = _sin_pi_series(nu0 - s, delta + Ieps, nmax)   # Sin[pi(nu-s+I eps)]
    fac = _trunc_eta(_trunc_eta(Eipinu * sinC, nmax) * inv(sinA), nmax)
    C2 = _scale_series(_trunc_eta(coeffUp * fac, nmax), iu)
    (C1, C2)
end

# SeriesTake-from-leading: keep `order` terms starting at eta^lead, i.e. O[eta]^(lead+order).
function _seriestake_from_lead(s::PNSeries, lead::Int, order::Int)
    newmax = lead + order
    # ensure the stored slots start at `lead` (the true leading order); the build may
    # have produced exact zeros below lead which we drop.
    actual_lead = _lead_eta(s)
    actual_lead < lead && (lead = actual_lead)    # safety: never cut a real term
    nslots = newmax - lead
    coeffs = Num[_coeff_eta(s, lead + i) for i in 0:(nslots - 1)]
    PNSeries(:η, s.x0, coeffs, lead, newmax, 1)
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    TeukolskyRadialPN(s, l, m, a, ω, (η, order); Normalization="Default") -> Dict

PN-expanded homogeneous radial Teukolsky solutions (Teukolsky`PN`TeukolskyRadialPN,
PN.wl 2977). Returns `Dict("In" => PNSeries, "Up" => PNSeries)` — each a series in
`η` (= `ηvar`) whose coefficients are `Num` in `r`, `ω` (and the imaginary-unit symbol
`imagunit_sym()`, EulerGamma, Pi, ...). Spec: ST LRR 6 (2003) 6 Eqs. 159/162/166.

Spin enters via κ=√(1−a²), τ=(ε−ma)/κ, εp=(ε+τ)/2, z=(κ−1+r)ω (PN.wl replsKerr l.913,
c["In"] l.2333, z[a,r] l.2330); a≠0 carries irrational κ as ComplexF64 coefficients.

SCOPE: Default normalization, In second term absent (order < InGap = 4l+2 if m·a==0
else 4l−1). `s in {0,-1,-2}`. Substitute `imagunit_sym() => im` to evaluate.
"""
function TeukolskyRadialPN(s::Integer, l::Integer, m::Integer, a::Real,
        ωvar::Symbol, porder::Tuple{Symbol,<:Integer}; Normalization::AbstractString = "Default")
    pvar, order = porder
    order >= 1 || throw(ArgumentError("order must be >= 1; got $order"))
    Normalization == "Default" || throw(ArgumentError("only Normalization=\"Default\" is ported; got $Normalization"))

    rsym = Symbolics.variable(:r)
    si, li, mi, oi = Int(s), Int(l), Int(m), Int(order)
    # compute the slice-2 eps-series inputs ONCE at the larger (Up) working order and
    # share them between the In and Up builds.
    workmax = max(_radial_work(si, li, oi, :In), _radial_work(si, li, oi, :Up))
    inputs = _radial_inputs(si, li, mi, a, oi, workmax)
    Rin = _radial_solution(si, li, mi, oi, :In; a = a, rsym = rsym, inputs = inputs)
    Rup = _radial_solution(si, li, mi, oi, :Up; a = a, rsym = rsym, inputs = inputs)

    # rename the expansion variable :η -> pvar if requested, and :ω -> ωvar
    Rin = _rename_radial(Rin, pvar)
    Rup = _rename_radial(Rup, pvar)
    if ωvar != :ω
        Rin = _sub_omega(Rin, ωvar)
        Rup = _sub_omega(Rup, ωvar)
    end
    Dict{String,PNSeries}("In" => Rin, "Up" => Rup)
end

_rename_radial(s::PNSeries, v::Symbol) = PNSeries(v, s.x0, copy(s.coeffs), s.nmin, s.nmax, s.den)
function _sub_omega(s::PNSeries, ωvar::Symbol)
    ω = Symbolics.variable(:ω); ωnew = Symbolics.variable(ωvar)
    PNSeries(s.var, s.x0, [Symbolics.substitute(c, Dict(ω => ωnew)) for c in s.coeffs], s.nmin, s.nmax, s.den)
end
# ---------------------------------------------------------------------------
# Amplitude eps→eta helpers (a==0). Amplitudes are functions of ε=2ωγ only.
# ---------------------------------------------------------------------------

# An eps-series (ES, ω-free numeric, e.g. aMST[n] or νMST) -> ω-FREE η-series: ε^k -> 2^k η^(3k).
# ω is carried by GRADE (an η^(3k) slot ↔ ε^k ↔ ω^k) and reattached at the very end by
# `_reattach_omega_amp`. Coeffs are FLOAT (ComplexF64-valued Num): the exact MST coeffs are
# Complex{Rational{BigInt}} with huge denominators, and chaining them through the amplitude
# Pochhammer/Gamma `inv` products is catastrophically slow; the oracle check is numeric
# (rtol ~1e-7), so float arithmetic in the amplitude layer is both correct-to-tolerance and
# orders of magnitude faster. (The radial layer keeps exact rationals; only the amplitudes
# convert.)  ω is carried by grade, so each slot stays a pure numeric constant.
function _amp_es_eta(es, nmax::Int)
    lead = MST._es_lead(es)
    (lead < 0 || lead > es.ord || 3 * lead >= nmax) && return _zero_eta(nmax)
    coeffs = Num[]
    nmin = 3 * lead
    for k in lead:es.ord
        ek = 3 * k
        ek >= nmax && break
        while nmin + length(coeffs) < ek
            push!(coeffs, Num(0))
        end
        push!(coeffs, _cnum_float(MST._es_coeff(es, k)) * (Num(2.0)^k))
    end
    isempty(coeffs) && return _zero_eta(nmax)
    PNSeries(:η, coeffs, nmin, nmax, 1)
end
const _amp_aMST_eta = _amp_es_eta

# Complex numeric -> Num via Float64 real/imag (the imag-unit symbol for the imaginary part).
_cnum_float(c) = Num(Float64(real(c))) + Num(Float64(imag(c))) * _IU()

# ε = 2 η³ (ω-free; ω carried by grade).
_eps_eta(nmax::Int) = _mono_eta(Num(2), 3, nmax)

# the symbol for Log[ω] (substitute => log(ωnum) to evaluate). Log[ε]=Log[2]+Log[ω] keeps
# the Log[ω] piece symbolic (ω-free build); the η^3-grade of ε is NOT an explicit log.
const _LOGW = Symbolics.variable(:Logω)
_logw() = _LOGW

# Reattach ω to a built ω-free amplitude η-series: slot-g coeff (g=nmin+i-1) -> ×ω^(g/3).
# Every η³ of an amplitude is one factor of ε=2ω, so the ω-power at η-slot g is exactly g/3
# (g divisible by 3 for an integer-graded amplitude; den must be 1). `Logω` stays symbolic.
function _reattach_omega_amp(s::PNSeries, ωsym::Num)
    s.den == 1 || throw(ArgumentError("_reattach_omega_amp needs den=1; got den=$(s.den)"))
    newcoeffs = Vector{Num}(undef, length(s.coeffs))
    for (i, c) in enumerate(s.coeffs)
        g = s.nmin + i - 1
        (g % 3 == 0) || (iszero(c) ? nothing : error("amplitude η-slot $g not ÷3 (ω-grade)"))
        newcoeffs[i] = iszero(c) ? Num(0) : c * ωsym^(div(g, 3))
    end
    PNSeries(:η, s.x0, newcoeffs, s.nmin, s.nmax, 1)
end

# 2^n as a Num for ANY integer n. Use plain Int/Rational{Int} (NOT BigInt) so Symbolics keeps
# fast machine-rational arithmetic — a BigInt-typed Num poisons every downstream product.
_pow2(n::Int) = n >= 0 ? Num(2^n) : Num(1 // (2^(-n)))

# (-1)^n / k!  as a FLOAT Num (k = |n_fac|).  Float (not BigInt rational) so it doesn't poison
# the float amplitude arithmetic. The factorials here are small (|n| ≤ ~6).
_signfac(sgn::Int, kfac::Int) = Num(Float64((-1)^sgn) / Float64(factorial(kfac)))

# i^k for integer k (= E^(iπ k/2)): cycles 1, i, -1, -i.
function _phase_pow_i(k::Int)
    iu = _IU()
    r = ((k % 4) + 4) % 4
    r == 0 ? Num(1) : r == 1 ? iu : r == 2 ? Num(-1) : -iu
end

# Build (nuser, eps, aMSTeta dict) for the amplitude sums at η-cutoff nmax (a==0 path; kept
# for the a==0 callers / tests). ALL ω-free (ω carried by grade); only aMST[n] within the
# η-window are kept (|n| with 3|n| < nmax) to bound the Pochhammer-product work.
_amp_inputs(s::Int, l::Int, m::Int, nmax::Int) =
    _amp_inputs(s, l, m, 0, nmax)[1:3]

# a-threaded variant: returns (nuser, eps, aMSTeta, kap, tauser, epspser). The spin shifts
# (PN.wl AAmplitude/KAmplitude `κ=Sqrt[1-a²]`, `τ=(ε-ma)/κ`, `εp=(ε+τ)/2`) are ω-free η-series:
#   ε  = 2 η³  (slot 3, ω^1)
#   τ  = (-m a/κ) η⁰ + (2/κ) η³     (grade-0 ω^0 piece + grade-3 ω^1 piece)
#   εp = (-m a/(2κ)) η⁰ + (1+1/κ) η³
# At a==0: κ=1, τ=ε (the η⁰ piece vanishes), εp=ε — recovering the a==0 amplitudes.
function _amp_inputs(s::Int, l::Int, m::Int, a::Real, nmax::Int)
    epsOrd = cld(nmax, 3) + 1
    nmaxn = cld(nmax, 3) + 1
    nu_es, aeps, kappa, av, _ = MST._mst_eps_series(s, l, m, a, epsOrd; nkeys = nmaxn)
    nuser = _amp_es_eta(nu_es, nmax)                   # ω-free η-series (slot 3k)
    eps = _eps_eta(nmax)
    ncap = cld(nmax, 3)
    aMSTeta = Dict{Int,PNSeries}()
    for (n, es) in aeps
        abs(n) > ncap && continue
        et = _amp_aMST_eta(es, nmax)
        _lead_eta(et) >= nmax && continue              # outside the η-window: drop entirely
        aMSTeta[n] = et
    end
    kap = _cnum(kappa)
    invkap = _cnum(inv(kappa))
    tau0 = _cnum(-m * av * inv(kappa))                 # τ grade-0 (0 at a==0)
    tau3 = _scale_series(eps, invkap)                  # τ grade-3 = (2/κ) η³  (ω^1)
    tauser = (iszero(tau0) ? _zero_eta(nmax) : _const_eta(tau0, nmax)) + tau3
    epsp0 = _cnum(-m * av * inv(kappa) / 2)            # εp grade-0
    epsp3 = Num(1) + invkap                            # εp grade-3 coeff (= 2 at a==0)
    epspser = (iszero(epsp0) ? _zero_eta(nmax) : _const_eta(epsp0, nmax)) +
              _mono_eta(epsp3, 3, nmax)
    (nuser, eps, aMSTeta, kap, tauser, epspser)
end

# Reduce powers of the imaginary-unit symbol 𝕚 in a Num: 𝕚^k -> 𝕚^(k mod 4) with the
# sign/factor (𝕚²=-1). Critical for the amplitude layer: 𝕚 is a FREE symbol, so without
# this every Pochhammer/Gamma product accumulates ever-higher 𝕚-powers that Symbolics never
# collapses, exploding the expression tree. We rewrite 𝕚^2->-1, 𝕚^3->-𝕚, 𝕚^4->1, ... by
# substituting 𝕚^2 => -1 enough times, then expand.
# Clear any `N / D` subterm whose denominator D is a PURE complex constant (a polynomial in 𝕚
# only) by rationalizing: N/D -> N·(1/D) with 1/D = _complex_inv(D) a canonical re'+im'𝕚.
# Such fractions arise at a≠0 from `inv(series)` whose leading coeff is purely imaginary (e.g.
# `inv(δ - iε)` ⇒ 1/(2𝕚)); at a==0 they cancel against the (exact) real prefactors, but the
# float-complex a≠0 prefactors leave them, and `_reduce_iu`'s `_degree_in` cannot read a
# fraction. This restores a polynomial-in-𝕚 form. Pure-real denominators are left untouched.
function _clear_iu_denoms(e::Num)
    v = Symbolics.value(e)
    Symbolics.iscall(v) || return e
    op = Symbolics.operation(v)
    args = Symbolics.arguments(v)
    if op === (/) && length(args) == 2
        num = _clear_iu_denoms(Num(args[1]))
        den = Num(args[2])
        if _is_iu_const(den) && imag(_cval(den)) != 0
            return num * _complex_inv(den)
        end
        return num / _clear_iu_denoms(den)
    elseif op === (+)
        return sum(_clear_iu_denoms(Num(a)) for a in args)
    elseif op === (*)
        return prod(_clear_iu_denoms(Num(a)) for a in args)
    elseif op === (^) && length(args) == 2
        return _clear_iu_denoms(Num(args[1]))^Num(args[2])
    else
        return e
    end
end

function _reduce_iu(c::Num)
    iu = _IU()
    e = Symbolics.expand(_clear_iu_denoms(Num(c)))
    deg = Series._degree_in(e, iu)
    deg <= 1 && return e
    # collect by 𝕚-power: coeff_k 𝕚^k -> coeff_k · (𝕚^(k mod 4) value).  𝕚^0=1,1=𝕚,2=-1,3=-𝕚.
    acc = Num(0)
    for k in 0:deg
        ck = Series._coeff_at(e, iu, k)
        iszero(ck) && continue
        r = k % 4
        acc += ck * (r == 0 ? Num(1) : r == 1 ? iu : r == 2 ? Num(-1) : -iu)
    end
    Symbolics.expand(acc)
end

# amplitude expand for SUM ACCUMULATION (the hot inner loop): plain flat expand only — cheap.
# (𝕚-power reduction is deferred to _final_amp at each amplitude's end; doing it per-term is
# far too slow because _degree_in/expand run on every partial sum.)
_expand_amp(s::PNSeries) = _expand_coeffs(s)

# Is `c` a PURE complex constant (a polynomial in the imag-unit 𝕚 only, no other symbols)?
_is_iu_const(c::Num) = isempty(setdiff(Symbolics.get_variables(Symbolics.value(c)), [Symbolics.value(_IU())]))

# inverse of a pure complex constant `b0 = re + im 𝕚`, returned canonically as re' + im' 𝕚
# (1/(re+im 𝕚) = (re - im 𝕚)/(re²+im²)). Keeps the amplitude layer free of symbolic 1/𝕚
# fractions (the generic series `inv` would otherwise form 1/(complex) which `_reduce_iu`
# cannot expand). Pure-real b0 ⇒ re'/0 collapses to the exact 1/re (no 𝕚).
function _complex_inv(b0::Num)
    z = _cval(b0)
    iz = inv(z)
    _cnum(iz)
end

# amplitude-layer series inverse robust to a COMPLEX-constant leading coefficient: factor the
# lead complex const out (invert it via _complex_inv), invert the normalized lead-1 series with
# the generic `inv` (now fraction-free), reattach. Only engages when the lead is a complex
# constant WITH a nonzero imaginary part (a≠0 path); a==0 leads are real rationals ⇒ plain inv
# (exact, byte-for-byte). The `1/(re+im𝕚)` symbolic fraction that `_reduce_iu` cannot expand is
# thereby never formed.
function _amp_inv(s::PNSeries)
    li = findfirst(c -> !iszero(c), s.coeffs)
    li === nothing && return inv(s)
    b0 = s.coeffs[li]
    (_is_iu_const(b0) && imag(_cval(b0)) != 0) || return inv(s)
    ib0 = _complex_inv(b0)
    # normalize so the lead is EXACTLY 1: scale by ib0, 𝕚-reduce, then OVERWRITE the lead with
    # exactly Num(1) (b0·ib0 = 1 by construction; float roundoff leaves a ~1e-16 𝕚 residue that
    # would make `inv` form 1/(1+ε𝕚) — a fraction `_reduce_iu` cannot expand).
    sc = [_reduce_iu(c) for c in _scale_series(s, ib0).coeffs]
    sc[li] = Num(1)
    snorm = PNSeries(s.var, s.x0, sc, s.nmin, s.nmax, s.den)
    _scale_series(inv(snorm), ib0)
end

# amplitude FINALIZE: flat expand + 𝕚-power reduction (𝕚^k -> 𝕚^(k mod 4)). Applied ONCE at
# the end of each amplitude so the returned coeffs are in canonical {1, 𝕚, -1, -𝕚}·real form.
_final_amp(s::PNSeries) =
    PNSeries(s.var, s.x0, [_reduce_iu(c) for c in s.coeffs], s.nmin, s.nmax, s.den)

# Pochhammer PH[base, n] over an η-series base (reuse slice-3 helper).
_PH(base::PNSeries, n::Int, nmax::Int) = _poch_seriesfull(base, n, nmax)

# ---------------------------------------------------------------------------
# A-amplitudes (PN.wl 1688-1737). a==0: κ=1, ε=2ωγ. nu0 = l.
#   A["+"] = coeff_+ Σ_n aMST[n]
#   A["-"] = coeff_- Σ_n (-1)^n PH[ν+1+s-iε,n]/PH[ν+1-s+iε,n] aMST[n]
# ---------------------------------------------------------------------------
# (all ω-free; ω reattached by the caller via _reattach_omega_amp.)
function _A_amp(sign::Symbol, s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int)
    iu = _IU()
    nu0 = l
    delta = nuser - _const_eta(nu0, nmax)              # νMST - l, η-lead 6
    Ieps = _scale_series(eps, iu)                       # i ε, η-lead 3
    Eeps = _exp_series(_scale_series(eps, -_PISYM * (1 // 2)), nmax)   # E^(-π ε/2)
    if sign === :plus
        ph0 = _phase_pow_i(nu0 + 1 + s)                 # E^(iπ(l+1+s)/2)
        Ephase = _scale_series(_exp_series(_scale_series(delta, iu * _PISYM * (1 // 2)), nmax), ph0)
        E2 = _scale_series(_exp_series(_scale_series(Ieps, -Num(log(2))), nmax), _pow2(-1 + s))  # 2^(-1+s-iε)
        gnum = _gamma_full_series(nu0 + 1 - s, delta + Ieps, nmax)   # Γ[ν+1-s+iε]
        gden = _gamma_full_series(nu0 + 1 + s, delta - Ieps, nmax)   # Γ[ν+1+s-iε]
        coeff = _trunc_eta(_trunc_eta(_trunc_eta(Eeps * Ephase, nmax) * E2, nmax) *
                           _trunc_eta(gnum * _amp_inv(gden), nmax), nmax)
        ssum = _zero_eta(nmax)
        for (n, an) in aMSTeta
            _lead_eta(an) >= nmax && continue
            ssum = ssum + an
        end
    else  # :minus
        ph0 = _phase_pow_i(-(nu0 + 1 + s))
        Ephase = _scale_series(_exp_series(_scale_series(delta, -iu * _PISYM * (1 // 2)), nmax), ph0)
        E2 = _scale_series(_exp_series(_scale_series(Ieps, Num(log(2))), nmax), _pow2(-1 - s))   # 2^(-1-s+iε)
        coeff = _trunc_eta(_trunc_eta(Eeps * Ephase, nmax) * E2, nmax)
        baseA = nuser + _const_eta(1 + s, nmax) - Ieps
        baseB = nuser + _const_eta(1 - s, nmax) + Ieps
        ssum = _zero_eta(nmax)
        for (n, an) in aMSTeta
            _lead_eta(an) >= nmax && continue
            ph = _trunc_eta(_PH(baseA, n, nmax) * _amp_inv(_PH(baseB, n, nmax)), nmax)
            ssum = ssum + _scale_series(_trunc_eta(ph * an, nmax), Num((-1)^n))
        end
    end
    _final_amp(_trunc_eta(coeff * ssum, nmax))
end

# ---------------------------------------------------------------------------
# K-amplitudes (PN.wl 2072-2242). a==0, ScriptR=0, κ=1, ε=2ωγ, τ=ε (so iτ=iε).
# Γ of νMST args reuse _gamma_full_series (handles Laurent/pole K≤0 ⇒ the pole-cancelling
# 1/η factors cancel in the ratios automatically, the MMA ExpandGamma role).
# ---------------------------------------------------------------------------

# Split a τ-like ω-free η-series into (grade-0 COMPLEX const, hi-part with lead>=1). At a==0
# τ=ε so the grade-0 const is 0 and hi=the full series; at a≠0 τ carries a grade-0 -ma/κ.
function _split_grade0(x::PNSeries)
    c = _lead_eta(x) == 0 && length(x.coeffs) >= 1 ? x.coeffs[1] : Num(0)
    hi = iszero(c) ? x : (x - _const_eta(c, x.nmax))
    (c, hi)
end

# K["ν"]  (PN.wl 2072-2100, r=0):
#   coeff = E^(iεκ) (εκ)^(s-ν) 2^(-ν) Γ[1-s-2iεp] Γ[2ν+2]/Γ[ν+1-s+iε] · Γ[1+2ν]/(Γ[1-s-iε+ν]Γ[1+ν-iτ])
#   sumUp   = Σ_n(-1)^n PH[1+s+iε+ν,n]PH[1+2ν,n]PH[1+ν+iτ,n]/(n!PH[1-s-iε+ν,n]PH[1+ν-iτ,n]) aMST[n]
#   sumDown = Σ_n(-1)^n/((-n)!PH[2ν+2,n]) PH[ν+1+s-iε,n]/PH[ν+1-s+iε,n] aMST[n]
# a==0: κ=1, τ=ε, εp=ε — the four spin args default to (1, ε, ε) byte-for-byte.
function _K_nu(s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int,
        kap::Num = Num(1), tauser::PNSeries = eps, epspser::PNSeries = eps)
    iu = _IU()
    nu0 = l
    delta = nuser - _const_eta(nu0, nmax)
    Ieps = _scale_series(eps, iu)
    Itau = _scale_series(tauser, iu)                     # i τ  (2-slot at a≠0)
    Iepsp = _scale_series(epspser, iu)                   # i εp
    tau_c, tau_hi = _split_grade0(Itau)                  # iτ = tau_c + tau_hi (hi lead>=1)
    # coeff
    # (εκ)^(s-ν) = (εκ)^(s-l) · exp(-δ log(εκ)); ω-free.
    eps_slm = _epskap_intpow(s - nu0, kap, nmax)                        # (εκ)^(s-l)
    logek = _log_epskap(kap, nmax)                                      # log(εκ)
    eps_delta = _exp_series(_scale_series(delta, -Num(1)) * logek, nmax) # (εκ)^(-δ)
    eps_pow = _trunc_eta(eps_slm * eps_delta, nmax)                     # (εκ)^(s-ν)
    Eieps = _exp_series(_scale_series(Ieps, kap), nmax)                 # E^(iεκ)
    two_mnu = _scale_series(_exp_series(_scale_series(delta, -Num(log(2))), nmax), _pow2(-nu0))  # 2^(-ν)
    m2Iepsp = _scale_series(Iepsp, -2)                                 # -2 i εp
    epsp_c, epsp_hi = _split_grade0(m2Iepsp)
    g1 = _gamma_full_series_shift(1 - s, epsp_c, epsp_hi, nmax)        # Γ[1-s-2iεp]
    g2 = _gamma_full_series(2 * nu0 + 2, _scale_series(delta, 2), nmax) # Γ[2ν+2]
    g3 = _gamma_full_series(nu0 + 1 - s, delta + Ieps, nmax)            # Γ[ν+1-s+iε]
    g4 = _gamma_full_series(2 * nu0 + 1, _scale_series(delta, 2), nmax) # Γ[1+2ν]
    g5 = _gamma_full_series(nu0 + 1 - s, delta - Ieps, nmax)            # Γ[1-s-iε+ν]
    g6 = _gamma_full_series_shift(nu0 + 1, -tau_c, delta - tau_hi, nmax) # Γ[1+ν-iτ]
    coeff = eps_pow
    coeff = _trunc_eta(coeff * Eieps, nmax)
    coeff = _trunc_eta(coeff * two_mnu, nmax)
    coeff = _trunc_eta(coeff * g1, nmax)
    coeff = _trunc_eta(coeff * g2, nmax)
    coeff = _trunc_eta(coeff * _amp_inv(g3), nmax)
    coeff = _trunc_eta(coeff * g4, nmax)
    coeff = _trunc_eta(coeff * _amp_inv(g5), nmax)
    coeff = _trunc_eta(coeff * _amp_inv(g6), nmax)
    # sums
    bU_a = nuser + _const_eta(1 + s, nmax) + Ieps        # 1+s+iε+ν
    bU_b = _scale_series(nuser, 2) + _const_eta(1, nmax) # 1+2ν
    bU_c = nuser + _const_eta(1, nmax) + Itau            # 1+ν+iτ
    bU_d = nuser + _const_eta(1 - s, nmax) - Ieps        # 1-s-iε+ν
    bU_e = nuser + _const_eta(1, nmax) - Itau            # 1+ν-iτ
    sumUp = _zero_eta(nmax)
    for n in 0:(cld(nmax, 3) + l + 2)
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && (n > 0 ? continue : nothing)
        num = _trunc_eta(_trunc_eta(_PH(bU_a, n, nmax) * _PH(bU_b, n, nmax), nmax) * _PH(bU_c, n, nmax), nmax)
        den = _trunc_eta(_PH(bU_d, n, nmax) * _PH(bU_e, n, nmax), nmax)
        t = _scale_series(_trunc_eta(_trunc_eta(num * _amp_inv(den), nmax) * an, nmax),
                          _signfac(n, n))
        _lead_eta(t) >= nmax && n > 0 && continue
        sumUp = sumUp + t
    end
    bD_a = nuser + _const_eta(1 + s, nmax) - Ieps        # ν+1+s-iε
    bD_b = nuser + _const_eta(1 - s, nmax) + Ieps        # ν+1-s+iε
    bD_c = _scale_series(nuser, 2) + _const_eta(2, nmax) # 2ν+2
    sumDown = _zero_eta(nmax)
    for n in -(cld(nmax, 3) + l + 2):0
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && n < 0 && continue
        ph = _trunc_eta(_PH(bD_a, n, nmax) * _amp_inv(_PH(bD_b, n, nmax)), nmax)
        # 1/((-n)! PH[2ν+2,n])
        fac = _signfac(n, -n)
        t = _scale_series(_trunc_eta(_trunc_eta(ph * _amp_inv(_PH(bD_c, n, nmax)), nmax) * an, nmax), fac)
        _lead_eta(t) >= nmax && n < 0 && continue
        sumDown = sumDown + t
    end
    sumUp = _expand_coeffs(sumUp); sumDown = _expand_coeffs(sumDown)
    _final_amp(_trunc_eta(_trunc_eta(coeff * sumUp, nmax) * _amp_inv(sumDown), nmax))
end

# K["Ratio"]  (PN.wl 2206-2242, r=0):
#   coeff = 2^(1+2ν) εκ (εκ)^(2ν) Γ[-2ν] Γ[1-s+iε+ν]/(Γ[-s+iε-ν] Γ[2+2ν])
#           · Γ[-1-2ν]Γ[1-s-iε+ν]Γ[1+ν-iτ]/(Γ[-s-iε-ν]Γ[1+2ν]Γ[-ν-iτ])
#   sumUpPH = Σ_n(-1)^(-n) PH[-1-2ν,-n]PH[s+iε-ν,-n]PH[-ν+iτ,-n]/((-n)!PH[-s-iε-ν,-n]PH[-ν-iτ,-n]) aMST[n]
#   sumUp2  = Σ_n(-1)^n/((-n)!PH[2ν+2,n]) PH[ν+1+s-iε,n]/PH[ν+1-s+iε,n] aMST[n]
#   sumDown = Σ_n(-1)^(-n) PH[s-iε-ν,-n]/((n)!PH[-2ν,-n]PH[-s+iε-ν,-n]) aMST[n]
#   sumDown2= Σ_n(-1)^n PH[1+s+iε+ν,n]PH[1+2ν,n]PH[1+ν+iτ,n]/(n!PH[1-s-iε+ν,n]PH[1+ν-iτ,n]) aMST[n]
# a==0: κ=1, τ=ε — the four spin args default to (1, ε, ε) byte-for-byte.
function _K_ratio(s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int,
        kap::Num = Num(1), tauser::PNSeries = eps, epspser::PNSeries = eps)
    iu = _IU()
    nu0 = l
    delta = nuser - _const_eta(nu0, nmax)
    Ieps = _scale_series(eps, iu)
    Itau = _scale_series(tauser, iu)                     # i τ
    tau_c, tau_hi = _split_grade0(Itau)
    # coeff
    two_2nup1 = _scale_series(_exp_series(_scale_series(delta, 2 * Num(log(2))), nmax), _pow2(1 + 2 * nu0)) # 2^(1+2ν)
    eps1 = _scale_series(eps, kap)                                      # εκ
    # (εκ)^(2ν) = (εκ)^(2l) exp(2δ log(εκ))  (ω-free)
    eps_2l = _epskap_intpow(2 * nu0, kap, nmax)
    logek = _log_epskap(kap, nmax)
    eps_2delta = _exp_series(_scale_series(delta, 2) * logek, nmax)
    eps_2nu = _trunc_eta(eps_2l * eps_2delta, nmax)
    gA = _gamma_full_series(-2 * nu0, _scale_series(delta, -2), nmax)   # Γ[-2ν]  (pole, K=-2l)
    gB = _gamma_full_series(nu0 + 1 - s, delta + Ieps, nmax)            # Γ[1-s+iε+ν]
    gC = _gamma_full_series(-nu0 - s, -delta + Ieps, nmax)              # Γ[-s+iε-ν]
    gD = _gamma_full_series(2 * nu0 + 2, _scale_series(delta, 2), nmax) # Γ[2+2ν]
    gE = _gamma_full_series(-2 * nu0 - 1, _scale_series(delta, -2), nmax) # Γ[-1-2ν] (pole)
    gF = _gamma_full_series(nu0 + 1 - s, delta - Ieps, nmax)            # Γ[1-s-iε+ν]
    gG = _gamma_full_series_shift(nu0 + 1, -tau_c, delta - tau_hi, nmax) # Γ[1+ν-iτ]
    gH = _gamma_full_series(-nu0 - s, -delta - Ieps, nmax)             # Γ[-s-iε-ν]
    gI = _gamma_full_series(2 * nu0 + 1, _scale_series(delta, 2), nmax) # Γ[1+2ν]
    gJ = _gamma_full_series_shift(-nu0, -tau_c, -delta - tau_hi, nmax) # Γ[-ν-iτ]
    coeff = two_2nup1
    coeff = _trunc_eta(coeff * eps1, nmax)
    coeff = _trunc_eta(coeff * eps_2nu, nmax)
    coeff = _trunc_eta(coeff * gA, nmax)
    coeff = _trunc_eta(coeff * gB, nmax)
    coeff = _trunc_eta(coeff * _amp_inv(gC), nmax)
    coeff = _trunc_eta(coeff * _amp_inv(gD), nmax)
    coeff = _trunc_eta(coeff * gE, nmax)
    coeff = _trunc_eta(coeff * gF, nmax)
    coeff = _trunc_eta(coeff * gG, nmax)
    coeff = _trunc_eta(coeff * _amp_inv(gH), nmax)
    coeff = _trunc_eta(coeff * _amp_inv(gI), nmax)
    coeff = _trunc_eta(coeff * _amp_inv(gJ), nmax)
    # sumUpPH (n from nMin to 0; uses PH[.,-n])
    pa = _scale_series(nuser, -2) + _const_eta(-1, nmax)   # -1-2ν
    pb = _const_eta(s, nmax) + Ieps - nuser                # s+iε-ν
    pc = _const_eta(0, nmax) + Itau - nuser                # -ν+iτ
    pd = _const_eta(-s, nmax) - Ieps - nuser               # -s-iε-ν
    pe = -nuser - Itau                                     # -ν-iτ
    sumUpPH = _zero_eta(nmax)
    for n in -(cld(nmax, 3) + l + 2):0
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && n < 0 && continue
        num = _trunc_eta(_trunc_eta(_PH(pa, -n, nmax) * _PH(pb, -n, nmax), nmax) * _PH(pc, -n, nmax), nmax)
        den = _trunc_eta(_PH(pd, -n, nmax) * _PH(pe, -n, nmax), nmax)
        fac = _signfac(-n, -n)
        t = _scale_series(_trunc_eta(_trunc_eta(num * _amp_inv(den), nmax) * an, nmax), fac)
        _lead_eta(t) >= nmax && n < 0 && continue
        sumUpPH = sumUpPH + t
    end
    # sumUp2 (n from nMin to 0)
    qa = nuser + _const_eta(1 + s, nmax) - Ieps           # ν+1+s-iε
    qb = nuser + _const_eta(1 - s, nmax) + Ieps           # ν+1-s+iε
    qc = _scale_series(nuser, 2) + _const_eta(2, nmax)    # 2ν+2
    sumUp2 = _zero_eta(nmax)
    for n in -(cld(nmax, 3) + l + 2):0
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && n < 0 && continue
        ph = _trunc_eta(_PH(qa, n, nmax) * _amp_inv(_PH(qb, n, nmax)), nmax)
        fac = _signfac(n, -n)
        t = _scale_series(_trunc_eta(_trunc_eta(ph * _amp_inv(_PH(qc, n, nmax)), nmax) * an, nmax), fac)
        _lead_eta(t) >= nmax && n < 0 && continue
        sumUp2 = sumUp2 + t
    end
    # sumDown (n from 0 to nMax; PH[.,-n])
    ra = _const_eta(s, nmax) - Ieps - nuser               # s-iε-ν
    rb = _scale_series(nuser, -2)                          # -2ν
    rc = _const_eta(-s, nmax) + Ieps - nuser              # -s+iε-ν
    sumDown = _zero_eta(nmax)
    for n in 0:(cld(nmax, 3) + l + 2)
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && n > 0 && continue
        den = _trunc_eta(_PH(rb, -n, nmax) * _PH(rc, -n, nmax), nmax)
        fac = _signfac(-n, n)
        t = _scale_series(_trunc_eta(_trunc_eta(_PH(ra, -n, nmax) * _amp_inv(den), nmax) * an, nmax), fac)
        _lead_eta(t) >= nmax && n > 0 && continue
        sumDown = sumDown + t
    end
    # sumDown2 (n from 0 to nMax) — same integrand as K["ν"] sumUp
    sa = nuser + _const_eta(1 + s, nmax) + Ieps
    sb = _scale_series(nuser, 2) + _const_eta(1, nmax)
    sc = nuser + _const_eta(1, nmax) + Itau               # 1+ν+iτ
    sd = nuser + _const_eta(1 - s, nmax) - Ieps
    se = nuser + _const_eta(1, nmax) - Itau               # 1+ν-iτ
    sumDown2 = _zero_eta(nmax)
    for n in 0:(cld(nmax, 3) + l + 2)
        an = get(aMSTeta, n, nothing); an === nothing && continue
        _lead_eta(an) >= nmax && n > 0 && continue
        num = _trunc_eta(_trunc_eta(_PH(sa, n, nmax) * _PH(sb, n, nmax), nmax) * _PH(sc, n, nmax), nmax)
        den = _trunc_eta(_PH(sd, n, nmax) * _PH(se, n, nmax), nmax)
        t = _scale_series(_trunc_eta(_trunc_eta(num * _amp_inv(den), nmax) * an, nmax),
                          _signfac(n, n))
        _lead_eta(t) >= nmax && n > 0 && continue
        sumDown2 = sumDown2 + t
    end
    # flat-expand the sums before inverting (inv reads coefficients; unexpanded 𝕚-trees are slow)
    sumUp2 = _expand_coeffs(sumUp2); sumDown = _expand_coeffs(sumDown)
    sumUpPH = _expand_coeffs(sumUpPH); sumDown2 = _expand_coeffs(sumDown2)
    r1 = _trunc_eta(sumUp2 * _amp_inv(sumDown), nmax)
    r2 = _trunc_eta(sumUpPH * _amp_inv(sumDown2), nmax)
    _final_amp(_trunc_eta(_trunc_eta(coeff * r1, nmax) * r2, nmax))
end

# Γ[K + x] for a CONSTANT-only (no-ν) argument K integer + x (x lead>=1): reuse helper.
_gamma_full_series_const(K::Int, x::PNSeries, nmax::Int) = _gamma_full_series(K, x, nmax)

# Γ[K + c + x] where `c` is a grade-0 COMPLEX shift (the τ/εp constant -ma/κ ·i etc.; ≠0 only
# at a≠0) and `x` has lead>=1.  The expansion base b = K + c is a complex NON-integer (the c
# shift moves K off any Γ pole), so Γ[b+x] = Γ[b] exp(Σ_{j>=1} d_j x^j), d_1 = ψ(b),
# d_j = (-1)^j ζ(j,b)/j  (Hurwitz ζ).  All numeric (ComplexF64) since a≠0 ⇒ irrational κ.
# At c==0 this falls back to the exact-rational integer-K helper (a==0 byte-for-byte).
function _gamma_full_series_shift(K::Int, c::Num, x::PNSeries, nmax::Int)
    cval = _cval(c)
    iszero(cval) && return _gamma_full_series(K, x, nmax)
    b = ComplexF64(K) + cval
    lead = _lead_eta(x)
    S = _zero_eta(nmax)
    xpow = _const_eta(1, nmax)
    j = 1
    while lead >= 1 && j * lead < nmax
        xpow = _trunc_eta(xpow * x, nmax)
        dj = j == 1 ? SpecialFunctions.digamma(b) :
             ((-1)^j / j) * SpecialFunctions.zeta(j, b)
        S = S + _scale_series(xpow, _cnum(dj))
        j += 1
    end
    gb = _cnum(SpecialFunctions.gamma(b))
    _scale_series(_exp_series(S, nmax), gb)
end

# (εκ)^p for integer p: εκ ω-free = 2κ at slot 3, so (εκ)^p = (2κ)^p η^(3p) (ω^p by grade).
# (2κ)^p is a clean complex const (κ=1 ⇒ 2^p exact; else _cnum to avoid a symbolic 1/𝕚 for p<0).
_epskap_intpow(p::Int, kap::Num, nmax::Int) =
    _mono_eta_signed(isequal(kap, Num(1)) ? _pow2(p) : _cnum(_cval(Num(2) * kap)^p), 3 * p, nmax)

# log(εκ) = log(ε) + log(κ) as an η-series CONSTANT (log ε = Log2+Logω symbolic; +log κ const).
_log_epskap(kap::Num, nmax::Int) = _log_eps(nmax) + _const_eta(_logkap(kap), nmax)

# ε^p = 2^p η^(3p) for integer p (ω-free; ω reattached by grade). p may be negative (Laurent).
_eps_intpow(p::Int, nmax::Int) = _mono_eta_signed(_pow2(p), 3 * p, nmax)

# η-monomial c η^q allowing NEGATIVE q (Laurent).  (the _mono_eta helper requires q>=0.)
function _mono_eta_signed(c, q::Int, nmax::Int)
    cN = Num(c)
    (iszero(cN) || q >= nmax) && return _zero_eta(nmax)
    PNSeries(:η, Num[cN], q, nmax, 1)
end

# log ε = Log[2 ω] = Log[2] + Log[ω]; the η^3-grade of ε is NOT an explicit log η (MMA's
# IgnoreExpansionParameter). ω-free: carried as the CONSTANT Log[2] + Logω (Logω symbolic).
_log_eps(nmax::Int) = _const_eta(Num(log(2)) + _logw(), nmax)

# ---------------------------------------------------------------------------
# B/C transmission + incidence amplitudes (PN.wl 1787-1848, 1930-1940).
# Spin shifts thread through κ (`kap`), τ (`tauser`), εp (`epspser`); at a==0 the a==0
# defaults (κ=1, τ=ε, εp=ε) recover the Schwarzschild amplitudes byte-for-byte.
# ---------------------------------------------------------------------------

# the a==0 spin defaults (κ=1, τ=ε, εp=ε) for the amplitude functions called WITHOUT spin args.
_spin0(eps::PNSeries) = (Num(1), eps, eps)

# Log[κ] as a Num PN-constant (κ a ComplexF64 PN-const). κ=1 ⇒ 0 (a==0).
_logkap(kap::Num) = isequal(kap, Num(1)) ? Num(0) : _cnum(log(_cval(kap)))

# exp of an η-series that MAY carry a grade-0 COMPLEX constant c0 (a≠0): exp(c0+hi) =
# exp(c0)·exp(hi), exp(c0) a scalar (numeric), hi lead>=1 handled by _exp_series. At a==0 the
# argument has lead>=1 (no constant) ⇒ plain _exp_series.
function _exp_const_series(x::PNSeries, nmax::Int)
    c0, hi = _split_grade0(x)
    iszero(c0) && return _exp_series(x, nmax)
    _scale_series(_exp_series(hi, nmax), _cnum(exp(_cval(c0))))
end
# the ComplexF64 value of a κ-like Num PN-constant (re + im·𝕚).
function _cval(c::Num)
    iu = _IU()
    re = Symbolics.value(Symbolics.substitute(c, Dict(iu => 0)))
    im = Symbolics.value(Symbolics.substitute(c - Num(re), Dict(iu => 1)))
    ComplexF64(re, im)
end

# C_trans (PN.wl 1930-1940): coeff = (ε/2)^(-1-2s) E^(i(ε Log ε - (1-κ)/2 ε)); a==0 ⇒ (1-κ)=0.
#   C_trans = coeff · A["-"].  `kap` is κ (a==0 ⇒ 1 ⇒ (1-κ)=0).
function _C_trans(s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int, kap::Num = Num(1))
    iu = _IU()
    halfeps_pow = _mono_eta_signed(Num(1), 3 * (-1 - 2s), nmax)   # (ε/2)^(-1-2s) = η^(3(-1-2s))
    # phase arg = ε Log ε - (1-κ)/2 ε; the -(1-κ)/2 ε term is a κ-const × ε (η^3, ω^1).
    phase = _trunc_eta(eps * _log_eps(nmax), nmax) -
            _scale_series(eps, (Num(1) - kap) / 2)
    Ephase = _exp_series(_scale_series(phase, iu), nmax)
    A = _A_amp(:minus, s, l, nuser, eps, aMSTeta, nmax)
    _final_amp(_trunc_eta(_trunc_eta(halfeps_pow * Ephase, nmax) * A, nmax))
end

# B_trans (PN.wl 1829-1848): coeff = (2κ)^(2s) E^(i κ εp (1+2 Log κ/(1+κ))); a==0 ⇒ κ=1,
#   εp=ε, Log κ=0 ⇒ coeff = 2^(2s) E^(i ε).  B_trans = coeff (Σ aMST)/K["ν"].
function _B_trans(s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int,
        kap::Num = Num(1), tauser::PNSeries = eps, epspser::PNSeries = eps)
    iu = _IU()
    # (2κ)^(2s): a const PN-factor (κ=1 ⇒ 2^(2s) real).  Computed as a clean complex const
    # (κ has 𝕚 ⇒ symbolic (2κ)^(2s) would form a 1/𝕚-fraction; _cnum(numeric) avoids it).
    coeff0 = isequal(kap, Num(1)) ? _pow2(2s) : _cnum(_cval(Num(2) * kap)^(2s))
    # E^(i κ εp (1 + 2 Log κ/(1+κ))): the scalar prefactor of εp (κ const), then exp the
    # εp η-series (which carries a grade-0 const at a≠0 ⇒ _exp_const_series).
    pref = isequal(kap, Num(1)) ? iu : _cnum(im * _cval(kap) * (1 + 2 * log(_cval(kap)) / (1 + _cval(kap))))
    expo = _scale_series(epspser, pref)
    coeff = _scale_series(_exp_const_series(expo, nmax), coeff0)
    ssum = _zero_eta(nmax)
    for (n, an) in aMSTeta
        _lead_eta(an) >= nmax && continue
        ssum = ssum + an
    end
    Knu = _K_nu(s, l, nuser, eps, aMSTeta, nmax, kap, tauser, epspser)
    _final_amp(_trunc_eta(_trunc_eta(coeff * ssum, nmax) * _amp_inv(Knu), nmax))
end

# B_inc (PN.wl 1787-1807): coeff = E^(-i(ε Log ε - (1-κ)/2 ε)) (ε/2)^(-1); a==0 ⇒ (1-κ)=0.
#   K1=1, K2coeff = -(i E^(-iπν) Sin[π(ν-s+iε)]/Sin[π(ν+s-iε)]) (ε ONLY, no τ),
#   K2 = K["Ratio"].  B_inc = A["+"] coeff (K1 + K2coeff K2).
function _B_inc(s::Int, l::Int, nuser::PNSeries, eps::PNSeries,
        aMSTeta::Dict{Int,PNSeries}, nmax::Int,
        kap::Num = Num(1), tauser::PNSeries = eps, epspser::PNSeries = eps)
    iu = _IU()
    nu0 = l
    delta = nuser - _const_eta(nu0, nmax)
    Ieps = _scale_series(eps, iu)
    halfeps_inv = _mono_eta_signed(Num(1), -3, nmax)              # (ε/2)^(-1) = η^-3 (ω-free)
    phase = _trunc_eta(eps * _log_eps(nmax), nmax) -
            _scale_series(eps, (Num(1) - kap) / 2)                # ε Log ε - (1-κ)/2 ε
    Ephase = _exp_series(_scale_series(phase, -iu), nmax)         # E^(-i(...))
    coeff = _trunc_eta(halfeps_inv * Ephase, nmax)
    # K2coeff = -(i E^(-iπν) Sin[π(ν-s+iε)]/Sin[π(ν+s-iε)])
    Eipinu = _scale_series(_exp_series(_scale_series(delta, -iu * _PISYM), nmax), Num((-1)^nu0))  # E^(-iπν)
    sinA = _sin_pi_series(nu0 - s, delta + Ieps, nmax)            # Sin[π(ν-s+iε)]
    sinB = _sin_pi_series(nu0 + s, delta - Ieps, nmax)            # Sin[π(ν+s-iε)]
    K2coeff = _scale_series(_trunc_eta(_trunc_eta(Eipinu * sinA, nmax) * _amp_inv(sinB), nmax), -iu)
    K2 = _K_ratio(s, l, nuser, eps, aMSTeta, nmax, kap, tauser, epspser)
    K1plus = _const_eta(1, nmax) + _trunc_eta(K2coeff * K2, nmax)
    A = _A_amp(:plus, s, l, nuser, eps, aMSTeta, nmax)
    _final_amp(_trunc_eta(_trunc_eta(A * coeff, nmax) * K1plus, nmax))
end


end # module RadialMode
