# Regularized l-sum (tail fit) → value ± uncertainty. Port of Punctures/regsum.jl (RegSum/NoRegSum,
# hR.m L50-60, Wardell et al. arXiv:2508.00087), REDESIGNED:
#
#   • the analytic ∞-tail digamma `ψ` is COLLAPSED off SpecialFunctions.digamma onto
#     `BigSpecial.digamma_big` (AHA — one digamma primitive, and the only one that carries
#     Complex{BigFloat}; here used on a real argument via real∘digamma_big∘complex). This makes the
#     uncertainty tail type-generic so the BigFloat l-tables inherit precision.
#
# Least-squares fit of the last five l-values to the fixed basis {b1,b2}, subtract the fit from every
# l, total the remainder, and report the analytic ∞-tail of the last fit term as the uncertainty.

import BigSpecial

# Digamma on a real argument, type-generic (BigFloat-ready): real part of the complex digamma_big.
_psi_real(x) = real(BigSpecial.digamma_big(complex(x)))

"""
    RegResult(value, uncertainty)

A regularized l-sum result: a central `value` with an `uncertainty` estimate. Mirrors the
Mathematica `Around[value, unc]` produced by `RegSum`/`NoRegSum` (hR.m L50-60). A 2-field struct
(not Measurements.jl) since only construction and the two fields are used downstream.
"""
struct RegResult{T<:Number}
    value::T
    uncertainty::T
end

# Fixed tail-fit basis (hR.m L52): b1, b2 as functions of l.
_regsum_b1(l) = 1 / ((2l - 1) * (2l + 3))
_regsum_b2(l) = 1 / ((2l - 3) * (2l - 1) * (2l + 3) * (2l + 5))

"""
    _regsum_tail(L)

Analytic remainder `Σ_{l=L}^∞ 1/((2l-3)(2l-1)(2l+3)(2l+5))`, the infinite tail of the last fit
basis function (hR.m L54). Partial fractions in `2l` give coefficients `{1/96, -1/48, 1/48, -1/96}`
on `{2l-3, 2l-1, 2l+3, 2l+5}` (summing to 0, so the per-term digamma divergences cancel), hence
`-(1/2) Σ a_k ψ(L + c_k/2)` with `ψ = BigSpecial.digamma_big`.
"""
function _regsum_tail(L)
    T = float(typeof(L))
    half = one(T) / 2
    return -half * (
        (one(T) / 96) * _psi_real(L - 3 * half) -
        (one(T) / 48) * _psi_real(L - half) +
        (one(T) / 48) * _psi_real(L + 3 * half) -
        (one(T) / 96) * _psi_real(L + 5 * half))
end

"""
    reg_sum(field) -> RegResult

Regularized l-sum of `field` (`field[l+1]` = value at multipole `l`, l = 0..lmax). Ports `RegSum`
(hR.m L50-54): least-squares fit of the last five l-values to the fixed basis {b1, b2}, subtract the
fit from every l, total the remainder, and report the analytic ∞-tail of the last fit term
(`c2 · Σ_{l=lmax+1}^∞ b2`) as the uncertainty.
"""
function reg_sum(field::AbstractVector{T}) where {T<:Number}
    lmax = length(field) - 1
    lmax >= 4 || throw(ArgumentError("reg_sum needs at least 5 l-values (lmax >= 4)"))
    R = float(real(T))
    ls = (lmax-4):lmax
    B = R[k == 1 ? _regsum_b1(l) : _regsum_b2(l) for l in ls, k in 1:2]
    y = collect(field[l+1] for l in ls)
    c = B \ y
    fitted = [c[1] * _regsum_b1(l) + c[2] * _regsum_b2(l) for l in 0:lmax]
    total = sum(field[l+1] - fitted[l+1] for l in 0:lmax)
    unc = (c[1] == 0 && c[2] == 0) ? zero(real(T)) :
          abs(c[2] * _regsum_tail(R(lmax + 1)))
    return RegResult(total, oftype(total, unc))
end

"""
    noreg_sum(field) -> RegResult

Plain l-sum of `field` with the last term as the uncertainty estimate. Ports `NoRegSum`
(hR.m L57-60): no tail fit; `value = Σ field`, `unc = |field[lmax]|`.
"""
function noreg_sum(field::AbstractVector{T}) where {T<:Number}
    lmax = length(field) - 1
    lmax >= 0 || throw(ArgumentError("noreg_sum needs at least one l-value"))
    return RegResult(sum(field), oftype(sum(field), abs(field[end])))
end
